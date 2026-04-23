from __future__ import annotations

import math
import os
import sqlite3
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Dict, Any

import pyarrow.parquet as pq

@dataclass
class WorkUnit:
    source_type: str
    source_path: str
    source_member: Optional[str] = None
    group_id: str = ""
    source_group: str = ""
    group_index: int = 0
    estimated_rows: int = 0
    estimated_bytes: int = 0
    start_row: Optional[int] = None
    max_rows: Optional[int] = None

@dataclass
class PlannerConfig:
    grouping_mode: str = "auto"
    cpu_workers: int | str = "auto"
    prefetch_queue_size: int | str = "auto"
    target_rows_per_group: int | str = "auto"
    target_bytes_per_group: int | str = "auto"
    chunk_size: int = 10000

class ScreeningPlanner:
    def __init__(self, config: PlannerConfig = PlannerConfig()):
        self.config = config

        self.grouping_mode = self.config.grouping_mode
        if self.grouping_mode not in ("auto", "source", "chunk", "none"):
            self.grouping_mode = "auto"

        self.cpu_workers = self._parse_auto(self.config.cpu_workers, default=max(2, min(os.cpu_count() or 4, 8)))
        if self.config.cpu_workers == "auto":
            # auto rule for cpu workers: max(2, physical_cores // 2), capped at 8
            cores = os.cpu_count() or 4
            self.cpu_workers = min(max(2, cores // 2), 8)

        self.prefetch_queue_size = self._parse_auto(self.config.prefetch_queue_size, default=self.cpu_workers * 2)
        self.target_rows_per_group = self._parse_auto(self.config.target_rows_per_group, default=100000)
        self.target_bytes_per_group = self._parse_auto(self.config.target_bytes_per_group, default=50 * 1024 * 1024)

    def _parse_auto(self, val: int | str, default: int) -> int:
        if str(val).lower() == "auto":
            return default
        return int(val)

    def plan(
        self,
        input_path: str | Path,
        archive_pattern: str,
        sqlite_table: Optional[str] = None,
        sqlite_query: Optional[str] = None,
        smiles_colname: str = "smiles",
    ) -> List[WorkUnit]:
        path = Path(input_path).expanduser().resolve()

        if path.is_file() and tarfile.is_tarfile(path):
            return self._plan_archive(path, archive_pattern)
        elif path.is_file() and path.suffix.lower() == ".parquet":
            return [self._build_parquet_unit(path, default_source_group=path.stem)]
        elif path.is_dir():
            parquet_files = sorted(candidate for candidate in path.rglob("*.parquet") if candidate.is_file())
            if parquet_files:
                return [
                    self._build_parquet_unit(candidate, default_source_group=candidate.parent.name or candidate.stem)
                    for candidate in parquet_files
                ]
            return self._plan_tabular(path)
        elif path.is_file() and path.suffix.lower() in {".sqlite", ".sqlite3", ".db", ".db3"}:
            return self._plan_sqlite(path, sqlite_table, sqlite_query, smiles_colname)
        else:
            return self._plan_tabular(path)

    def _plan_archive(self, path: Path, archive_pattern: str) -> List[WorkUnit]:
        units = []
        import fnmatch
        with tarfile.open(path, "r:*") as archive:
            members = [
                m for m in archive.getmembers()
                if m.isfile() and fnmatch.fnmatch(Path(m.name).name, archive_pattern)
            ]

            for m in sorted(members, key=lambda x: x.name):
                # approximate rows based on file size
                estimated_rows = max(1, m.size // 100)

                group_id_base = self._infer_archive_source_group(m.name)

                # NO row-split for archive to avoid fake split reading whole files multiple times
                gid = group_id_base if self.grouping_mode != "none" else "all"
                units.append(WorkUnit(
                    source_type="archive",
                    source_path=str(path),
                    source_member=m.name,
                    group_id=gid,
                    source_group=group_id_base,
                    group_index=0,
                    estimated_rows=estimated_rows,
                    estimated_bytes=m.size
                ))
        return units

    def _plan_sqlite(self, path: Path, sqlite_table: Optional[str], sqlite_query: Optional[str], smiles_colname: str) -> List[WorkUnit]:
        units = []
        with sqlite3.connect(path) as conn:
            if sqlite_query:
                # hard to estimate rows for arbitrary query without running it as count
                import pandas as pd
                try:
                    count = pd.read_sql_query(f"SELECT count(*) as c FROM ({sqlite_query})", conn).iloc[0]["c"]
                except Exception:
                    count = 1000  # fallback
                table = sqlite_table or "query"
                units.extend(self._split_unit("sqlite", str(path), table, table, count, count * 100))
            else:
                tables = []
                if sqlite_table:
                    tables = [sqlite_table]
                else:
                    import pandas as pd
                    all_t = pd.read_sql_query(
                        "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'",
                        conn,
                    )["name"].tolist()
                    for t in all_t:
                        cols = pd.read_sql_query(f'PRAGMA table_info("{t}")', conn)["name"].tolist()
                        if smiles_colname in cols:
                            tables.append(t)

                for t in tables:
                    import pandas as pd
                    count = pd.read_sql_query(f"SELECT count(*) as c FROM \"{t}\"", conn).iloc[0]["c"]

                    if self.grouping_mode in ("auto", "chunk") and count > self.target_rows_per_group:
                        units.extend(self._split_unit("sqlite", str(path), t, t, count, count * 100))
                    else:
                        gid = t if self.grouping_mode != "none" else "all"
                        units.append(WorkUnit(
                            source_type="sqlite",
                            source_path=str(path),
                            source_member=t,
                            group_id=gid,
                            source_group=t,
                            group_index=0,
                            estimated_rows=count,
                            estimated_bytes=count * 100
                        ))
        return units

    def _plan_tabular(self, path: Path) -> List[WorkUnit]:
        size = path.stat().st_size
        estimated_rows = self._count_tabular_rows(path)

        group_id_base = self._basename_without_suffixes(path)

        if self.grouping_mode in ("auto", "chunk") and (
            estimated_rows > self.target_rows_per_group or
            size > self.target_bytes_per_group
        ):
            return self._split_unit(
                source_type="tabular",
                source_path=str(path),
                source_member=None,
                group_id_base=group_id_base,
                total_rows=estimated_rows,
                total_bytes=size
            )
        else:
            gid = group_id_base if self.grouping_mode != "none" else "all"
            return [WorkUnit(
                source_type="tabular",
                source_path=str(path),
                source_member=None,
                group_id=gid,
                source_group=group_id_base,
                group_index=0,
                estimated_rows=estimated_rows,
                estimated_bytes=size
            )]

    def _count_tabular_rows(self, path: Path) -> int:
        """Count data rows in a tabular file, excluding the header line.

        The planner uses these row counts as authoritative slice boundaries for
        ``skiprows``/``nrows`` in ``process_work_unit()``, so rough byte-based
        estimates are not safe here.
        """
        newline_count = 0
        has_bytes = False
        last_byte = b""

        with path.open("rb") as handle:
            while True:
                chunk = handle.read(1024 * 1024)
                if not chunk:
                    break
                has_bytes = True
                newline_count += chunk.count(b"\n")
                last_byte = chunk[-1:]

        if not has_bytes:
            return 0

        total_lines = newline_count
        if last_byte not in {b"\n", b"\r"}:
            total_lines += 1

        # Screening tables are expected to have a single header row.
        return max(0, total_lines - 1)

    def _split_unit(self, source_type: str, source_path: str, source_member: Optional[str], group_id_base: str, total_rows: int, total_bytes: int) -> List[WorkUnit]:
        units = []

        if self.grouping_mode == "none":
            group_id_base = "all"

        chunks_by_rows = math.ceil(total_rows / max(1, self.target_rows_per_group))
        chunks_by_bytes = math.ceil(total_bytes / max(1, self.target_bytes_per_group))
        num_chunks = max(1, chunks_by_rows, chunks_by_bytes)

        rows_per_chunk = math.ceil(total_rows / num_chunks)
        bytes_per_chunk = math.ceil(total_bytes / num_chunks)

        for i in range(num_chunks):
            start = i * rows_per_chunk
            count = min(rows_per_chunk, total_rows - start)
            if count <= 0:
                continue

            gid = f"{group_id_base}_part{i+1}" if self.grouping_mode != "none" else "all"

            units.append(WorkUnit(
                source_type=source_type,
                source_path=source_path,
                source_member=source_member,
                group_id=gid,
                source_group=group_id_base if self.grouping_mode != "none" else "all",
                group_index=i,
                estimated_rows=count,
                estimated_bytes=bytes_per_chunk,
                start_row=start,
                max_rows=count
            ))

        return units

    def _infer_archive_source_group(self, member_name: str) -> str:
        parts = Path(member_name).parts
        if len(parts) >= 2:
            return parts[1]
        return Path(member_name).stem

    def _build_parquet_unit(self, path: Path, default_source_group: str) -> WorkUnit:
        parquet_file = pq.ParquetFile(path)
        estimated_rows = parquet_file.metadata.num_rows
        source_group = self._infer_parquet_source_group(
            parquet_file, default_source_group=default_source_group
        )
        group_id = path.stem if self.grouping_mode != "none" else "all"
        return WorkUnit(
            source_type="parquet",
            source_path=str(path.resolve()),
            group_id=group_id,
            source_group=source_group,
            group_index=0,
            estimated_rows=estimated_rows,
            estimated_bytes=path.stat().st_size,
        )

    def _infer_parquet_source_group(self, parquet_file: pq.ParquetFile, default_source_group: str) -> str:
        if "source_group" in parquet_file.schema.names and parquet_file.metadata.num_rows > 0:
            first_row_group = parquet_file.read_row_group(0, columns=["source_group"])
            if first_row_group.num_rows > 0:
                value = first_row_group.column("source_group")[0].as_py()
                if value not in (None, ""):
                    return str(value)
        return default_source_group

    def _basename_without_suffixes(self, path: Path) -> str:
        name = path.name
        for suffix in reversed(path.suffixes):
            if name.endswith(suffix):
                name = name[: -len(suffix)]
        return name
