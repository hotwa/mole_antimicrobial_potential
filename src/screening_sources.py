"""Source adapters yielding normalized chunks for batch screening."""

from __future__ import annotations

import fnmatch
import sqlite3
import tarfile
from io import TextIOWrapper
from pathlib import Path
from typing import Iterator

import pandas as pd
import pyarrow.parquet as pq

from src.screening_planner import WorkUnit

DEFAULT_ARCHIVE_PATTERN = "*_scheme_b_unique_products.csv"
DEFAULT_ARCHIVE_SMILES_COL = "product_smiles_canonical"
DEFAULT_ARCHIVE_CHEM_ID_COL = "example_combo_id"
DEFAULT_INPUT_CHUNK_SIZE = 10000


def _basename_without_suffixes(path: Path) -> str:
    name = path.name
    for suffix in reversed(path.suffixes):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


def _infer_tabular_separator(path: Path) -> str:
    suffixes = [suffix.lower() for suffix in path.suffixes]
    if suffixes[-2:] == [".tsv", ".gz"] or path.suffix.lower() in {".tsv", ".txt"}:
        return "\t"
    return ","


def _uniquify_chem_ids(frame: pd.DataFrame, seen_ids: set[str]) -> pd.DataFrame:
    # We mutate in place or copy? Better copy if deduplicating across chunks.
    # Actually, deduplication of chem_ids globally across chunks might be expensive
    # but we will just ensure they don't overlap by renaming.
    counts: dict[str, int] = {}
    unique_ids: list[str] = []
    for raw_id in frame["chem_id"].astype(str):
        if raw_id not in seen_ids and raw_id not in counts:
            unique_ids.append(raw_id)
            seen_ids.add(raw_id)
            counts[raw_id] = 1
        else:
            cand = f"{raw_id}__1"
            idx = 1
            while cand in seen_ids or cand in counts:
                idx += 1
                cand = f"{raw_id}__{idx}"
            unique_ids.append(cand)
            seen_ids.add(cand)
            counts[cand] = 1

    deduped = frame.copy()
    deduped["chem_id"] = unique_ids
    return deduped


def _load_tabular_input(
    path: Path,
    smiles_colname: str,
    chem_id_colname: str,
    chunk_size: int,
) -> Iterator[pd.DataFrame]:
    reader = pd.read_csv(path, sep=_infer_tabular_separator(path), chunksize=chunk_size)
    try:
        global_row_idx = 1
        for frame in reader:
            if smiles_colname not in frame.columns:
                raise ValueError(f"Missing SMILES column '{smiles_colname}' in {path}")

            if chem_id_colname not in frame.columns:
                frame[chem_id_colname] = [f"mol{global_row_idx + i}" for i in range(len(frame))]

            frame = frame.rename(columns={smiles_colname: "smiles", chem_id_colname: "chem_id"})
            frame["source_file"] = str(path)
            frame["source_group"] = _basename_without_suffixes(path)
            frame["source_row"] = range(global_row_idx, global_row_idx + len(frame))
            global_row_idx += len(frame)
            yield frame
    finally:
        reader.close()


def _infer_archive_source_group(member_name: str) -> str:
    parts = Path(member_name).parts
    if len(parts) >= 2:
        return parts[1]
    return Path(member_name).stem


def _load_archive_input(
    path: Path,
    archive_pattern: str,
    archive_smiles_colname: str,
    archive_chem_id_colname: str,
    fallback_smiles_colname: str,
    fallback_chem_id_colname: str,
    chunk_size: int,
) -> Iterator[pd.DataFrame]:
    with tarfile.open(path, "r:*") as archive:
        members = [
            member
            for member in archive.getmembers()
            if member.isfile() and fnmatch.fnmatch(Path(member.name).name, archive_pattern)
        ]
        if not members:
            raise ValueError(f"No files matched pattern '{archive_pattern}' inside {path}")

        for member in sorted(members, key=lambda item: item.name):
            extracted = archive.extractfile(member)
            if extracted is None:
                continue

            # Since tar members are usually small, read the whole member
            frame = pd.read_csv(TextIOWrapper(extracted, encoding="utf-8"))
            smiles_col = archive_smiles_colname if archive_smiles_colname in frame.columns else fallback_smiles_colname
            chem_id_col = archive_chem_id_colname if archive_chem_id_colname in frame.columns else fallback_chem_id_colname

            if smiles_col not in frame.columns:
                raise ValueError(f"Missing SMILES column in archive member {member.name}")

            if chem_id_col not in frame.columns:
                member_stem = _basename_without_suffixes(Path(member.name))
                frame[chem_id_col] = [f"{member_stem}__{i + 1}" for i in range(len(frame))]

            frame = frame.rename(columns={smiles_col: "smiles", chem_id_col: "chem_id"})
            frame["source_file"] = member.name
            frame["source_group"] = _infer_archive_source_group(member.name)
            frame["source_row"] = range(1, len(frame) + 1)
            frame["source_member"] = member.name

            # Yield chunks if the member is large
            for start in range(0, len(frame), chunk_size):
                yield frame.iloc[start : start + chunk_size].copy()


def _load_sqlite_input(
    path: Path,
    smiles_colname: str,
    chem_id_colname: str,
    sqlite_table: str | None,
    sqlite_query: str | None,
    chunk_size: int,
) -> Iterator[pd.DataFrame]:
    with sqlite3.connect(path) as connection:
        if sqlite_query:
            query = sqlite_query
            source_group = sqlite_table or "query"
        else:
            if sqlite_table:
                source_group = sqlite_table
            else:
                tables = pd.read_sql_query(
                    "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'",
                    connection,
                )["name"].tolist()
                candidates: list[str] = []
                for table_name in tables:
                    columns = pd.read_sql_query(f'PRAGMA table_info("{table_name}")', connection)["name"].tolist()
                    if smiles_colname in columns:
                        candidates.append(table_name)
                if not candidates:
                    raise ValueError(
                        f"No SQLite table in {path} contains the required SMILES column '{smiles_colname}'"
                    )
                if len(candidates) > 1:
                    raise ValueError(
                        f"Multiple SQLite tables in {path} contain SMILES column '{smiles_colname}': "
                        f"{candidates}. Pass --sqlite-table or --sqlite-query."
                    )
                source_group = candidates[0]
            query = f'SELECT * FROM "{source_group}"'

        global_row_idx = 1
        for frame in pd.read_sql_query(query, connection, chunksize=chunk_size):
            if smiles_colname not in frame.columns:
                raise ValueError(f"Missing SMILES column '{smiles_colname}' in {path}")

            if chem_id_colname not in frame.columns:
                frame[chem_id_colname] = [f"{source_group}__{global_row_idx + i}" for i in range(len(frame))]

            frame = frame.rename(columns={smiles_colname: "smiles", chem_id_colname: "chem_id"})
            frame["source_file"] = str(path)
            frame["source_group"] = str(source_group)
            frame["source_row"] = range(global_row_idx, global_row_idx + len(frame))
            frame["source_member"] = str(source_group)
            global_row_idx += len(frame)
            yield frame


def _load_parquet_input(
    path: Path,
    smiles_colname: str,
    chem_id_colname: str,
    chunk_size: int,
) -> Iterator[pd.DataFrame]:
    parquet_file = pq.ParquetFile(path)
    schema_names = set(parquet_file.schema.names)
    smiles_column = "smiles" if "smiles" in schema_names else smiles_colname
    chem_id_column = "chem_id" if "chem_id" in schema_names else chem_id_colname
    source_group_column = "source_group" if "source_group" in schema_names else None

    if smiles_column not in schema_names:
        raise ValueError(f"Missing SMILES column '{smiles_colname}' in {path}")

    columns = [smiles_column]
    if chem_id_column in schema_names:
        columns.append(chem_id_column)
    if source_group_column:
        columns.append(source_group_column)

    row_offset = 0
    default_source_group = _basename_without_suffixes(path)
    for batch in parquet_file.iter_batches(batch_size=chunk_size, columns=columns):
        frame = batch.to_pandas()
        if chem_id_column not in frame.columns:
            frame[chem_id_column] = [f"{default_source_group}__{row_offset + 1 + i}" for i in range(len(frame))]

        frame = frame.rename(columns={smiles_column: "smiles", chem_id_column: "chem_id"})
        frame["source_file"] = str(path)
        if source_group_column and source_group_column in frame.columns:
            frame["source_group"] = frame[source_group_column].astype(str)
        else:
            frame["source_group"] = default_source_group
        frame["source_row"] = range(row_offset + 1, row_offset + 1 + len(frame))
        row_offset += len(frame)
        yield frame


def stream_screening_input(
    input_path: str | Path,
    smiles_colname: str = "smiles",
    chem_id_colname: str = "chem_id",
    archive_pattern: str = DEFAULT_ARCHIVE_PATTERN,
    archive_smiles_colname: str = DEFAULT_ARCHIVE_SMILES_COL,
    archive_chem_id_colname: str = DEFAULT_ARCHIVE_CHEM_ID_COL,
    sqlite_table: str | None = None,
    sqlite_query: str | None = None,
    dedupe_smiles: bool = True,
    chunk_size: int = DEFAULT_INPUT_CHUNK_SIZE,
) -> Iterator[pd.DataFrame]:
    """Iterate over batches of an input source."""

    path = Path(input_path).expanduser().resolve()

    if path.is_file() and tarfile.is_tarfile(path):
        iterator = _load_archive_input(
            path=path,
            archive_pattern=archive_pattern,
            archive_smiles_colname=archive_smiles_colname,
            archive_chem_id_colname=archive_chem_id_colname,
            fallback_smiles_colname=smiles_colname,
            fallback_chem_id_colname=chem_id_colname,
            chunk_size=chunk_size,
        )
    elif path.is_file() and path.suffix.lower() in {".sqlite", ".sqlite3", ".db", ".db3"}:
        iterator = _load_sqlite_input(
            path=path,
            smiles_colname=smiles_colname,
            chem_id_colname=chem_id_colname,
            sqlite_table=sqlite_table,
            sqlite_query=sqlite_query,
            chunk_size=chunk_size,
        )
    elif path.is_file() and path.suffix.lower() == ".parquet":
        iterator = _load_parquet_input(
            path=path,
            smiles_colname=smiles_colname,
            chem_id_colname=chem_id_colname,
            chunk_size=chunk_size,
        )
    else:
        iterator = _load_tabular_input(
            path=path,
            smiles_colname=smiles_colname,
            chem_id_colname=chem_id_colname,
            chunk_size=chunk_size,
        )

    seen_smiles: set[str] = set()
    seen_ids: set[str] = set()
    global_order = 0

    for frame in iterator:
        frame = frame.dropna(subset=["smiles"]).copy()
        frame["smiles"] = frame["smiles"].astype(str)
        frame["chem_id"] = frame["chem_id"].astype(str)

        if dedupe_smiles:
            # Drop duplicates internally
            frame = frame.drop_duplicates(subset=["smiles"], keep="first")
            # Drop duplicates against seen chunks
            mask = ~frame["smiles"].isin(seen_smiles)
            frame = frame[mask].copy()
            seen_smiles.update(frame["smiles"])

        if frame.empty:
            continue

        frame = _uniquify_chem_ids(frame, seen_ids)
        frame["input_order"] = range(global_order, global_order + len(frame))
        global_order += len(frame)

        yield frame.reset_index(drop=True)


def process_work_unit(
    unit: WorkUnit,
    smiles_colname: str,
    chem_id_colname: str,
    archive_smiles_colname: str,
    archive_chem_id_colname: str,
    chunk_size: int = DEFAULT_INPUT_CHUNK_SIZE,
) -> Iterator[pd.DataFrame]:
    """Parse a single WorkUnit into a normalized DataFrame iterator."""


    if unit.source_type == "archive":
        # extract just this member from tar
        with tarfile.open(unit.source_path, "r:*") as archive:
            member = archive.getmember(unit.source_member)
            extracted = archive.extractfile(member)
            if extracted:
                row_offset = unit.start_row or 0
                for frame in pd.read_csv(TextIOWrapper(extracted, encoding="utf-8"), chunksize=chunk_size):
                    smiles_col = archive_smiles_colname if archive_smiles_colname in frame.columns else smiles_colname
                    chem_id_col = archive_chem_id_colname if archive_chem_id_colname in frame.columns else chem_id_colname

                    if smiles_col not in frame.columns:
                        raise ValueError(f"Missing SMILES column in archive member {unit.source_member}")

                    if chem_id_col not in frame.columns:
                        member_stem = _basename_without_suffixes(Path(unit.source_member))
                        frame[chem_id_col] = [f"{member_stem}__{row_offset + i + 1}" for i in range(len(frame))]

                    frame = frame.rename(columns={smiles_col: "smiles", chem_id_col: "chem_id"})
                    frame["source_file"] = unit.source_member
                    frame["source_group"] = unit.source_group or unit.group_id
                    frame["source_member"] = unit.source_member

                    if unit.start_row is not None and unit.max_rows is not None:
                        # Since we do not use skiprows for tar, slice manually or handle
                        # Assuming start_row/max_rows isn't used for archive, we just assign source_row
                        pass

                    frame["source_row"] = range(row_offset + 1, row_offset + 1 + len(frame))
                    row_offset += len(frame)

                    frame = frame.dropna(subset=["smiles"]).copy()
                    if not frame.empty:
                        frame["smiles"] = frame["smiles"].astype(str)
                        frame["chem_id"] = frame["chem_id"].astype(str)
                        yield frame

    elif unit.source_type == "parquet":
        parquet_file = pq.ParquetFile(unit.source_path)
        schema_names = set(parquet_file.schema.names)
        smiles_column = "smiles" if "smiles" in schema_names else smiles_colname
        chem_id_column = "chem_id" if "chem_id" in schema_names else chem_id_colname
        source_group_column = "source_group" if "source_group" in schema_names else None

        if smiles_column not in schema_names:
            raise ValueError(f"Missing SMILES column '{smiles_colname}' in {unit.source_path}")

        columns = [smiles_column]
        if chem_id_column in schema_names:
            columns.append(chem_id_column)
        if source_group_column:
            columns.append(source_group_column)

        row_offset = 0
        for batch in parquet_file.iter_batches(batch_size=chunk_size, columns=columns):
            frame = batch.to_pandas()
            if chem_id_column not in frame.columns:
                frame[chem_id_column] = [f"{Path(unit.source_path).stem}__{row_offset + 1 + i}" for i in range(len(frame))]

            frame = frame.rename(columns={smiles_column: "smiles", chem_id_column: "chem_id"})
            frame["source_file"] = str(unit.source_path)
            if source_group_column and source_group_column in frame.columns:
                frame["source_group"] = frame[source_group_column].astype(str)
            else:
                frame["source_group"] = unit.source_group or unit.group_id
            if unit.source_member:
                frame["source_member"] = str(unit.source_member)
            frame["source_row"] = range(row_offset + 1, row_offset + 1 + len(frame))
            row_offset += len(frame)

            frame = frame.dropna(subset=["smiles"]).copy()
            if not frame.empty:
                frame["smiles"] = frame["smiles"].astype(str)
                frame["chem_id"] = frame["chem_id"].astype(str)
                yield frame

    elif unit.source_type == "sqlite":
        with sqlite3.connect(unit.source_path) as conn:
            query = f'SELECT * FROM "{unit.source_member}"'
            if unit.start_row is not None and unit.max_rows is not None:
                query += f" LIMIT {unit.max_rows} OFFSET {unit.start_row}"

            row_offset = unit.start_row or 0
            for frame in pd.read_sql_query(query, conn, chunksize=chunk_size):
                if smiles_colname not in frame.columns:
                    raise ValueError(f"Missing SMILES column '{smiles_colname}' in {unit.source_path}")

                if chem_id_colname not in frame.columns:
                    frame[chem_id_colname] = [f"{unit.source_member}__{row_offset + 1 + i}" for i in range(len(frame))]

                frame = frame.rename(columns={smiles_colname: "smiles", chem_id_colname: "chem_id"})
                frame["source_file"] = str(unit.source_path)
                frame["source_group"] = unit.source_group or unit.group_id
                frame["source_member"] = unit.source_member
                frame["source_row"] = range(row_offset + 1, row_offset + 1 + len(frame))
                row_offset += len(frame)

                frame = frame.dropna(subset=["smiles"]).copy()
                if not frame.empty:
                    frame["smiles"] = frame["smiles"].astype(str)
                    frame["chem_id"] = frame["chem_id"].astype(str)
                    yield frame

    elif unit.source_type == "tabular":
        kwargs = {"sep": _infer_tabular_separator(Path(unit.source_path)), "chunksize": chunk_size}

        if unit.start_row is not None and unit.max_rows is not None and unit.start_row > 0:
            header_frame = pd.read_csv(unit.source_path, nrows=0, sep=kwargs["sep"])
            kwargs["skiprows"] = unit.start_row + 1
            kwargs["nrows"] = unit.max_rows
            kwargs["names"] = header_frame.columns
            kwargs["header"] = None
        elif unit.max_rows is not None:
            kwargs["nrows"] = unit.max_rows

        row_offset = unit.start_row or 0
        reader = pd.read_csv(unit.source_path, **kwargs)
        try:
            for frame in reader:
                if smiles_colname not in frame.columns:
                    raise ValueError(f"Missing SMILES column '{smiles_colname}' in {unit.source_path}")

                if chem_id_colname not in frame.columns:
                    frame[chem_id_colname] = [f"mol{row_offset + 1 + i}" for i in range(len(frame))]

                frame = frame.rename(columns={smiles_colname: "smiles", chem_id_colname: "chem_id"})
                frame["source_file"] = str(unit.source_path)
                frame["source_group"] = unit.source_group or unit.group_id
                frame["source_row"] = range(row_offset + 1, row_offset + 1 + len(frame))
                row_offset += len(frame)

                frame = frame.dropna(subset=["smiles"]).copy()
                if not frame.empty:
                    frame["smiles"] = frame["smiles"].astype(str)
                    frame["chem_id"] = frame["chem_id"].astype(str)
                    yield frame
        finally:
            reader.close()
