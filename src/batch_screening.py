"""Batch screening helpers for TSV/CSV inputs and archive bundles."""

from __future__ import annotations

import fnmatch
import sqlite3
import tarfile
from io import TextIOWrapper
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from src.models import MoleculeInfo, MoleculeInput
from src.service import get_predictor

DEFAULT_ARCHIVE_PATTERN = "*_scheme_b_unique_products.csv"
DEFAULT_ARCHIVE_SMILES_COL = "product_smiles_canonical"
DEFAULT_ARCHIVE_CHEM_ID_COL = "example_combo_id"
DEFAULT_CHUNK_SIZE = 256


def _basename_without_suffixes(path: Path) -> str:
    name = path.name
    for suffix in reversed(path.suffixes):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


def _looks_like_archive(path: Path) -> bool:
    return path.is_file() and tarfile.is_tarfile(path)


def _looks_like_sqlite(path: Path) -> bool:
    return path.is_file() and path.suffix.lower() in {".sqlite", ".sqlite3", ".db", ".db3"}


def _infer_tabular_separator(path: Path) -> str:
    suffixes = [suffix.lower() for suffix in path.suffixes]
    if suffixes[-2:] == [".tsv", ".gz"] or path.suffix.lower() in {".tsv", ".txt"}:
        return "\t"
    return ","


def _load_tabular_input(
    path: Path,
    smiles_colname: str,
    chem_id_colname: str,
) -> pd.DataFrame:
    frame = pd.read_csv(path, sep=_infer_tabular_separator(path))
    if smiles_colname not in frame.columns:
        raise ValueError(f"Missing SMILES column '{smiles_colname}' in {path}")

    if chem_id_colname not in frame.columns:
        frame[chem_id_colname] = [f"mol{i + 1}" for i in range(len(frame))]

    frame = frame.rename(columns={smiles_colname: "smiles", chem_id_colname: "chem_id"})
    frame["source_file"] = str(path)
    frame["source_group"] = _basename_without_suffixes(path)
    frame["source_row"] = range(1, len(frame) + 1)
    frame["input_order"] = range(len(frame))
    return frame


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
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
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
            frame["input_order"] = range(len(frame))
            frames.append(frame)

    if not frames:
        raise ValueError(f"No prediction rows could be loaded from {path}")

    return pd.concat(frames, ignore_index=True)


def _load_sqlite_input(
    path: Path,
    smiles_colname: str,
    chem_id_colname: str,
    sqlite_table: str | None,
    sqlite_query: str | None,
) -> pd.DataFrame:
    with sqlite3.connect(path) as connection:
        if sqlite_query:
            frame = pd.read_sql_query(sqlite_query, connection)
            source_group = sqlite_table or "query"
        else:
            if sqlite_table:
                frame = pd.read_sql_query(f'SELECT * FROM "{sqlite_table}"', connection)
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
                sqlite_table = candidates[0]
                frame = pd.read_sql_query(f'SELECT * FROM "{sqlite_table}"', connection)
                source_group = sqlite_table

    if smiles_colname not in frame.columns:
        raise ValueError(f"Missing SMILES column '{smiles_colname}' in {path}")

    if chem_id_colname not in frame.columns:
        frame[chem_id_colname] = [f"{source_group}__{i + 1}" for i in range(len(frame))]

    frame = frame.rename(columns={smiles_colname: "smiles", chem_id_colname: "chem_id"})
    frame["source_file"] = str(path)
    frame["source_group"] = str(source_group)
    frame["source_row"] = range(1, len(frame) + 1)
    frame["source_member"] = str(source_group)
    frame["input_order"] = range(len(frame))
    return frame


def _uniquify_chem_ids(frame: pd.DataFrame) -> pd.DataFrame:
    if not frame["chem_id"].duplicated().any():
        return frame

    counts: dict[str, int] = {}
    unique_ids: list[str] = []
    for raw_id in frame["chem_id"].astype(str):
        count = counts.get(raw_id, 0)
        counts[raw_id] = count + 1
        unique_ids.append(raw_id if count == 0 else f"{raw_id}__{count + 1}")

    deduped = frame.copy()
    deduped["chem_id"] = unique_ids
    return deduped


def load_screening_input(
    input_path: str | Path,
    smiles_colname: str = "smiles",
    chem_id_colname: str = "chem_id",
    archive_pattern: str = DEFAULT_ARCHIVE_PATTERN,
    archive_smiles_colname: str = DEFAULT_ARCHIVE_SMILES_COL,
    archive_chem_id_colname: str = DEFAULT_ARCHIVE_CHEM_ID_COL,
    sqlite_table: str | None = None,
    sqlite_query: str | None = None,
    dedupe_smiles: bool = True,
) -> pd.DataFrame:
    """Load a tabular screening input or archive bundle into a normalized frame."""

    path = Path(input_path).expanduser().resolve()
    if _looks_like_archive(path):
        frame = _load_archive_input(
            path=path,
            archive_pattern=archive_pattern,
            archive_smiles_colname=archive_smiles_colname,
            archive_chem_id_colname=archive_chem_id_colname,
            fallback_smiles_colname=smiles_colname,
            fallback_chem_id_colname=chem_id_colname,
        )
    elif _looks_like_sqlite(path):
        frame = _load_sqlite_input(
            path=path,
            smiles_colname=smiles_colname,
            chem_id_colname=chem_id_colname,
            sqlite_table=sqlite_table,
            sqlite_query=sqlite_query,
        )
    else:
        frame = _load_tabular_input(path, smiles_colname, chem_id_colname)

    frame = frame.dropna(subset=["smiles"]).copy()
    frame["smiles"] = frame["smiles"].astype(str)
    frame["chem_id"] = frame["chem_id"].astype(str)

    source_counts = frame["smiles"].value_counts().to_dict()
    frame["input_occurrence_count"] = frame["smiles"].map(source_counts).fillna(1).astype(int)

    if dedupe_smiles:
        frame = frame.sort_values("input_order").drop_duplicates(subset=["smiles"], keep="first")
    frame = _uniquify_chem_ids(frame)
    frame = frame.reset_index(drop=True)
    frame["input_order"] = range(len(frame))
    return frame


def _build_molecules(frame: pd.DataFrame) -> list[MoleculeInfo]:
    return [
        MoleculeInfo(smiles=row.smiles, chem_id=row.chem_id)
        for row in frame[["smiles", "chem_id"]].itertuples(index=False)
    ]


async def screen_frame(
    frame: pd.DataFrame,
    aggregate_scores: bool,
    app_threshold: float,
    min_nkill: int,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> pd.DataFrame:
    """Run batch prediction on a normalized frame and attach metadata."""

    if frame.empty:
        raise ValueError("No screening rows available")

    predictor = get_predictor()
    await predictor.ensure_loaded()

    results: list[pd.DataFrame] = []
    for start in range(0, len(frame), chunk_size):
        chunk = frame.iloc[start : start + chunk_size].copy()
        request = MoleculeInput(
            molecules=_build_molecules(chunk),
            aggregate_scores=aggregate_scores,
            app_threshold=app_threshold,
            min_nkill=min_nkill,
        ).normalize()
        items = await predictor.predict(request)
        chunk_result = pd.DataFrame(items)
        if aggregate_scores:
            chunk_result = chunk_result.merge(chunk, on="chem_id", how="left")
        else:
            split = chunk_result["pred_id"].astype(str).str.rsplit(":", n=1, expand=True)
            chunk_result["chem_id"] = split[0]
            chunk_result["strain_name"] = split[1]
            chunk_result = chunk_result.merge(chunk, on="chem_id", how="left")
        results.append(chunk_result)

    combined = pd.concat(results, ignore_index=True)
    if "source_group" not in combined.columns:
        combined["source_group"] = "input"
    if "input_order" in combined.columns:
        sort_columns = ["input_order"]
        if "strain_name" in combined.columns:
            sort_columns.append("strain_name")
        combined = combined.sort_values(sort_columns, kind="stable")
    return combined.reset_index(drop=True)


async def screen_path(
    input_path: str | Path,
    smiles_colname: str = "smiles",
    chem_id_colname: str = "chem_id",
    archive_pattern: str = DEFAULT_ARCHIVE_PATTERN,
    archive_smiles_colname: str = DEFAULT_ARCHIVE_SMILES_COL,
    archive_chem_id_colname: str = DEFAULT_ARCHIVE_CHEM_ID_COL,
    sqlite_table: str | None = None,
    sqlite_query: str | None = None,
    dedupe_smiles: bool = True,
    aggregate_scores: bool = True,
    app_threshold: float = 0.04374140128493309,
    min_nkill: int = 10,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = load_screening_input(
        input_path=input_path,
        smiles_colname=smiles_colname,
        chem_id_colname=chem_id_colname,
        archive_pattern=archive_pattern,
        archive_smiles_colname=archive_smiles_colname,
        archive_chem_id_colname=archive_chem_id_colname,
        sqlite_table=sqlite_table,
        sqlite_query=sqlite_query,
        dedupe_smiles=dedupe_smiles,
    )
    predicted = await screen_frame(
        frame=frame,
        aggregate_scores=aggregate_scores,
        app_threshold=app_threshold,
        min_nkill=min_nkill,
        chunk_size=chunk_size,
    )
    return frame, predicted
