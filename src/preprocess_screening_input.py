from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

REQUIRED_PARQUET_COLUMNS = ["smiles", "chem_id", "source_group", "source_file", "source_row"]


def _infer_separator(path: Path) -> str:
    suffixes = [suffix.lower() for suffix in path.suffixes]
    if suffixes[-2:] == [".tsv", ".gz"] or path.suffix.lower() in {".tsv", ".txt"}:
        return "\t"

    with path.open("r", encoding="utf-8", newline="") as handle:
        sample = handle.read(4096)
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t;").delimiter
    except csv.Error:
        return ","


def _normalize_shard(
    frame: pd.DataFrame,
    *,
    smiles_colname: str,
    chem_id_colname: str,
    source_group: str,
    source_file: str,
    row_offset: int,
    chem_counter: int,
) -> tuple[pd.DataFrame, int]:
    shard = pd.DataFrame(
        {
            "smiles": frame[smiles_colname],
            "source_row": range(row_offset + 1, row_offset + 1 + len(frame)),
        }
    )
    if chem_id_colname in frame.columns:
        shard["chem_id"] = frame[chem_id_colname]

    shard = shard.dropna(subset=["smiles"]).reset_index(drop=True)

    if "chem_id" not in shard.columns:
        shard["chem_id"] = [f"{source_group}__{chem_counter + i}" for i in range(len(shard))]
        chem_counter += len(shard)

    missing_mask = shard["chem_id"].isna() | (shard["chem_id"].astype(str).str.len() == 0)
    if missing_mask.any():
        replacement_count = int(missing_mask.sum())
        replacements = [f"{source_group}__{chem_counter + i}" for i in range(replacement_count)]
        shard.loc[missing_mask, "chem_id"] = replacements
        chem_counter += replacement_count

    shard["smiles"] = shard["smiles"].astype(str)
    shard["chem_id"] = shard["chem_id"].astype(str)
    shard["source_group"] = source_group
    shard["source_file"] = source_file

    return shard[REQUIRED_PARQUET_COLUMNS], chem_counter


def preprocess_to_parquet(
    input_path: str | Path,
    output_dir: str | Path,
    smiles_colname: str,
    chem_id_colname: str,
    source_group: str,
    rows_per_shard: int,
    row_group_size: int,
) -> dict[str, Any]:
    input_path = Path(input_path).expanduser().resolve()
    output_dir = Path(output_dir).expanduser().resolve()
    target_dir = output_dir / source_group
    target_dir.mkdir(parents=True, exist_ok=True)

    if rows_per_shard <= 0:
        raise ValueError("rows_per_shard must be > 0")
    if row_group_size <= 0:
        raise ValueError("row_group_size must be > 0")

    separator = _infer_separator(input_path)
    rows_written = 0
    shard_count = 0
    chem_counter = 1
    row_offset = 0

    reader = pd.read_csv(input_path, sep=separator, chunksize=rows_per_shard)
    try:
        for frame in reader:
            if smiles_colname not in frame.columns:
                raise ValueError(f"Missing SMILES column '{smiles_colname}' in {input_path}")

            shard, chem_counter = _normalize_shard(
                frame,
                smiles_colname=smiles_colname,
                chem_id_colname=chem_id_colname,
                source_group=source_group,
                source_file=str(input_path),
                row_offset=row_offset,
                chem_counter=chem_counter,
            )
            row_offset += len(frame)

            if shard.empty:
                continue

            shard_count += 1
            shard_path = target_dir / f"shard_{shard_count:04d}.parquet"
            table = pa.Table.from_pandas(shard, preserve_index=False)
            pq.write_table(table, shard_path, row_group_size=row_group_size)
            rows_written += len(shard)
    finally:
        reader.close()

    return {
        "input_path": str(input_path),
        "output_dir": str(output_dir),
        "source_group": source_group,
        "rows_written": rows_written,
        "shard_count": shard_count,
        "row_group_size": row_group_size,
        "rows_per_shard": rows_per_shard,
    }


def preprocess_many_to_parquet(
    input_paths: list[str | Path],
    output_dir: str | Path,
    smiles_colname: str,
    chem_id_colname: str,
    rows_per_shard: int,
    row_group_size: int,
) -> dict[str, Any]:
    manifests = []
    for raw_input_path in input_paths:
        input_path = Path(raw_input_path)
        manifests.append(
            preprocess_to_parquet(
                input_path=input_path,
                output_dir=output_dir,
                smiles_colname=smiles_colname,
                chem_id_colname=chem_id_colname,
                source_group=input_path.stem,
                rows_per_shard=rows_per_shard,
                row_group_size=row_group_size,
            )
        )

    return {
        "source_groups": sorted(item["source_group"] for item in manifests),
        "inputs": manifests,
    }
