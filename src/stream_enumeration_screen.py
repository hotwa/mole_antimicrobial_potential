"""Streaming scaffold/fragment enumeration with resumable hit-only screening."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import pandas as pd
from rdkit import Chem

from src.models import MoleculeInfo
from src.reinvent4_workflow import sanitize_name, validate_scaffold_file, validate_scaffold_smiles
from src.service import get_scheduler

DEFAULT_SCAFFOLD_FILE = (
    Path(__file__).resolve().parent.parent
    / "workflows"
    / "reinvent4"
    / "inputs"
    / "scaffolds"
    / "mother_scaffold.template.smi"
)
DEFAULT_POSITION_TO_ATTACHMENT = {
    "pos4": 0,
    "pos3": 1,
    "pos13": 2,
    "pos12": 3,
}
HIT_COLUMNS = [
    "chem_id",
    "smiles",
    "broad_spectrum",
    "ginhib_total",
    "apscore_total",
    "global_combination_index",
    "shard_id",
    "scaffold_slug",
]


def _utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


@dataclass(frozen=True)
class IndexSpace:
    scaffold_count: int
    pos3_count: int
    pos4_count: int
    pos12_count: int
    pos13_count: int

    @property
    def per_scaffold_combinations(self) -> int:
        return self.pos3_count * self.pos4_count * self.pos12_count * self.pos13_count

    @property
    def total_combinations(self) -> int:
        return self.scaffold_count * self.per_scaffold_combinations


@dataclass(frozen=True)
class ScaffoldEntry:
    scaffold_idx: int
    scaffold_slug: str
    scaffold_smiles: str
    source_path: str | None = None


@dataclass
class StreamScreenSummary:
    output_dir: Path
    run_state_path: Path
    shard_manifest_path: Path
    attempted_count: int
    hit_count: int
    completed_shards: int
    start_index: int
    stop_index: int
    total_combinations: int


def encode_global_index(parts: Sequence[int], space: IndexSpace) -> int:
    if len(parts) != 5:
        raise ValueError("expected 5 indices: scaffold, pos3, pos4, pos12, pos13")
    scaffold_idx, pos3_idx, pos4_idx, pos12_idx, pos13_idx = [int(value) for value in parts]
    if not (0 <= scaffold_idx < space.scaffold_count):
        raise ValueError("scaffold index out of range")
    if not (0 <= pos3_idx < space.pos3_count):
        raise ValueError("pos3 index out of range")
    if not (0 <= pos4_idx < space.pos4_count):
        raise ValueError("pos4 index out of range")
    if not (0 <= pos12_idx < space.pos12_count):
        raise ValueError("pos12 index out of range")
    if not (0 <= pos13_idx < space.pos13_count):
        raise ValueError("pos13 index out of range")
    index = scaffold_idx
    index = (index * space.pos3_count) + pos3_idx
    index = (index * space.pos4_count) + pos4_idx
    index = (index * space.pos12_count) + pos12_idx
    index = (index * space.pos13_count) + pos13_idx
    return index


def decode_global_index(index: int, space: IndexSpace) -> tuple[int, int, int, int, int]:
    index = int(index)
    if index < 0 or index >= space.total_combinations:
        raise ValueError("global combination index out of range")
    scaffold_idx, remainder = divmod(index, space.per_scaffold_combinations)
    pos3_idx, remainder = divmod(remainder, space.pos4_count * space.pos12_count * space.pos13_count)
    pos4_idx, remainder = divmod(remainder, space.pos12_count * space.pos13_count)
    pos12_idx, pos13_idx = divmod(remainder, space.pos13_count)
    return scaffold_idx, pos3_idx, pos4_idx, pos12_idx, pos13_idx


def _read_first_scaffold(path: Path) -> str:
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line and not line.startswith("#"):
            return line
    raise ValueError(f"No scaffold entries found in {path}")


def _load_scaffolds(
    *,
    scaffold_file: str | Path | None,
    scaffold_dir: str | Path | None,
    scaffold_catalog: str | Path | None,
) -> list[ScaffoldEntry]:
    entries: list[ScaffoldEntry] = []
    next_idx = 0

    if scaffold_catalog:
        catalog_path = Path(scaffold_catalog).expanduser().resolve()
        separator = "\t" if catalog_path.suffix.lower() in {".tsv", ".txt"} else ","
        frame = pd.read_csv(catalog_path, sep=separator)
        if "scaffold_slug" not in frame.columns:
            raise ValueError("scaffold catalog requires scaffold_slug column")
        if "scaffold_smiles" not in frame.columns and "scaffold_file" not in frame.columns:
            raise ValueError("scaffold catalog requires scaffold_smiles or scaffold_file")
        for row in frame.itertuples(index=False):
            slug = sanitize_name(str(getattr(row, "scaffold_slug")))
            smiles = getattr(row, "scaffold_smiles", None)
            source_path = None
            if smiles in (None, ""):
                source_path = str(Path(getattr(row, "scaffold_file")).expanduser().resolve())
                smiles = _read_first_scaffold(Path(source_path))
            validate_scaffold_smiles(str(smiles))
            entries.append(
                ScaffoldEntry(
                    scaffold_idx=next_idx,
                    scaffold_slug=slug,
                    scaffold_smiles=str(smiles),
                    source_path=source_path,
                )
            )
            next_idx += 1
    elif scaffold_dir:
        scaffold_root = Path(scaffold_dir).expanduser().resolve()
        for path in sorted(scaffold_root.glob("*.smi")):
            validate_scaffold_file(path)
            entries.append(
                ScaffoldEntry(
                    scaffold_idx=next_idx,
                    scaffold_slug=sanitize_name(path.stem),
                    scaffold_smiles=_read_first_scaffold(path),
                    source_path=str(path),
                )
            )
            next_idx += 1
    else:
        path = Path(scaffold_file or DEFAULT_SCAFFOLD_FILE).expanduser().resolve()
        validate_scaffold_file(path)
        entries.append(
            ScaffoldEntry(
                scaffold_idx=0,
                scaffold_slug=sanitize_name(path.stem),
                scaffold_smiles=_read_first_scaffold(path),
                source_path=str(path),
            )
        )

    if not entries:
        raise ValueError("No scaffold entries resolved")
    return entries


def _load_fragment_smiles(path: str | Path) -> list[str]:
    frame = pd.read_csv(path)
    for column_name in ("fragment_smiles_canonical", "fragment_smiles_plain", "fragment_smiles_labeled_example"):
        if column_name in frame.columns:
            series = frame[column_name].dropna().astype(str)
            return series.tolist()
    raise ValueError(f"Could not find fragment SMILES column in {path}")


def _find_scaffold_attachment_index(mol: Chem.Mol, attachment_label: int) -> int:
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() == attachment_label:
            return atom.GetIdx()
    raise ValueError(f"Missing scaffold attachment point [*:{attachment_label}]")


def _find_fragment_attachment_index(mol: Chem.Mol) -> int:
    matches = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]
    if len(matches) != 1:
        raise ValueError("Fragments must contain exactly one dummy atom")
    return matches[0]


def _attach_fragment(scaffold: Chem.Mol, attachment_label: int, fragment_smiles: str) -> Chem.Mol:
    fragment = Chem.MolFromSmiles(fragment_smiles)
    if fragment is None:
        raise ValueError(f"Invalid fragment SMILES: {fragment_smiles}")
    scaffold_idx = _find_scaffold_attachment_index(scaffold, attachment_label)
    fragment_idx = _find_fragment_attachment_index(fragment)

    scaffold_atom = scaffold.GetAtomWithIdx(scaffold_idx)
    fragment_atom = fragment.GetAtomWithIdx(fragment_idx)
    if len(scaffold_atom.GetNeighbors()) != 1 or len(fragment_atom.GetNeighbors()) != 1:
        raise ValueError("Attachment dummy atoms must be terminal")

    scaffold_neighbor = scaffold_atom.GetNeighbors()[0].GetIdx()
    fragment_neighbor = fragment_atom.GetNeighbors()[0].GetIdx()
    scaffold_bond = scaffold.GetBondBetweenAtoms(scaffold_idx, scaffold_neighbor)
    fragment_bond = fragment.GetBondBetweenAtoms(fragment_idx, fragment_neighbor)
    bond_type = Chem.BondType.SINGLE
    if scaffold_bond is not None and scaffold_bond.GetBondType() != Chem.BondType.UNSPECIFIED:
        bond_type = scaffold_bond.GetBondType()
    elif fragment_bond is not None and fragment_bond.GetBondType() != Chem.BondType.UNSPECIFIED:
        bond_type = fragment_bond.GetBondType()

    combined = Chem.RWMol(Chem.CombineMols(scaffold, fragment))
    fragment_offset = scaffold.GetNumAtoms()
    combined.AddBond(scaffold_neighbor, fragment_neighbor + fragment_offset, bond_type)
    for atom_index in sorted([scaffold_idx, fragment_idx + fragment_offset], reverse=True):
        combined.RemoveAtom(atom_index)

    result = combined.GetMol()
    Chem.SanitizeMol(result)
    return result


def _assemble_molecule(scaffold_smiles: str, fragments: Mapping[str, str]) -> str:
    scaffold = Chem.MolFromSmiles(scaffold_smiles)
    if scaffold is None:
        raise ValueError(f"Invalid scaffold SMILES: {scaffold_smiles}")
    current = scaffold
    for position_name in ("pos4", "pos3", "pos13", "pos12"):
        current = _attach_fragment(current, DEFAULT_POSITION_TO_ATTACHMENT[position_name], fragments[position_name])
    return Chem.MolToSmiles(current, canonical=True, isomericSmiles=True)


def _shard_id_for_range(start_idx: int, end_idx: int) -> str:
    return f"shard_{start_idx:020d}_{end_idx:020d}"


def _hits_output_path(output_dir: Path, shard_id: str) -> Path:
    return output_dir / "hits" / f"{shard_id}.parquet"


def _load_manifest(path: Path) -> dict[str, dict[str, Any]]:
    if not path.exists():
        return {}
    records: dict[str, dict[str, Any]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        record = json.loads(line)
        records[str(record["shard_id"])] = record
    return records


def _write_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(f"{path.name}.tmp")
    tmp_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    tmp_path.replace(path)


def _write_manifest_atomic(path: Path, records: Mapping[str, Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(f"{path.name}.tmp")
    ordered = [records[key] for key in sorted(records)]
    rendered = "".join(json.dumps(record, ensure_ascii=False) + "\n" for record in ordered)
    tmp_path.write_text(rendered, encoding="utf-8")
    tmp_path.replace(path)


def _read_optional_json(path: str | Path | None) -> Any:
    if not path:
        return None
    resolved = Path(path).expanduser().resolve()
    if not resolved.exists():
        return None
    return json.loads(resolved.read_text(encoding="utf-8"))


def _read_optional_tabular_preview(path: str | Path | None) -> list[dict[str, Any]] | None:
    if not path:
        return None
    resolved = Path(path).expanduser().resolve()
    if not resolved.exists():
        return None
    separator = "\t" if resolved.suffix.lower() in {".tsv", ".txt"} else ","
    frame = pd.read_csv(resolved, sep=separator)
    return frame.head(10).to_dict(orient="records")


def _build_parameter_snapshot(
    *,
    scaffold_file: str | Path | None,
    scaffold_dir: str | Path | None,
    scaffold_catalog: str | Path | None,
    ordinary_library: str | Path,
    pos13_library: str | Path,
    run_state_source: str | Path | None,
    chunk_manifest_source: str | Path | None,
    start_index: int,
    stop_index: int,
    shard_size: int,
    prediction_batch_size: int,
    app_threshold: float,
    min_nkill: int,
    classifier_backend: str,
    deterministic_representation: bool,
    space: IndexSpace,
    scaffolds: Sequence[ScaffoldEntry],
) -> dict[str, Any]:
    return {
        "scaffold_file": str(Path(scaffold_file or DEFAULT_SCAFFOLD_FILE).expanduser().resolve()) if scaffold_file or not scaffold_dir and not scaffold_catalog else None,
        "scaffold_dir": str(Path(scaffold_dir).expanduser().resolve()) if scaffold_dir else None,
        "scaffold_catalog": str(Path(scaffold_catalog).expanduser().resolve()) if scaffold_catalog else None,
        "ordinary_library": str(Path(ordinary_library).expanduser().resolve()),
        "pos13_library": str(Path(pos13_library).expanduser().resolve()),
        "run_state_source": str(Path(run_state_source).expanduser().resolve()) if run_state_source else None,
        "chunk_manifest_source": str(Path(chunk_manifest_source).expanduser().resolve()) if chunk_manifest_source else None,
        "start_index": int(start_index),
        "stop_index": int(stop_index),
        "shard_size": int(shard_size),
        "prediction_batch_size": int(prediction_batch_size),
        "app_threshold": float(app_threshold),
        "min_nkill": int(min_nkill),
        "classifier_backend": str(classifier_backend),
        "deterministic_representation": bool(deterministic_representation),
        "counts": asdict(space),
        "per_scaffold_total_combinations": int(space.per_scaffold_combinations),
        "total_combinations": int(space.total_combinations),
        "scaffolds": [asdict(entry) for entry in scaffolds],
    }


def _validate_existing_run_state(path: Path, parameter_snapshot: Mapping[str, Any]) -> None:
    if not path.exists():
        return
    existing = json.loads(path.read_text(encoding="utf-8"))
    if existing.get("parameters") != parameter_snapshot:
        raise ValueError(f"Existing run state at {path} does not match the current parameter snapshot")


def _plan_shards(start_index: int, stop_index: int, shard_size: int) -> list[tuple[str, int, int]]:
    if shard_size <= 0:
        raise ValueError("shard_size must be positive")
    return [
        (_shard_id_for_range(start_idx, min(stop_index, start_idx + shard_size)), start_idx, min(stop_index, start_idx + shard_size))
        for start_idx in range(start_index, stop_index, shard_size)
    ]


def _build_manifest_record(*, shard_id: str, start_idx: int, end_idx: int, output_path: Path) -> dict[str, Any]:
    return {
        "shard_id": shard_id,
        "start_idx": start_idx,
        "end_idx": end_idx,
        "last_committed_idx": start_idx - 1,
        "status": "pending",
        "attempted_count": 0,
        "hit_count": 0,
        "output_path": str(output_path),
        "updated_at": _utc_now(),
        "error": None,
    }


def _materialize_batch(
    *,
    start_idx: int,
    end_idx: int,
    space: IndexSpace,
    scaffolds: Sequence[ScaffoldEntry],
    ordinary_fragments: Sequence[str],
    pos13_fragments: Sequence[str],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for global_index in range(start_idx, end_idx):
        scaffold_idx, pos3_idx, pos4_idx, pos12_idx, pos13_idx = decode_global_index(global_index, space)
        scaffold = scaffolds[scaffold_idx]
        smiles = _assemble_molecule(
            scaffold.scaffold_smiles,
            {
                "pos3": ordinary_fragments[pos3_idx],
                "pos4": ordinary_fragments[pos4_idx],
                "pos12": ordinary_fragments[pos12_idx],
                "pos13": pos13_fragments[pos13_idx],
            },
        )
        rows.append(
            {
                "chem_id": f"{scaffold.scaffold_slug}__g{global_index}",
                "smiles": smiles,
                "global_combination_index": global_index,
                "scaffold_idx": scaffold_idx,
                "scaffold_slug": scaffold.scaffold_slug,
            }
        )
    return rows


async def stream_enumeration_screen(
    *,
    output_dir: str | Path,
    scaffold_file: str | Path | None = None,
    scaffold_dir: str | Path | None = None,
    scaffold_catalog: str | Path | None = None,
    ordinary_library: str | Path,
    pos13_library: str | Path,
    run_state_source: str | Path | None = None,
    chunk_manifest_source: str | Path | None = None,
    start_index: int = 0,
    stop_index: int | None = None,
    shard_size: int = 100000,
    prediction_batch_size: int = 1024,
    scheduler: Any | None = None,
    app_threshold: float = 0.04374140128493309,
    min_nkill: int = 10,
    classifier_backend: str = "auto",
    num_graph_workers: int | str = "auto",
    graph_batch_size: int = 1024,
    prefetch_batches: int = 2,
    deterministic_representation: bool = False,
    enable_profiling: bool = False,
    fail_after_shards: int | None = None,
) -> StreamScreenSummary:
    output_path = Path(output_dir).expanduser().resolve()
    output_path.mkdir(parents=True, exist_ok=True)
    run_state_path = output_path / "run_state.json"
    shard_manifest_path = output_path / "shard_manifest.jsonl"

    scaffolds = _load_scaffolds(
        scaffold_file=scaffold_file,
        scaffold_dir=scaffold_dir,
        scaffold_catalog=scaffold_catalog,
    )
    ordinary_fragments = _load_fragment_smiles(ordinary_library)
    pos13_fragments = _load_fragment_smiles(pos13_library)
    space = IndexSpace(
        scaffold_count=len(scaffolds),
        pos3_count=len(ordinary_fragments),
        pos4_count=len(ordinary_fragments),
        pos12_count=len(ordinary_fragments),
        pos13_count=len(pos13_fragments),
    )
    total_combinations = space.total_combinations
    effective_stop = total_combinations if stop_index is None else int(stop_index)
    if start_index < 0 or effective_stop < 0 or start_index > effective_stop or effective_stop > total_combinations:
        raise ValueError("start/stop index range is invalid for the resolved index space")

    parameter_snapshot = _build_parameter_snapshot(
        scaffold_file=scaffold_file,
        scaffold_dir=scaffold_dir,
        scaffold_catalog=scaffold_catalog,
        ordinary_library=ordinary_library,
        pos13_library=pos13_library,
        run_state_source=run_state_source,
        chunk_manifest_source=chunk_manifest_source,
        start_index=start_index,
        stop_index=effective_stop,
        shard_size=shard_size,
        prediction_batch_size=prediction_batch_size,
        app_threshold=app_threshold,
        min_nkill=min_nkill,
        classifier_backend=classifier_backend,
        deterministic_representation=deterministic_representation,
        space=space,
        scaffolds=scaffolds,
    )
    _validate_existing_run_state(run_state_path, parameter_snapshot)

    manifest_records = _load_manifest(shard_manifest_path)
    for shard_id, shard_start, shard_end in _plan_shards(start_index, effective_stop, shard_size):
        output_file = _hits_output_path(output_path, shard_id)
        manifest_records.setdefault(
            shard_id,
            _build_manifest_record(shard_id=shard_id, start_idx=shard_start, end_idx=shard_end, output_path=output_file),
        )

    scheduler = scheduler or get_scheduler(
        num_graph_workers=num_graph_workers,
        graph_batch_size=graph_batch_size,
        prefetch_batches=prefetch_batches,
        deterministic_representation=deterministic_representation,
    )

    attempted_total = 0
    hit_total = 0
    completed_shards = 0

    def _write_run_state(status: str, error: str | None = None) -> None:
        completed_records = [record for record in manifest_records.values() if record["status"] == "completed"]
        payload = {
            "status": status,
            "updated_at": _utc_now(),
            "parameters": parameter_snapshot,
            "progress": {
                "completed_shards": len(completed_records),
                "attempted_count": sum(int(record.get("attempted_count", 0)) for record in completed_records),
                "hit_count": sum(int(record.get("hit_count", 0)) for record in completed_records),
                "last_completed_idx": max((int(record.get("last_committed_idx", -1)) for record in completed_records), default=start_index - 1),
            },
            "provenance": {
                "upstream_run_state": _read_optional_json(run_state_source),
                "chunk_manifest_preview": _read_optional_tabular_preview(chunk_manifest_source),
            },
            "error": error,
        }
        _write_json_atomic(run_state_path, payload)

    _write_run_state("running")
    _write_manifest_atomic(shard_manifest_path, manifest_records)

    try:
        for shard_id, shard_start, shard_end in _plan_shards(start_index, effective_stop, shard_size):
            record = manifest_records[shard_id]
            output_file = Path(record["output_path"])
            if record.get("status") == "completed" and output_file.exists():
                attempted_total += int(record.get("attempted_count", 0))
                hit_total += int(record.get("hit_count", 0))
                completed_shards += 1
                continue

            temp_output = output_file.with_name(f"{output_file.name}.tmp")
            if temp_output.exists():
                temp_output.unlink()

            record.update(
                {
                    "status": "in_progress",
                    "last_committed_idx": shard_start - 1,
                    "attempted_count": 0,
                    "hit_count": 0,
                    "updated_at": _utc_now(),
                    "error": None,
                }
            )
            _write_manifest_atomic(shard_manifest_path, manifest_records)
            _write_run_state("running")

            shard_hits: list[pd.DataFrame] = []
            shard_attempted = 0
            for batch_start in range(shard_start, shard_end, prediction_batch_size):
                batch_end = min(shard_end, batch_start + prediction_batch_size)
                batch_rows = _materialize_batch(
                    start_idx=batch_start,
                    end_idx=batch_end,
                    space=space,
                    scaffolds=scaffolds,
                    ordinary_fragments=ordinary_fragments,
                    pos13_fragments=pos13_fragments,
                )
                shard_attempted += len(batch_rows)
                molecule_rows = [MoleculeInfo(smiles=row["smiles"], chem_id=row["chem_id"]) for row in batch_rows]
                predicted = await scheduler.predict_molecules(
                    molecules=molecule_rows,
                    aggregate_scores=True,
                    app_threshold=app_threshold,
                    min_nkill=min_nkill,
                    num_graph_workers=num_graph_workers,
                    graph_batch_size=graph_batch_size,
                    prefetch_batches=prefetch_batches,
                    enable_profiling=enable_profiling,
                    deterministic_representation=deterministic_representation,
                )
                predicted_frame = pd.DataFrame(predicted)
                metadata_frame = pd.DataFrame(batch_rows)
                merged = predicted_frame.merge(metadata_frame, on="chem_id", how="left")
                hits = merged[(merged["broad_spectrum"] == 1) | (merged["ginhib_total"] >= min_nkill)].copy()
                if not hits.empty:
                    hits["shard_id"] = shard_id
                    shard_hits.append(hits[HIT_COLUMNS].reset_index(drop=True))

            hits_frame = pd.concat(shard_hits, ignore_index=True) if shard_hits else pd.DataFrame(columns=HIT_COLUMNS)
            temp_output.parent.mkdir(parents=True, exist_ok=True)
            hits_frame.to_parquet(temp_output, index=False)
            temp_output.replace(output_file)

            record.update(
                {
                    "status": "completed",
                    "last_committed_idx": shard_end - 1,
                    "attempted_count": shard_attempted,
                    "hit_count": len(hits_frame),
                    "updated_at": _utc_now(),
                    "error": None,
                }
            )
            _write_manifest_atomic(shard_manifest_path, manifest_records)
            _write_run_state("running")

            attempted_total += shard_attempted
            hit_total += len(hits_frame)
            completed_shards += 1

            if fail_after_shards is not None and completed_shards >= fail_after_shards:
                raise RuntimeError("Simulated failure after shard commit")

    except BaseException as exc:
        terminal_status = "interrupted" if isinstance(exc, KeyboardInterrupt) else "failed"
        failing_record = next(
            (record for record in manifest_records.values() if record.get("status") == "in_progress"),
            None,
        )
        if failing_record is not None:
            failing_record.update(
                {
                    "status": terminal_status,
                    "updated_at": _utc_now(),
                    "error": str(exc),
                }
            )
        _write_manifest_atomic(shard_manifest_path, manifest_records)
        _write_run_state(terminal_status, error=str(exc))
        raise

    _write_manifest_atomic(shard_manifest_path, manifest_records)
    _write_run_state("completed")
    return StreamScreenSummary(
        output_dir=output_path,
        run_state_path=run_state_path,
        shard_manifest_path=shard_manifest_path,
        attempted_count=attempted_total,
        hit_count=hit_total,
        completed_shards=completed_shards,
        start_index=start_index,
        stop_index=effective_stop,
        total_combinations=total_combinations,
    )
