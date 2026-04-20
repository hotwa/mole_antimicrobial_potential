"""Helpers for integrating REINVENT4 with the local antimicrobial scorer."""

from __future__ import annotations

import csv
import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path
from string import Template
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence
from urllib import error, request

DEFAULT_APP_THRESHOLD = 0.04374140128493309
DEFAULT_MIN_NKILL = 10
DEFAULT_TAU = 0.02
DEFAULT_REINVENT_RUNNER = "reinvent"
DEFAULT_DEVICE = "cuda:0"
DEFAULT_BATCH_SIZE = 64
DEFAULT_SIGMA = 128
DEFAULT_LEARNING_RATE = 0.0001
DEFAULT_MAX_SCORE = 0.8
DEFAULT_MIN_STEPS = 25
DEFAULT_MAX_STEPS = 500
DEFAULT_BUCKET_SIZE = 25
DEFAULT_DIVERSITY_MINSCORE = 0.4
DEFAULT_TIMEOUT_SECONDS = 120
DEFAULT_TOP_N = 50

ATTACHMENT_POINT_PATTERN = re.compile(r"\[\*:([0-9]+)\]")


def normalize_site_reward_spec(data: Any, source: str = "<memory>") -> Dict[str, Any] | None:
    """Normalize the optional experimental site_reward block."""
    if data is None:
        return None
    if not isinstance(data, Mapping):
        raise ValueError(f"site_reward must be a JSON object: {source}")

    normalized = {
        "enabled": bool(data.get("enabled", False)),
        "range_min": parse_int(data.get("range_min"), 4),
        "range_max": parse_int(data.get("range_max"), 27),
        "alpha": parse_float(data.get("alpha"), 1.5),
        "beta": parse_float(data.get("beta"), 2.5),
        "coverage_weight": parse_float(data.get("coverage_weight"), 0.7),
        "balance_weight": parse_float(data.get("balance_weight"), 0.3),
        "lambda": parse_float(data.get("lambda"), 0.85),
    }
    scaffold_smiles = data.get("scaffold_smiles")
    if scaffold_smiles is not None:
        cleaned = str(scaffold_smiles).strip()
        if not cleaned:
            raise ValueError("site_reward.scaffold_smiles must be a non-empty string")
        normalized["scaffold_smiles"] = cleaned
    if normalized["range_min"] < 0 or normalized["range_max"] < normalized["range_min"]:
        raise ValueError("site_reward requires 0 <= range_min <= range_max")
    if normalized["alpha"] <= 0 or normalized["beta"] <= 0:
        raise ValueError("site_reward alpha and beta must be positive")
    if normalized["coverage_weight"] < 0 or normalized["balance_weight"] < 0:
        raise ValueError("site_reward weights must be non-negative")
    if (normalized["coverage_weight"] + normalized["balance_weight"]) <= 0:
        raise ValueError("site_reward coverage_weight + balance_weight must be positive")
    if not (0 <= normalized["lambda"] <= 1):
        raise ValueError("site_reward lambda must be within [0, 1]")
    return normalized


def repo_root() -> Path:
    """Return the repository root."""
    return Path(__file__).resolve().parent.parent


def parse_env_file(env_path: str | Path) -> Dict[str, str]:
    """Parse a minimal .env style file."""
    path = Path(env_path).expanduser().resolve()
    values: Dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("export "):
            line = line[len("export ") :].strip()
        if "=" not in line:
            raise ValueError(f"Invalid env line in {path}: {raw_line}")
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip().strip("'").strip('"')
        values[key] = os.path.expandvars(value)
    return values


def resolve_path(value: str | Path, base_dir: str | Path | None = None) -> Path:
    """Resolve a path relative to an optional base directory."""
    path = Path(value).expanduser()
    if not path.is_absolute() and base_dir is not None:
        path = Path(base_dir).expanduser() / path
    return path.resolve()


def parse_float(value: Any, default: float) -> float:
    """Parse a float with a default."""
    if value in (None, ""):
        return default
    return float(value)


def parse_int(value: Any, default: int) -> int:
    """Parse an int with a default."""
    if value in (None, ""):
        return default
    return int(value)


def load_runtime_settings(env_path: str | Path) -> Dict[str, Any]:
    """Load and normalize runtime settings for REINVENT4 runs."""
    env_file = resolve_path(env_path)
    values = parse_env_file(env_file)
    required = ["REINVENT4_ROOT", "REINVENT4_PRIOR_FILE", "MOLE_SCORE_API_URL"]
    missing = [key for key in required if key not in values]
    if missing:
        raise ValueError(f"Missing required settings in {env_file}: {missing}")

    reinvent_root = resolve_path(values["REINVENT4_ROOT"], env_file.parent)
    prior_file = resolve_path(values["REINVENT4_PRIOR_FILE"], env_file.parent)
    agent_file = resolve_path(values.get("REINVENT4_AGENT_FILE", prior_file), env_file.parent)
    if not reinvent_root.is_dir():
        raise ValueError(f"REINVENT4_ROOT does not exist or is not a directory: {reinvent_root}")
    if not prior_file.is_file():
        raise ValueError(f"REINVENT4_PRIOR_FILE does not exist: {prior_file}")
    if not agent_file.is_file():
        raise ValueError(f"REINVENT4_AGENT_FILE does not exist: {agent_file}")

    settings: Dict[str, Any] = {
        "env_file": env_file,
        "reinvent_root": reinvent_root,
        "reinvent_runner": values.get("REINVENT4_RUNNER", DEFAULT_REINVENT_RUNNER),
        "prior_file": prior_file,
        "agent_file": agent_file,
        "score_url": values["MOLE_SCORE_API_URL"],
        "api_health_url": values.get("MOLE_API_HEALTH_URL", infer_health_url(values["MOLE_SCORE_API_URL"])),
        "device": values.get("REINVENT_DEVICE", DEFAULT_DEVICE),
        "batch_size": parse_int(values.get("REINVENT_BATCH_SIZE"), DEFAULT_BATCH_SIZE),
        "sigma": parse_float(values.get("REINVENT_SIGMA"), DEFAULT_SIGMA),
        "learning_rate": parse_float(values.get("REINVENT_LEARNING_RATE"), DEFAULT_LEARNING_RATE),
        "max_score": parse_float(values.get("REINVENT_MAX_SCORE"), DEFAULT_MAX_SCORE),
        "min_steps": parse_int(values.get("REINVENT_MIN_STEPS"), DEFAULT_MIN_STEPS),
        "max_steps": parse_int(values.get("REINVENT_MAX_STEPS"), DEFAULT_MAX_STEPS),
        "bucket_size": parse_int(values.get("REINVENT_BUCKET_SIZE"), DEFAULT_BUCKET_SIZE),
        "diversity_minscore": parse_float(
            values.get("REINVENT_DIVERSITY_MINSCORE"), DEFAULT_DIVERSITY_MINSCORE
        ),
        "request_timeout_seconds": parse_int(
            values.get("MOLE_REQUEST_TIMEOUT_SECONDS"), DEFAULT_TIMEOUT_SECONDS
        ),
    }
    return settings


def infer_health_url(score_url: str) -> str:
    """Infer the FastAPI /health URL from the /score URL."""
    if score_url.endswith("/score"):
        return f"{score_url[:-len('/score')]}/health"
    return f"{score_url.rstrip('/')}/health"


def sanitize_name(name: str) -> str:
    """Convert free text into a filesystem-friendly name."""
    lowered = name.strip().lower()
    cleaned = re.sub(r"[^a-z0-9]+", "_", lowered)
    return cleaned.strip("_") or "experiment"


def normalize_objective_spec(data: Mapping[str, Any], source: str = "<memory>") -> Dict[str, Any]:
    """Normalize a REINVENT4 objective spec from an in-memory mapping."""
    if not isinstance(data, Mapping):
        raise ValueError(f"Objective spec must be a JSON object: {source}")
    mode = data.get("mode")
    if mode not in {"single_strain", "broad_spectrum_soft"}:
        raise ValueError(f"Unsupported objective mode in {source}: {mode}")

    strain = data.get("strain")
    strains = data.get("strains")
    strain_index = data.get("strain_index")
    strain_indices = data.get("strain_indices")
    weights = data.get("weights")
    if strain and strains:
        raise ValueError("Objective spec must not contain both 'strain' and 'strains'")
    if strain_index is not None and strain_indices is not None:
        raise ValueError("Objective spec must not contain both 'strain_index' and 'strain_indices'")
    if strain and strain_index is not None:
        raise ValueError("Objective spec must not contain both 'strain' and 'strain_index'")
    if strains and strain_indices is not None:
        raise ValueError("Objective spec must not contain both 'strains' and 'strain_indices'")

    if strain_index is not None:
        strain = strain_name_from_index(strain_index)
    if strain_indices is not None:
        strains = strain_names_from_indices(strain_indices)

    if mode == "single_strain" and not (strain or strains):
        raise ValueError("single_strain objective requires 'strain', 'strains', 'strain_index', or 'strain_indices'")
    if mode == "broad_spectrum_soft" and weights is not None:
        raise ValueError("broad_spectrum_soft objective does not support weights")

    if strains is not None:
        if not isinstance(strains, list) or not strains:
            raise ValueError("'strains' must be a non-empty list")
        if len(set(strains)) != len(strains):
            raise ValueError("'strains' must contain unique entries")
    if weights is not None:
        if not isinstance(weights, list) or not weights:
            raise ValueError("'weights' must be a non-empty list")
        if strains is None:
            raise ValueError("'weights' require 'strains'")
        if len(weights) != len(strains):
            raise ValueError("'weights' length must match 'strains' length")
        if any(float(weight) < 0 for weight in weights):
            raise ValueError("'weights' must be non-negative")
        if sum(float(weight) for weight in weights) <= 0:
            raise ValueError("'weights' must sum to a positive value")

    normalized: Dict[str, Any] = {
        "mode": mode,
        "app_threshold": parse_float(data.get("app_threshold"), DEFAULT_APP_THRESHOLD),
        "min_nkill": parse_int(data.get("min_nkill"), DEFAULT_MIN_NKILL),
        "tau": parse_float(data.get("tau"), DEFAULT_TAU),
    }
    if strain:
        normalized["strain"] = str(strain).strip()
    if strains:
        normalized["strains"] = [str(item).strip() for item in strains]
    if weights:
        normalized["weights"] = [float(weight) for weight in weights]
    site_reward = normalize_site_reward_spec(data.get("site_reward"), source=source)
    if site_reward is not None:
        normalized["site_reward"] = site_reward
    normalized["label"] = data.get("label", default_objective_label(normalized))
    return normalized


def load_objective_spec(objective_path: str | Path) -> Dict[str, Any]:
    """Load and validate a REINVENT4 objective spec."""
    path = resolve_path(objective_path)
    data = json.loads(path.read_text(encoding="utf-8"))
    return normalize_objective_spec(data, source=str(path))


def load_long_run_spec(spec_path: str | Path) -> Dict[str, Any]:
    """Load and validate a long-running RL experiment spec."""
    path = resolve_path(spec_path)
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Long-run spec must be a JSON object: {path}")

    experiment_name = str(data.get("experiment_name", "")).strip()
    if not experiment_name:
        raise ValueError(f"Long-run spec requires 'experiment_name': {path}")

    objective = normalize_objective_spec(data.get("objective", {}), source=f"{path}:objective")
    chunk_steps = parse_int(data.get("chunk_steps"), DEFAULT_BUCKET_SIZE)
    min_total_steps = parse_int(data.get("min_total_steps"), chunk_steps)
    max_total_steps = parse_int(data.get("max_total_steps"), chunk_steps)
    success_high = parse_float(data.get("success_high"), 0.8)
    success_top5_mean = parse_float(data.get("success_top5_mean"), 0.6)
    success_consecutive_chunks = parse_int(data.get("success_consecutive_chunks"), 3)
    monitor_top_n = parse_int(data.get("monitor_top_n"), 20)

    if chunk_steps <= 0:
        raise ValueError("'chunk_steps' must be a positive integer")
    if min_total_steps <= 0 or max_total_steps <= 0:
        raise ValueError("'min_total_steps' and 'max_total_steps' must be positive integers")
    if min_total_steps > max_total_steps:
        raise ValueError("'min_total_steps' must be <= 'max_total_steps'")
    if min_total_steps % chunk_steps != 0 or max_total_steps % chunk_steps != 0:
        raise ValueError("'min_total_steps' and 'max_total_steps' must be divisible by 'chunk_steps'")
    if success_consecutive_chunks <= 0:
        raise ValueError("'success_consecutive_chunks' must be a positive integer")
    if monitor_top_n <= 0:
        raise ValueError("'monitor_top_n' must be a positive integer")

    commensal_monitor_indices = data.get("commensal_monitor_indices", [])
    commensal_monitor_strains: List[str] = []
    if commensal_monitor_indices:
        if not isinstance(commensal_monitor_indices, list):
            raise ValueError("'commensal_monitor_indices' must be a list")
        commensal_monitor_strains = strain_names_from_indices(commensal_monitor_indices)

    normalized = {
        "experiment_name": experiment_name,
        "objective": objective,
        "chunk_steps": chunk_steps,
        "min_total_steps": min_total_steps,
        "max_total_steps": max_total_steps,
        "success_high": success_high,
        "success_top5_mean": success_top5_mean,
        "success_consecutive_chunks": success_consecutive_chunks,
        "monitor_top_n": monitor_top_n,
        "commensal_monitor_indices": list(commensal_monitor_indices),
        "commensal_monitor_strains": commensal_monitor_strains,
        "max_chunks": max_total_steps // chunk_steps,
        "min_chunks": min_total_steps // chunk_steps,
    }
    return normalized


def strain_index_file() -> Path:
    """Return the curated strain index table path."""
    return repo_root() / "workflows" / "reinvent4" / "inputs" / "strain_index.tsv"


def read_strain_index_mapping() -> Dict[int, str]:
    """Read the index -> strain_name mapping table."""
    path = strain_index_file()
    if path.is_file():
        mapping: Dict[int, str] = {}
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if "index" not in reader.fieldnames or "strain_name" not in reader.fieldnames:
                raise ValueError(f"Invalid strain index table header in {path}")
            for row in reader:
                raw_index = (row.get("index") or "").strip()
                raw_name = (row.get("strain_name") or "").strip()
                if not raw_index or not raw_name:
                    continue
                mapping[int(raw_index)] = raw_name
        if mapping:
            return mapping

    strains = read_all_strains_from_data()
    return {index: strain_name for index, strain_name in enumerate(strains, start=1)}


def parse_strain_index(value: Any) -> int:
    """Parse a user-provided strain index alias."""
    if isinstance(value, bool):
        raise ValueError("strain index must be a positive integer")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid strain index: {value!r}") from exc
    if parsed <= 0:
        raise ValueError("strain index must be a positive integer")
    return parsed


def strain_name_from_index(value: Any) -> str:
    """Resolve a 1-based strain index to the canonical strain name."""
    index = parse_strain_index(value)
    mapping = read_strain_index_mapping()
    if index not in mapping:
        raise ValueError(f"Unknown strain index: {index}")
    return mapping[index]


def strain_names_from_indices(values: Sequence[Any]) -> List[str]:
    """Resolve a list of 1-based strain indices to canonical strain names."""
    if not isinstance(values, list) or not values:
        raise ValueError("'strain_indices' must be a non-empty list")
    names = [strain_name_from_index(value) for value in values]
    if len(set(names)) != len(names):
        raise ValueError("'strain_indices' must resolve to unique strains")
    return names


def default_objective_label(objective: Mapping[str, Any]) -> str:
    """Create a short label for an objective spec."""
    if objective["mode"] == "single_strain":
        if "strain" in objective:
            return f"single_{sanitize_name(objective['strain'])}"
        return f"multi_{len(objective.get('strains', []))}_strain"
    if objective.get("strains"):
        return f"broad_subset_{len(objective['strains'])}"
    return "broad_all_strains"


def validate_scaffold_smiles(scaffold: str) -> List[int]:
    """Validate a LibInvent scaffold string and return attachment labels."""
    labels = [int(match.group(1)) for match in ATTACHMENT_POINT_PATTERN.finditer(scaffold)]
    if not labels:
        raise ValueError("Scaffold must contain numbered attachment points like [*:0]")
    if scaffold.count("*") != len(labels):
        raise ValueError("All attachment points must be numbered as [*:n]")
    if len(labels) != len(set(labels)):
        raise ValueError("Attachment point labels must be unique")
    expected = list(range(len(labels)))
    if sorted(labels) != expected:
        raise ValueError(f"Attachment point labels must form a contiguous set starting at 0: {expected}")

    try:
        from rdkit import Chem
    except ImportError:  # pragma: no cover - optional
        return labels

    mol = Chem.MolFromSmiles(scaffold)
    if mol is None:
        raise ValueError(f"Scaffold is not a valid SMILES: {scaffold}")
    seen_attachment_targets = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue
        neighbours = atom.GetNeighbors()
        if len(neighbours) != 1:
            raise ValueError("Attachment points must be terminal dummy atoms")
        neighbour_idx = neighbours[0].GetIdx()
        if neighbour_idx in seen_attachment_targets:
            raise ValueError("Multiple attachment points on the same atom are not supported by LibInvent")
        seen_attachment_targets.add(neighbour_idx)
    return labels


def validate_scaffold_file(scaffold_path: str | Path) -> List[Dict[str, Any]]:
    """Validate every non-comment scaffold in a file."""
    path = resolve_path(scaffold_path)
    records: List[Dict[str, Any]] = []
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        scaffold = raw_line.strip()
        if not scaffold or scaffold.startswith("#"):
            continue
        labels = validate_scaffold_smiles(scaffold)
        records.append(
            {
                "line_number": line_number,
                "scaffold": scaffold,
                "attachment_points": labels,
                "attachment_count": len(labels),
            }
        )
    if not records:
        raise ValueError(f"No scaffold entries found in {path}")
    return records


def toml_literal(value: Any) -> str:
    """Render a Python value as a TOML-compatible literal."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    return json.dumps(str(value))


def render_template_file(template_path: str | Path, context: Mapping[str, Any]) -> str:
    """Render a template file using string.Template."""
    template = Template(resolve_path(template_path).read_text(encoding="utf-8"))
    return template.substitute({key: str(value) for key, value in context.items()})


def build_reinvent_context(
    settings: Mapping[str, Any],
    objective: Mapping[str, Any],
    scaffold_path: str | Path,
    objective_path: str | Path,
    run_dir: str | Path,
) -> Dict[str, Any]:
    """Build the variable context used to render REINVENT4 config templates."""
    run_path = resolve_path(run_dir)
    bridge_script = repo_root() / "workflows" / "reinvent4" / "scripts" / "score_bridge.py"
    checkpoint_file = run_path / "stage1.chkpt"
    context = {
        "device": toml_literal(settings["device"]),
        "tb_logdir": toml_literal(str(run_path / "tb_logs")),
        "prior_file": toml_literal(str(settings["prior_file"])),
        "agent_file": toml_literal(str(settings["agent_file"])),
        "smiles_file": toml_literal(str(resolve_path(scaffold_path))),
        "summary_csv_prefix": toml_literal(str(run_path / "rl_run")),
        "batch_size": settings["batch_size"],
        "sigma": settings["sigma"],
        "learning_rate": settings["learning_rate"],
        "bucket_size": settings["bucket_size"],
        "diversity_minscore": settings["diversity_minscore"],
        "max_score": settings["max_score"],
        "min_steps": settings["min_steps"],
        "max_steps": settings["max_steps"],
        "chkpt_file": toml_literal(str(checkpoint_file)),
        "score_component_name": toml_literal("Antimicrobial Reward"),
        "score_component_weight": 1.0,
        "score_bridge_executable": toml_literal(sys.executable),
        "score_bridge_args": toml_literal(
            " ".join(
                shlex.quote(arg)
                for arg in [
                    str(bridge_script),
                    "--score-url",
                    settings["score_url"],
                    "--objective-file",
                    str(resolve_path(objective_path)),
                    "--scaffold-file",
                    str(resolve_path(scaffold_path)),
                    "--audit-file",
                    str(run_path / "bridge_audit.jsonl"),
                    "--request-timeout",
                    str(settings["request_timeout_seconds"]),
                ]
            )
        ),
        "score_property": toml_literal("predictions"),
        "objective_label": toml_literal(objective["label"]),
    }
    return context


def build_score_request(smiles: Sequence[str], objective: Mapping[str, Any]) -> Dict[str, Any]:
    """Create the request body expected by POST /score."""
    chem_ids = [f"mol{i + 1}" for i in range(len(smiles))]
    objective_body = {"mode": objective["mode"], "tau": objective["tau"]}
    if "strain" in objective:
        objective_body["strain"] = objective["strain"]
    if "strains" in objective:
        objective_body["strains"] = objective["strains"]
    if "weights" in objective:
        objective_body["weights"] = objective["weights"]
    if "site_reward" in objective:
        objective_body["site_reward"] = dict(objective["site_reward"])
    return {
        "smiles": list(smiles),
        "chem_id": chem_ids,
        "objective": objective_body,
        "app_threshold": objective["app_threshold"],
        "min_nkill": objective["min_nkill"],
    }


def order_score_items(items: Sequence[Mapping[str, Any]], chem_ids: Sequence[str]) -> List[Dict[str, Any]]:
    """Reorder /score response items to match the input chem_id order."""
    keyed = {}
    for item in items:
        chem_id = item.get("chem_id")
        if chem_id in keyed:
            raise ValueError(f"Duplicate chem_id in score response: {chem_id}")
        keyed[chem_id] = dict(item)
    ordered = []
    missing = []
    for chem_id in chem_ids:
        if chem_id not in keyed:
            missing.append(chem_id)
            continue
        ordered.append(keyed[chem_id])
    if missing:
        raise ValueError(f"Score response missing chem_ids: {missing}")
    return ordered


def external_process_payload(
    score_response: Mapping[str, Any], chem_ids: Sequence[str]
) -> Dict[str, Any]:
    """Convert a /score response into REINVENT4 ExternalProcess JSON."""
    ordered_items = order_score_items(score_response.get("items", []), chem_ids)
    objective_mode = score_response.get("objective", {}).get("mode")
    return {
        "version": 1,
        "payload": {
            "predictions": [float(item["score"]) for item in ordered_items],
            "chem_ids": list(chem_ids),
            "score_details": ordered_items,
            # ExternalProcess metadata must be per-SMILES lists. Scalar strings
            # silently truncate the associated results object in REINVENT4.
            "objective_mode": [objective_mode] * len(ordered_items),
        },
    }


def http_json_request(
    url: str,
    payload: Mapping[str, Any],
    timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS,
) -> Dict[str, Any]:
    """POST a JSON payload and parse a JSON response."""
    body = json.dumps(payload).encode("utf-8")
    req = request.Request(
        url,
        data=body,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    try:
        with request.urlopen(req, timeout=timeout_seconds) as response:
            return json.loads(response.read().decode("utf-8"))
    except error.HTTPError as exc:  # pragma: no cover - exercised via higher level integration
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"HTTP {exc.code} calling {url}: {detail}") from exc
    except error.URLError as exc:  # pragma: no cover - exercised via higher level integration
        raise RuntimeError(f"Failed to call {url}: {exc}") from exc


def http_health_check(url: str, timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS) -> Dict[str, Any]:
    """GET a health endpoint and parse JSON."""
    req = request.Request(url, headers={"Accept": "application/json"}, method="GET")
    try:
        with request.urlopen(req, timeout=timeout_seconds) as response:
            return json.loads(response.read().decode("utf-8"))
    except error.HTTPError as exc:  # pragma: no cover
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"HTTP {exc.code} calling {url}: {detail}") from exc
    except error.URLError as exc:  # pragma: no cover
        raise RuntimeError(f"Failed to call {url}: {exc}") from exc


def build_reinvent_command(settings: Mapping[str, Any], config_path: str | Path) -> List[str]:
    """Build the subprocess command used to launch REINVENT4."""
    runner = settings["reinvent_runner"]
    args = shlex.split(runner)
    return args + [str(resolve_path(config_path))]


def run_reinvent_command(
    settings: Mapping[str, Any],
    config_path: str | Path,
    dry_run: bool = False,
) -> subprocess.CompletedProcess[Any]:
    """Run or preview a REINVENT4 command."""
    command = build_reinvent_command(settings, config_path)
    if dry_run:
        return subprocess.CompletedProcess(command, 0)
    return subprocess.run(
        command,
        check=True,
        cwd=str(settings["reinvent_root"]),
    )


def create_run_directory(
    output_root: str | Path, experiment_name: str, run_id: Optional[str] = None
) -> Path:
    """Create a deterministic run directory."""
    root = resolve_path(output_root)
    if run_id is None:
        from datetime import datetime

        run_id = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    run_dir = root / sanitize_name(experiment_name) / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


def read_reinvent_csv(csv_path: str | Path) -> List[Dict[str, str]]:
    """Read a REINVENT4 CSV output file."""
    path = resolve_path(csv_path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _score_sort_value(row: Mapping[str, str]) -> float:
    for key in ("total_score", "Antimicrobial Reward", "Antimicrobial Reward_raw"):
        if key in row and row[key] not in (None, ""):
            return float(row[key])
    return 0.0


def extract_top_unique_smiles(
    rows: Sequence[Mapping[str, str]], top_n: int = DEFAULT_TOP_N
) -> List[Dict[str, Any]]:
    """Select the top-N unique molecules from a REINVENT4 CSV output."""
    if not rows:
        raise ValueError("No rows available in REINVENT4 CSV output")
    sorted_rows = sorted(rows, key=_score_sort_value, reverse=True)
    seen = set()
    selected: List[Dict[str, Any]] = []
    for row in sorted_rows:
        smiles = row.get("SMILES") or row.get("smiles")
        if not smiles or smiles in seen:
            continue
        seen.add(smiles)
        selected.append(
            {
                "smiles": smiles,
                "reinvent_total_score": _score_sort_value(row),
                "step": row.get("step", ""),
                "input_scaffold": row.get("Input_Scaffold", ""),
            }
        )
        if len(selected) >= top_n:
            break
    return selected


def read_all_strains_from_data() -> List[str]:
    """Read the strain list from the local screening data header."""
    import gzip

    data_path = repo_root() / "data" / "01.prepare_training_data" / "maier_screening_results.tsv.gz"
    with gzip.open(data_path, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    return header[1:]


def summarize_predictions(
    selected: Sequence[Mapping[str, Any]],
    score_response: Mapping[str, Any],
    aggregate_response: Mapping[str, Any],
) -> List[Dict[str, Any]]:
    """Merge REINVENT ranking, /score output, and aggregate /predict output."""
    score_by_chem = {item["chem_id"]: item for item in score_response.get("items", [])}
    aggregate_by_chem = {item["chem_id"]: item for item in aggregate_response.get("items", [])}
    summary_rows: List[Dict[str, Any]] = []
    for index, item in enumerate(selected, start=1):
        chem_id = f"mol{index}"
        score_item = score_by_chem[chem_id]
        aggregate_item = aggregate_by_chem[chem_id]
        summary_rows.append(
            {
                "chem_id": chem_id,
                "smiles": item["smiles"],
                "reinvent_total_score": item["reinvent_total_score"],
                "reward_score": score_item["score"],
                "objective_mode": score_item["objective_mode"],
                "selected_probabilities_json": json.dumps(
                    score_item.get("selected_probabilities", {}), ensure_ascii=False, sort_keys=True
                ),
                "soft_inhibition_count": score_item.get("soft_inhibition_count", ""),
                "apscore_total": aggregate_item.get("apscore_total", ""),
                "ginhib_total": aggregate_item.get("ginhib_total", ""),
                "broad_spectrum": aggregate_item.get("broad_spectrum", ""),
            }
        )
    return summary_rows


def write_json_file(path: str | Path, payload: Any) -> None:
    """Write JSON with UTF-8 encoding."""
    output_path = Path(path)
    output_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def write_tsv(path: str | Path, rows: Sequence[Mapping[str, Any]]) -> None:
    """Write a TSV file from a list of dict rows."""
    output_path = Path(path)
    if not rows:
        output_path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0].keys())
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def score_value_from_row(row: Mapping[str, Any], key: str = "Antimicrobial Reward") -> float:
    """Extract a numeric reward-like field from a REINVENT4 CSV row."""
    raw_value = row.get(key)
    if raw_value in (None, ""):
        raw_value = row.get("Score", row.get("total_score", 0.0))
    return float(raw_value)


def zero_score_fraction(rows: Sequence[Mapping[str, Any]], key: str = "Score") -> float:
    """Return the fraction of rows with an exact zero score-like value."""
    if not rows:
        return 0.0
    zero_count = 0
    for row in rows:
        raw_value = row.get(key, row.get("total_score", 0.0))
        if raw_value in (None, ""):
            raw_value = 0.0
        if float(raw_value) == 0.0:
            zero_count += 1
    return zero_count / len(rows)


def unique_smiles_ratio(rows: Sequence[Mapping[str, Any]]) -> float:
    """Return the unique-SMILES ratio within a REINVENT4 CSV batch."""
    if not rows:
        return 0.0
    smiles = [row.get("SMILES") or row.get("smiles") for row in rows]
    normalized = [item for item in smiles if item]
    if not normalized:
        return 0.0
    return len(set(normalized)) / len(normalized)


def unique_murcko_count(smiles_list: Sequence[str]) -> int:
    """Count unique Murcko scaffolds for a set of SMILES strings."""
    if not smiles_list:
        return 0
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold
    except ImportError:  # pragma: no cover - optional dependency guard
        return 0

    scaffolds = set()
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
        if scaffold:
            scaffolds.add(scaffold)
    return len(scaffolds)


def chunk_reward_summary(
    rows: Sequence[Mapping[str, Any]],
    reward_key: str = "Antimicrobial Reward",
    top_n: int = 5,
) -> Dict[str, float]:
    """Summarize a chunk of REINVENT4 rows for long-run RL monitoring."""
    if not rows:
        raise ValueError("No REINVENT4 rows available for chunk summary")
    rewards = [score_value_from_row(row, key=reward_key) for row in rows]
    sorted_rewards = sorted(rewards, reverse=True)
    top_values = sorted_rewards[: min(top_n, len(sorted_rewards))]
    return {
        "chunk_max_reward": max(rewards),
        "chunk_mean_reward": sum(rewards) / len(rewards),
        "chunk_top5_mean_reward": sum(top_values) / len(top_values),
        "zero_score_fraction": zero_score_fraction(rows),
        "unique_smiles_ratio": unique_smiles_ratio(rows),
    }
