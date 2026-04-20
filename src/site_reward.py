"""Site-aware structural reward helpers for REINVENT4 experiments."""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence

try:
    from rdkit import Chem
    from rdkit import RDLogger
    from rdkit.Chem import rdRGroupDecomposition

    RDLogger.DisableLog("rdApp.*")
    RDKIT_AVAILABLE = True
except ImportError:  # pragma: no cover - optional dependency
    RDKIT_AVAILABLE = False


def _sigmoid(value: float) -> float:
    return 1.0 / (1.0 + math.exp(-value))


def _heavy_atom_count(fragment_smiles: str) -> int | None:
    if not RDKIT_AVAILABLE:
        raise RuntimeError("RDKit is required for site_reward calculations")
    mol = Chem.MolFromSmiles(fragment_smiles)
    if mol is None:
        return None
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)


def _std(values: Sequence[float]) -> float:
    if not values:
        return 0.0
    mean_value = sum(values) / len(values)
    variance = sum((value - mean_value) ** 2 for value in values) / len(values)
    return math.sqrt(variance)


def normalize_scaffold_for_rgroup(scaffold_smiles: str) -> str:
    """Normalize a labeled scaffold into the unlabeled canonical form used for decomposition."""
    if not RDKIT_AVAILABLE:
        raise RuntimeError("RDKit is required for site_reward calculations")
    mol = Chem.MolFromSmiles(scaffold_smiles)
    if mol is None:
        raise ValueError(f"Invalid scaffold_smiles for site_reward: {scaffold_smiles}")
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def scaffold_from_file(scaffold_path: str | Path) -> str:
    """Read the first non-comment scaffold from a .smi file."""
    path = Path(scaffold_path).expanduser().resolve()
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        scaffold = raw_line.strip()
        if scaffold and not scaffold.startswith("#"):
            return scaffold
    raise ValueError(f"No scaffold entries found in {path}")


def attachment_count(scaffold_smiles: str) -> int:
    """Count attachment points in a labeled scaffold."""
    if not RDKIT_AVAILABLE:
        raise RuntimeError("RDKit is required for site_reward calculations")
    mol = Chem.MolFromSmiles(scaffold_smiles)
    if mol is None:
        raise ValueError(f"Invalid scaffold_smiles for site_reward: {scaffold_smiles}")
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0)


def decompose_rgroups(smiles: str, scaffold_smiles: str) -> Dict[str, Any]:
    """Decompose a molecule into R-groups against a fixed scaffold."""
    if not RDKIT_AVAILABLE:
        raise RuntimeError("RDKit is required for site_reward calculations")

    normalized_scaffold = normalize_scaffold_for_rgroup(scaffold_smiles)
    scaffold = Chem.MolFromSmiles(normalized_scaffold)
    molecule = Chem.MolFromSmiles(smiles)
    if scaffold is None:
        raise ValueError(f"Invalid normalized scaffold for site_reward: {normalized_scaffold}")
    if molecule is None:
        raise ValueError(f"Invalid molecule SMILES for site_reward: {smiles}")

    params = rdRGroupDecomposition.RGroupDecompositionParameters()
    params.removeAllHydrogenRGroups = False
    params.removeAllHydrogenRGroupsAndLabels = False
    decomposition = rdRGroupDecomposition.RGroupDecomposition([scaffold], params)
    added_index = decomposition.Add(molecule)
    if added_index == -1 or not decomposition.Process():
        return {
            "success": False,
            "normalized_scaffold_smiles": normalized_scaffold,
            "attachment_count": attachment_count(scaffold_smiles),
            "rgroup_smiles": [],
            "site_heavy_atoms": [],
            "error": "rgroup_decomposition_failed",
        }

    rows = decomposition.GetRGroupsAsRows(asSmiles=True)
    if not rows:
        return {
            "success": False,
            "normalized_scaffold_smiles": normalized_scaffold,
            "attachment_count": attachment_count(scaffold_smiles),
            "rgroup_smiles": [],
            "site_heavy_atoms": [],
            "error": "rgroup_rows_empty",
        }

    row = rows[0]
    fragment_keys = sorted(
        (key for key in row if key.startswith("R")),
        key=lambda item: int(item[1:]),
    )
    fragments = [str(row[key]) for key in fragment_keys]
    heavy_atoms = [_heavy_atom_count(fragment) for fragment in fragments]
    return {
        "success": True,
        "normalized_scaffold_smiles": normalized_scaffold,
        "attachment_count": attachment_count(scaffold_smiles),
        "rgroup_smiles": fragments,
        "site_heavy_atoms": heavy_atoms,
        "error": None,
    }


def compute_site_reward(
    smiles: str,
    config: Mapping[str, Any],
) -> Dict[str, Any]:
    """Compute the experimental site_reward for a molecule."""
    lower = int(config.get("range_min", 4))
    upper = int(config.get("range_max", 27))
    alpha = float(config.get("alpha", 1.5))
    beta = float(config.get("beta", 2.5))
    coverage_weight = float(config.get("coverage_weight", 0.7))
    balance_weight = float(config.get("balance_weight", 0.3))

    decomposition = decompose_rgroups(smiles, str(config["scaffold_smiles"]))
    site_heavy_atoms: List[int | None] = list(decomposition["site_heavy_atoms"])
    expected_sites = int(decomposition["attachment_count"])
    while len(site_heavy_atoms) < expected_sites:
        site_heavy_atoms.append(None)

    per_site_scores: List[float] = []
    clipped_values: List[float] = []
    for count in site_heavy_atoms[:expected_sites]:
        if count is None:
            per_site_scores.append(0.0)
            clipped_values.append(float(lower))
            continue
        per_site_scores.append(
            _sigmoid((count - lower) / alpha) * _sigmoid((upper - count) / beta)
        )
        clipped_values.append(float(min(max(count, lower), upper)))

    coverage = sum(per_site_scores) / len(per_site_scores) if per_site_scores else 0.0
    if clipped_values:
        clipped_mean = sum(clipped_values) / len(clipped_values)
        balance = math.exp(-_std(clipped_values) / (clipped_mean + 1e-6))
    else:
        balance = 0.0
    site_reward = (coverage_weight * coverage) + (balance_weight * balance)

    return {
        "site_reward": site_reward,
        "site_coverage": coverage,
        "site_balance": balance,
        "site_heavy_atoms": site_heavy_atoms[:expected_sites],
        "site_rgroups": list(decomposition["rgroup_smiles"]),
        "site_decomposition_success": bool(decomposition["success"]),
        "site_decomposition_error": decomposition["error"],
        "site_reward_parameters": {
            "range_min": lower,
            "range_max": upper,
            "alpha": alpha,
            "beta": beta,
            "coverage_weight": coverage_weight,
            "balance_weight": balance_weight,
            "lambda": float(config.get("lambda", 0.85)),
            "normalized_scaffold_smiles": decomposition["normalized_scaffold_smiles"],
        },
    }
