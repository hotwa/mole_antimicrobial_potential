"""Test equivalence between old and cached assembly implementations."""
import sys
sys.path.insert(0, '.')

from rdkit import Chem
from src.stream_enumeration_screen import (
    _load_scaffolds, _load_fragment_smiles, IndexSpace,
    _assemble_molecule, _attach_fragment_cached,
    _compute_fragment_attachment_info,
    DEFAULT_POSITION_TO_ATTACHMENT
)

# Load data
scaffolds = _load_scaffolds(
    scaffold_file=None,
    scaffold_dir=None,
    scaffold_catalog=None,
)
ordinary_fragments = _load_fragment_smiles('data/05.stream_tasks/thought2_stream_screen_2026-04-23/thought2_stream_screen_2026-04-23/validation_output/thought2_enumeration/input_libraries/shared_ordinary_library.csv')
pos13_fragments = _load_fragment_smiles('data/05.stream_tasks/thought2_stream_screen_2026-04-23/thought2_stream_screen_2026-04-23/validation_output/thought2_enumeration/input_libraries/pos13_sugar_library.csv')

print(f"Scaffolds: {len(scaffolds)}")
print(f"Ordinary fragments: {len(ordinary_fragments)}")
print(f"Pos13 fragments: {len(pos13_fragments)}")

# Pre-compute scaffold mols
scaffold_mols = []
for scaffold in scaffolds:
    mol = Chem.MolFromSmiles(scaffold.scaffold_smiles)
    if mol is None:
        raise ValueError(f"Bad scaffold SMILES: {scaffold.scaffold_smiles}")
    scaffold_mols.append(mol)

# Pre-compute ordinary fragment mols and attachment info
ordinary_fragment_mols = []
ordinary_fragment_attachment_infos = []
for smiles in ordinary_fragments:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Bad fragment SMILES: {smiles}")
    ordinary_fragment_mols.append(mol)
    ordinary_fragment_attachment_infos.append(_compute_fragment_attachment_info(mol))

# Pre-compute pos13 fragment mols and attachment info
pos13_fragment_mols = []
pos13_fragment_attachment_infos = []
for smiles in pos13_fragments:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Bad pos13 fragment SMILES: {smiles}")
    pos13_fragment_mols.append(mol)
    pos13_fragment_attachment_infos.append(_compute_fragment_attachment_info(mol))

# Test equivalence for random combinations
import random
random.seed(42)

N = 1000
mismatches = []

print(f"\nTesting {N} random combinations...")

for i in range(N):
    scaffold_idx = random.randint(0, len(scaffolds) - 1)
    pos3_idx = random.randint(0, len(ordinary_fragments) - 1)
    pos4_idx = random.randint(0, len(ordinary_fragments) - 1)
    pos12_idx = random.randint(0, len(ordinary_fragments) - 1)
    pos13_idx = random.randint(0, len(pos13_fragments) - 1)

    scaffold = scaffolds[scaffold_idx]

    # Old implementation
    old_smiles = _assemble_molecule(
        scaffold.scaffold_smiles,
        {
            "pos3": ordinary_fragments[pos3_idx],
            "pos4": ordinary_fragments[pos4_idx],
            "pos12": ordinary_fragments[pos12_idx],
            "pos13": pos13_fragments[pos13_idx],
        },
    )

    # New implementation
    scaffold_mol = scaffold_mols[scaffold_idx]

    current = scaffold_mol
    for position_name, frag_idx, frag_mol, frag_info in [
        ("pos4", pos4_idx, ordinary_fragment_mols[pos4_idx], ordinary_fragment_attachment_infos[pos4_idx]),
        ("pos3", pos3_idx, ordinary_fragment_mols[pos3_idx], ordinary_fragment_attachment_infos[pos3_idx]),
        ("pos13", pos13_idx, pos13_fragment_mols[pos13_idx], pos13_fragment_attachment_infos[pos13_idx]),
        ("pos12", pos12_idx, ordinary_fragment_mols[pos12_idx], ordinary_fragment_attachment_infos[pos12_idx]),
    ]:
        attachment_label = DEFAULT_POSITION_TO_ATTACHMENT[position_name]
        current = _attach_fragment_cached(
            current, attachment_label, frag_mol, frag_info
        )

    new_smiles = Chem.MolToSmiles(current, canonical=True, isomericSmiles=True)

    if old_smiles != new_smiles:
        mismatches.append({
            "scaffold_idx": scaffold_idx,
            "pos3_idx": pos3_idx,
            "pos4_idx": pos4_idx,
            "pos12_idx": pos12_idx,
            "pos13_idx": pos13_idx,
            "old_smiles": old_smiles,
            "new_smiles": new_smiles,
        })

if mismatches:
    print(f"\nFound {len(mismatches)} mismatches!")
    for mismatch in mismatches[:10]:  # Show first 10
        print(f"  scaffold={mismatch['scaffold_idx']}, "
              f"pos3={mismatch['pos3_idx']}, pos4={mismatch['pos4_idx']}, "
              f"pos12={mismatch['pos12_idx']}, pos13={mismatch['pos13_idx']}")
        print(f"    Old: {mismatch['old_smiles']}")
        print(f"    New: {mismatch['new_smiles']}")
else:
    print(f"\nAll {N} combinations match!")

# Performance comparison
print("\n=== Performance comparison ===")
import time

N_perf = 1000

# Old implementation
start = time.perf_counter()
for i in range(N_perf):
    scaffold_idx = i % len(scaffolds)
    pos3_idx = i % len(ordinary_fragments)
    pos4_idx = (i + 1) % len(ordinary_fragments)
    pos12_idx = (i + 2) % len(ordinary_fragments)
    pos13_idx = i % len(pos13_fragments)

    scaffold = scaffolds[scaffold_idx]
    old_smiles = _assemble_molecule(
        scaffold.scaffold_smiles,
        {
            "pos3": ordinary_fragments[pos3_idx],
            "pos4": ordinary_fragments[pos4_idx],
            "pos12": ordinary_fragments[pos12_idx],
            "pos13": pos13_fragments[pos13_idx],
        },
    )
elapsed_old = time.perf_counter() - start

# New implementation
start = time.perf_counter()
for i in range(N_perf):
    scaffold_idx = i % len(scaffolds)
    pos3_idx = i % len(ordinary_fragments)
    pos4_idx = (i + 1) % len(ordinary_fragments)
    pos12_idx = (i + 2) % len(ordinary_fragments)
    pos13_idx = i % len(pos13_fragments)

    scaffold_mol = scaffold_mols[scaffold_idx]

    current = scaffold_mol
    for position_name, frag_idx, frag_mol, frag_info in [
        ("pos4", pos4_idx, ordinary_fragment_mols[pos4_idx], ordinary_fragment_attachment_infos[pos4_idx]),
        ("pos3", pos3_idx, ordinary_fragment_mols[pos3_idx], ordinary_fragment_attachment_infos[pos3_idx]),
        ("pos13", pos13_idx, pos13_fragment_mols[pos13_idx], pos13_fragment_attachment_infos[pos13_idx]),
        ("pos12", pos12_idx, ordinary_fragment_mols[pos12_idx], ordinary_fragment_attachment_infos[pos12_idx]),
    ]:
        attachment_label = DEFAULT_POSITION_TO_ATTACHMENT[position_name]
        current = _attach_fragment_cached(
            current, attachment_label, frag_mol, frag_info
        )

    new_smiles = Chem.MolToSmiles(current, canonical=True, isomericSmiles=True)
elapsed_new = time.perf_counter() - start

print(f"Old implementation: {elapsed_old/N_perf*1e6:.1f} μs/call ({N_perf/elapsed_old:.0f} calls/s)")
print(f"New implementation: {elapsed_new/N_perf*1e6:.1f} μs/call ({N_perf/elapsed_new:.0f} calls/s)")
print(f"Speedup: {elapsed_old/elapsed_new:.2f}x")
