"""Test equivalence between old and new assembly implementations."""
import sys
sys.path.insert(0, '.')

from rdkit import Chem
from src.stream_enumeration_screen import (
    _load_scaffolds, _load_fragment_smiles, IndexSpace,
    _assemble_molecule, _assemble_molecule_molzip,
    _convert_scaffold_to_isotope, _label_fragment_with_isotope,
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

# Pre-compute scaffold in isotope format
scaffold_mols = []
for scaffold in scaffolds:
    iso_smiles = _convert_scaffold_to_isotope(scaffold.scaffold_smiles)
    mol = Chem.MolFromSmiles(iso_smiles)
    if mol is None:
        raise ValueError(f"Bad scaffold SMILES: {scaffold.scaffold_smiles}")
    scaffold_mols.append(mol)

# Pre-compute ordinary fragments with isotope labels
ordinary_fragment_mols = []
for smiles in ordinary_fragments:
    mol_pos4 = _label_fragment_with_isotope(smiles, 1)
    mol_pos3 = _label_fragment_with_isotope(smiles, 2)
    mol_pos12 = _label_fragment_with_isotope(smiles, 4)
    ordinary_fragment_mols.append({
        "pos4": mol_pos4,
        "pos3": mol_pos3,
        "pos12": mol_pos12,
    })

# Pre-compute pos13 fragments with isotope label 3
pos13_fragment_mols = []
for smiles in pos13_fragments:
    mol = _label_fragment_with_isotope(smiles, 3)
    pos13_fragment_mols.append(mol)

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
    ordinary_frag_pos3 = ordinary_fragment_mols[pos3_idx]
    ordinary_frag_pos4 = ordinary_fragment_mols[pos4_idx]
    ordinary_frag_pos12 = ordinary_fragment_mols[pos12_idx]
    pos13_frag = pos13_fragment_mols[pos13_idx]

    fragments_by_pos = {
        "pos4": ordinary_frag_pos4["pos4"],
        "pos3": ordinary_frag_pos3["pos3"],
        "pos12": ordinary_frag_pos12["pos12"],
        "pos13": pos13_frag,
    }

    new_smiles = _assemble_molecule_molzip(scaffold_mol, fragments_by_pos)

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
    ordinary_frag_pos3 = ordinary_fragment_mols[pos3_idx]
    ordinary_frag_pos4 = ordinary_fragment_mols[pos4_idx]
    ordinary_frag_pos12 = ordinary_fragment_mols[pos12_idx]
    pos13_frag = pos13_fragment_mols[pos13_idx]

    fragments_by_pos = {
        "pos4": ordinary_frag_pos4["pos4"],
        "pos3": ordinary_frag_pos3["pos3"],
        "pos12": ordinary_frag_pos12["pos12"],
        "pos13": pos13_frag,
    }

    new_smiles = _assemble_molecule_molzip(scaffold_mol, fragments_by_pos)
elapsed_new = time.perf_counter() - start

print(f"Old implementation: {elapsed_old/N_perf*1e6:.1f} μs/call ({N_perf/elapsed_old:.0f} calls/s)")
print(f"New implementation: {elapsed_new/N_perf*1e6:.1f} μs/call ({N_perf/elapsed_new:.0f} calls/s)")
print(f"Speedup: {elapsed_old/elapsed_new:.2f}x")
