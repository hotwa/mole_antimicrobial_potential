"""Debug valence error."""
import sys
sys.path.insert(0, '.')

from rdkit import Chem
from src.stream_enumeration_screen import (
    _load_scaffolds, _load_fragment_smiles,
    _assemble_molecule, _assemble_molecule_molzip,
    _convert_scaffold_to_isotope, _label_fragment_with_isotope,
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

# Test specific combination that causes valence error
# From the error, it seems to be scaffold_idx=0, some fragment combination
# Let's test a few combinations to find the problematic one
import random
random.seed(42)

print("\nTesting combinations to find valence error...")

for i in range(100):
    scaffold_idx = 0  # Only 1 scaffold
    pos3_idx = random.randint(0, len(ordinary_fragments) - 1)
    pos4_idx = random.randint(0, len(ordinary_fragments) - 1)
    pos12_idx = random.randint(0, len(ordinary_fragments) - 1)
    pos13_idx = random.randint(0, len(pos13_fragments) - 1)

    scaffold = scaffolds[scaffold_idx]

    try:
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
            print(f"Mismatch at combination {i}:")
            print(f"  pos3={pos3_idx}, pos4={pos4_idx}, pos12={pos12_idx}, pos13={pos13_idx}")
            print(f"  Old: {old_smiles}")
            print(f"  New: {new_smiles}")
    except Exception as e:
        print(f"Error at combination {i}:")
        print(f"  pos3={pos3_idx}, pos4={pos4_idx}, pos12={pos12_idx}, pos13={pos13_idx}")
        print(f"  Error: {e}")
        print(f"  Fragments:")
        print(f"    pos3: {ordinary_fragments[pos3_idx]}")
        print(f"    pos4: {ordinary_fragments[pos4_idx]}")
        print(f"    pos12: {ordinary_fragments[pos12_idx]}")
        print(f"    pos13: {pos13_fragments[pos13_idx]}")
        break
