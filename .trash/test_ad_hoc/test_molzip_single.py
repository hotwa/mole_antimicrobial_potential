"""Test molzip for single position."""
from rdkit import Chem
from rdkit.Chem import rdmolops

# Test scaffold with 1 attachment point
scaffold_smiles = "[*:0]CC"

# Test fragment
fragment_smiles = "*CC"

def label_fragment(fragment_smiles: str, label: int) -> Chem.Mol:
    """Replace * with [*:label] for atom map labeling."""
    labeled = fragment_smiles.replace("*", f"[*:{label}]", 1)
    mol = Chem.MolFromSmiles(labeled)
    if mol is None:
        raise ValueError(f"Bad fragment: {fragment_smiles} -> {labeled}")
    return mol

# Test 1: Label fragment
print("=== Test 1: Label fragment ===")
labeled_fragment = label_fragment(fragment_smiles, 0)
print(f"Fragment: {fragment_smiles} -> {Chem.MolToSmiles(labeled_fragment)}")

# Test 2: Try molzip
print("\n=== Test 2: Try molzip ===")
scaffold = Chem.MolFromSmiles(scaffold_smiles)
print(f"Scaffold: {Chem.MolToSmiles(scaffold)}")

# Check atoms
print("Scaffold atoms:")
for atom in scaffold.GetAtoms():
    print(f"  Atom {atom.GetIdx()}: {atom.GetSymbol()}, map={atom.GetAtomMapNum()}, isotope={atom.GetIsotope()}")

print("Fragment atoms:")
for atom in labeled_fragment.GetAtoms():
    print(f"  Atom {atom.GetIdx()}: {atom.GetSymbol()}, map={atom.GetAtomMapNum()}, isotope={atom.GetIsotope()}")

# Try molzip
params = rdmolops.MolzipParams()
params.label = rdmolops.MolzipLabel.AtomMapNumber
try:
    result = rdmolops.molzip(scaffold, labeled_fragment, params)
    print(f"Molzip result: {Chem.MolToSmiles(result)}")
except Exception as e:
    print(f"Molzip failed: {e}")

# Test 3: Try with isotope label
print("\n=== Test 3: Try with isotope label ===")
scaffold_iso = Chem.MolFromSmiles("[1*]CC")
fragment_iso = Chem.MolFromSmiles("[1*]CC")
print(f"Scaffold (isotope): {Chem.MolToSmiles(scaffold_iso)}")
print(f"Fragment (isotope): {Chem.MolToSmiles(fragment_iso)}")

print("Scaffold atoms (isotope):")
for atom in scaffold_iso.GetAtoms():
    print(f"  Atom {atom.GetIdx()}: {atom.GetSymbol()}, map={atom.GetAtomMapNum()}, isotope={atom.GetIsotope()}")

print("Fragment atoms (isotope):")
for atom in fragment_iso.GetAtoms():
    print(f"  Atom {atom.GetIdx()}: {atom.GetSymbol()}, map={atom.GetAtomMapNum()}, isotope={atom.GetIsotope()}")

params_iso = rdmolops.MolzipParams()
params_iso.label = rdmolops.MolzipLabel.Isotope
try:
    result_iso = rdmolops.molzip(scaffold_iso, fragment_iso, params_iso)
    print(f"Molzip result (isotope): {Chem.MolToSmiles(result_iso)}")
except Exception as e:
    print(f"Molzip failed (isotope): {e}")
