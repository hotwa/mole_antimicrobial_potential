"""Test molzip approach for fragment attachment."""
from rdkit import Chem
from rdkit.Chem import rdmolops

# Test scaffold with 4 attachment points
scaffold_smiles = "[*:0]CC([*:1])C([*:2])C[*:3]"

# Test fragments
fragment_smiles = {
    "pos4": "CC",
    "pos3": "NN",
    "pos13": "OO",
    "pos12": "FF",
}

# Position labels
POSITION_LABELS = {
    "pos4": 0,
    "pos3": 1,
    "pos13": 2,
    "pos12": 3,
}

ORDER = ("pos4", "pos3", "pos13", "pos12")

def label_fragment_with_isotope(fragment_smiles: str, label: int) -> Chem.Mol:
    """Replace * with [label*] for isotope labeling."""
    labeled = fragment_smiles.replace("*", f"[{label}*]", 1)
    mol = Chem.MolFromSmiles(labeled)
    if mol is None:
        raise ValueError(f"Bad fragment: {fragment_smiles} -> {labeled}")
    return mol

def label_fragment_with_atommap(fragment_smiles: str, label: int) -> Chem.Mol:
    """Replace * with [:label] for atom map labeling."""
    labeled = fragment_smiles.replace("*", f"[:{label}]", 1)
    mol = Chem.MolFromSmiles(labeled)
    if mol is None:
        raise ValueError(f"Bad fragment: {fragment_smiles} -> {labeled}")
    return mol

# Test 1: Try molzip with isotope label
print("=== Test 1: Try molzip with isotope label ===")
scaffold = Chem.MolFromSmiles(scaffold_smiles)
print(f"Scaffold: {Chem.MolToSmiles(scaffold)}")

# Create fragment with isotope label
pos = "pos4"
label = POSITION_LABELS[pos]
fragment = label_fragment_with_isotope(fragment_smiles[pos], label)
print(f"Fragment ({pos}): {Chem.MolToSmiles(fragment)}")

# Try molzip
params = rdmolops.MolzipParams()
params.label = rdmolops.MolzipLabel.Isotope
try:
    result = rdmolops.molzip(scaffold, fragment, params)
    print(f"Molzip result: {Chem.MolToSmiles(result)}")
except Exception as e:
    print(f"Molzip failed: {e}")

# Test 2: Try molzip with atom map label
print("\n=== Test 2: Try molzip with atom map label ===")
fragment2 = label_fragment_with_atommap(fragment_smiles[pos], label)
print(f"Fragment2 ({pos}): {Chem.MolToSmiles(fragment2)}")

params2 = rdmolops.MolzipParams()
params2.label = rdmolops.MolzipLabel.AtomMap
try:
    result2 = rdmolops.molzip(scaffold, fragment2, params2)
    print(f"Molzip result2: {Chem.MolToSmiles(result2)}")
except Exception as e:
    print(f"Molzip failed: {e}")

# Test 3: Try with a simpler scaffold
print("\n=== Test 3: Try with a simpler scaffold ===")
scaffold_simple = Chem.MolFromSmiles("[1*]CC[2*]")
print(f"Simple scaffold: {Chem.MolToSmiles(scaffold_simple)}")

fragment_simple = Chem.MolFromSmiles("[1*]NN")
print(f"Simple fragment: {Chem.MolToSmiles(fragment_simple)}")

params3 = rdmolops.MolzipParams()
params3.label = rdmolops.MolzipLabel.Isotope
try:
    result3 = rdmolops.molzip(scaffold_simple, fragment_simple, params3)
    print(f"Molzip result3: {Chem.MolToSmiles(result3)}")
except Exception as e:
    print(f"Molzip failed: {e}")

# Test 4: Try with atom map
print("\n=== Test 4: Try with atom map ===")
scaffold_map = Chem.MolFromSmiles("[:1]CC[:2]")
print(f"Scaffold with map: {Chem.MolToSmiles(scaffold_map)}")

fragment_map = Chem.MolFromSmiles("[:1]NN")
print(f"Fragment with map: {Chem.MolToSmiles(fragment_map)}")

params4 = rdmolops.MolzipParams()
params4.label = rdmolops.MolzipLabel.AtomMap
try:
    result4 = rdmolops.molzip(scaffold_map, fragment_map, params4)
    print(f"Molzip result4: {Chem.MolToSmiles(result4)}")
except Exception as e:
    print(f"Molzip failed: {e}")
