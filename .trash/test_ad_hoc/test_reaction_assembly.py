"""Test ReactionFromSmarts assembly approach."""
from rdkit import Chem
from rdkit.Chem import rdChemReactions

# Test scaffold with 4 attachment points
scaffold_smiles = "[*:0]CC([*:1])C([*:2])C[*:3]"

# Test fragments
fragment_smiles = {
    "pos4": "CC",  # will be labeled as [0*]CC
    "pos3": "NN",  # will be labeled as [1*]NN
    "pos13": "OO",  # will be labeled as [2*]OO
    "pos12": "FF",  # will be labeled as [3*]FF
}

# Position labels
POSITION_LABELS = {
    "pos4": 0,
    "pos3": 1,
    "pos13": 2,
    "pos12": 3,
}

ORDER = ("pos4", "pos3", "pos13", "pos12")

def label_fragment(fragment_smiles: str, label: int) -> Chem.Mol:
    """Add atom map number to fragment dummy atom."""
    # Replace * with [label*]
    labeled = fragment_smiles.replace("*", f"[{label}*]", 1)
    mol = Chem.MolFromSmiles(labeled)
    if mol is None:
        raise ValueError(f"Bad fragment: {fragment_smiles} -> {labeled}")
    return mol

# Test 1: Label fragments
print("=== Test 1: Label fragments ===")
for pos, smiles in fragment_smiles.items():
    label = POSITION_LABELS[pos]
    mol = label_fragment(smiles, label)
    print(f"{pos}: {smiles} -> {Chem.MolToSmiles(mol)}")

# Test 2: Create reactions
print("\n=== Test 2: Create reactions ===")
REACTIONS = {}
for pos, label in POSITION_LABELS.items():
    # Try different SMARTS patterns
    # Pattern 1: [*:1]-[label*].[label*]-[*:2]>>[*:1]-[*:2]
    rxn_smarts = f"[*:1]-[{label}*].[{label}*]-[*:2]>>[*:1]-[*:2]"
    rxn = rdChemReactions.ReactionFromSmarts(rxn_smarts)
    if rxn is not None:
        REACTIONS[pos] = rxn
        print(f"{pos}: {rxn_smarts} -> OK")
    else:
        print(f"{pos}: {rxn_smarts} -> FAILED")

# Test 3: Run reaction
print("\n=== Test 3: Run reaction ===")
scaffold = Chem.MolFromSmiles(scaffold_smiles)
print(f"Scaffold: {Chem.MolToSmiles(scaffold)}")

# Try pos4 reaction
pos = "pos4"
label = POSITION_LABELS[pos]
fragment = label_fragment(fragment_smiles[pos], label)
print(f"Fragment ({pos}): {Chem.MolToSmiles(fragment)}")

rxn = REACTIONS[pos]
products = rxn.RunReactants((scaffold, fragment))
if products and products[0]:
    product = products[0][0]
    print(f"Product: {Chem.MolToSmiles(product)}")
else:
    print("Reaction failed!")

# Test 4: Try different SMARTS pattern
print("\n=== Test 4: Try different SMARTS pattern ===")
rxn_smarts2 = f"[*:1][{label}*].[{label}*][*:2]>>[*:1]-[*:2]"
rxn2 = rdChemReactions.ReactionFromSmarts(rxn_smarts2)
if rxn2:
    products2 = rxn2.RunReactants((scaffold, fragment))
    if products2 and products2[0]:
        product2 = products2[0][0]
        print(f"Product with pattern 2: {Chem.MolToSmiles(product2)}")
    else:
        print("Reaction 2 failed!")
else:
    print(f"Failed to create reaction: {rxn_smarts2}")

# Test 5: Try molzip
print("\n=== Test 5: Try molzip ===")
from rdkit.Chem import rdmolops
try:
    # Create labeled fragment for molzip
    frag_labeled = Chem.MolFromSmiles(f"[{label}*]CC")
    print(f"Fragment for molzip: {Chem.MolToSmiles(frag_labeled)}")

    # Try molzip with isotope label
    params = rdmolops.MolzipParams()
    params.label = rdmolops.MolzipLabel.Isotope
    result = rdmolops.molzip(scaffold, frag_labeled, params)
    print(f"Molzip result: {Chem.MolToSmiles(result)}")
except Exception as e:
    print(f"Molzip failed: {e}")
