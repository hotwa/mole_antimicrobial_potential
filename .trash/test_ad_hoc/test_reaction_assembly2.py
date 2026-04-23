"""Test ReactionFromSmarts assembly - version 2."""
from rdkit import Chem
from rdkit.Chem import rdChemReactions

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

def label_fragment(fragment_smiles: str, label: int) -> Chem.Mol:
    """Add atom map number to fragment dummy atom."""
    # Replace * with [label*]
    labeled = fragment_smiles.replace("*", f"[{label}*]", 1)
    mol = Chem.MolFromSmiles(labeled)
    if mol is None:
        raise ValueError(f"Bad fragment: {fragment_smiles} -> {labeled}")
    return mol

# Test 1: Create simpler reactions
print("=== Test 1: Create simpler reactions ===")
REACTIONS = {}
for pos, label in POSITION_LABELS.items():
    # Simpler pattern: just match dummy atoms and combine
    rxn_smarts = f"[{label}*].[{label}*]>>[{label}*]-[{label}*]"
    rxn = rdChemReactions.ReactionFromSmarts(rxn_smarts)
    if rxn is not None:
        REACTIONS[pos] = rxn
        print(f"{pos}: {rxn_smarts} -> OK")
    else:
        print(f"{pos}: {rxn_smarts} -> FAILED")

# Test 2: Run reaction with labeled fragment
print("\n=== Test 2: Run reaction with labeled fragment ===")
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

# Test 3: Try even simpler pattern
print("\n=== Test 3: Try even simpler pattern ===")
rxn_smarts3 = f"[{label}*].[{label}*]>>[{label}*]-[{label}*]"
rxn3 = rdChemReactions.ReactionFromSmarts(rxn_smarts3)
if rxn3:
    products3 = rxn3.RunReactants((scaffold, fragment))
    if products3 and products3[0]:
        product3 = products3[0][0]
        print(f"Product with pattern 3: {Chem.MolToSmiles(product3)}")
    else:
        print("Reaction 3 failed!")
else:
    print(f"Failed to create reaction: {rxn_smarts3}")

# Test 4: Try with atom map 1 and 2
print("\n=== Test 4: Try with atom map 1 and 2 ===")
rxn_smarts4 = "[*:1].[*:2]>>[*:1]-[*:2]"
rxn4 = rdChemReactions.ReactionFromSmarts(rxn_smarts4)
if rxn4:
    # Create fragment with atom map 2
    frag2 = Chem.MolFromSmiles("[2*]CC")
    print(f"Fragment2: {Chem.MolToSmiles(frag2)}")
    products4 = rxn4.RunReactants((scaffold, frag2))
    if products4 and products4[0]:
        product4 = products4[0][0]
        print(f"Product with pattern 4: {Chem.MolToSmiles(product4)}")
    else:
        print("Reaction 4 failed!")
else:
    print(f"Failed to create reaction: {rxn_smarts4}")

# Test 5: Try with atom map 1 only
print("\n=== Test 5: Try with atom map 1 only ===")
rxn_smarts5 = "[*:1].[*:1]>>[*:1]-[*:1]"
rxn5 = rdChemReactions.ReactionFromSmarts(rxn_smarts5)
if rxn5:
    # Create fragment with atom map 1
    frag1 = Chem.MolFromSmiles("[1*]CC")
    print(f"Fragment1: {Chem.MolToSmiles(frag1)}")
    products5 = rxn5.RunReactants((scaffold, frag1))
    if products5 and products5[0]:
        product5 = products5[0][0]
        print(f"Product with pattern 5: {Chem.MolToSmiles(product5)}")
    else:
        print("Reaction 5 failed!")
else:
    print(f"Failed to create reaction: {rxn_smarts5}")
