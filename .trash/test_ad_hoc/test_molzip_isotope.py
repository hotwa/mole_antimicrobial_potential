"""Test molzip with isotope label for 4 positions."""
from rdkit import Chem
from rdkit.Chem import rdmolops

# Test scaffold with 4 attachment points
scaffold_smiles = "[*:0]CC([*:1])C([*:2])C[*:3]"

# Test fragments (one per position) - must contain * for dummy atom
fragment_smiles = {
    "pos4": "*CC",
    "pos3": "*NN",
    "pos13": "*OO",
    "pos12": "*CC",
}

# Position labels (use isotope numbers)
POSITION_LABELS = {
    "pos4": 1,
    "pos3": 2,
    "pos13": 3,
    "pos12": 4,
}

ORDER = ("pos4", "pos3", "pos13", "pos12")

def convert_scaffold_to_isotope(scaffold_smiles: str) -> str:
    """Convert scaffold from atom map to isotope format."""
    # Replace [*:0] with [1*], [*:1] with [2*], etc.
    result = scaffold_smiles
    for i in range(4):
        result = result.replace(f"[*:{i}]", f"[{i+1}*]")
    return result

def label_fragment_with_isotope(fragment_smiles: str, isotope: int) -> Chem.Mol:
    """Replace * with [isotope*] for isotope labeling."""
    labeled = fragment_smiles.replace("*", f"[{isotope}*]", 1)
    mol = Chem.MolFromSmiles(labeled)
    if mol is None:
        raise ValueError(f"Bad fragment: {fragment_smiles} -> {labeled}")
    return mol

def assemble_with_molzip(scaffold: Chem.Mol, fragments_by_pos: dict[str, Chem.Mol]) -> str:
    """Assemble molecule using molzip with isotope label."""
    cur = scaffold
    for pos in ORDER:
        frag = fragments_by_pos[pos]
        # Create molzip params with Isotope
        params = rdmolops.MolzipParams()
        params.label = rdmolops.MolzipLabel.Isotope
        # Zip current molecule with fragment
        cur = rdmolops.molzip(cur, frag, params)
    return Chem.MolToSmiles(cur, canonical=True, isomericSmiles=True)

# Test 1: Convert scaffold to isotope format
print("=== Test 1: Convert scaffold to isotope format ===")
scaffold_iso_smiles = convert_scaffold_to_isotope(scaffold_smiles)
print(f"Original scaffold: {scaffold_smiles}")
print(f"Isotope scaffold:  {scaffold_iso_smiles}")

# Test 2: Label fragments with isotope
print("\n=== Test 2: Label fragments with isotope ===")
labeled_fragments = {}
for pos, smiles in fragment_smiles.items():
    isotope = POSITION_LABELS[pos]
    mol = label_fragment_with_isotope(smiles, isotope)
    labeled_fragments[pos] = mol
    print(f"{pos}: {smiles} -> {Chem.MolToSmiles(mol)}")

# Test 3: Assemble molecule
print("\n=== Test 3: Assemble molecule ===")
scaffold = Chem.MolFromSmiles(scaffold_iso_smiles)
print(f"Scaffold: {Chem.MolToSmiles(scaffold)}")

result_smiles = assemble_with_molzip(scaffold, labeled_fragments)
print(f"Result: {result_smiles}")

# Test 4: Compare with current implementation
print("\n=== Test 4: Compare with current implementation ===")
# Simulate current implementation
def _find_scaffold_attachment_index(mol, attachment_label):
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() == attachment_label:
            return atom.GetIdx()
    raise ValueError(f"Missing scaffold attachment point [*:{attachment_label}]")

def _find_fragment_attachment_index(mol):
    matches = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]
    if len(matches) != 1:
        raise ValueError("Fragments must contain exactly one dummy atom")
    return matches[0]

def _attach_fragment(scaffold, attachment_label, fragment_smiles):
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

# Run current implementation (using original scaffold with atom map)
scaffold_original = Chem.MolFromSmiles(scaffold_smiles)
current = scaffold_original
for pos in ORDER:
    label = POSITION_LABELS[pos] - 1  # Convert to 0-based
    current = _attach_fragment(current, label, fragment_smiles[pos])
current_smiles = Chem.MolToSmiles(current, canonical=True, isomericSmiles=True)

print(f"Current implementation: {current_smiles}")
print(f"Molzip implementation:  {result_smiles}")
print(f"Match: {current_smiles == result_smiles}")

# Test 5: Performance comparison
print("\n=== Test 5: Performance comparison ===")
import time

N = 1000

# Current implementation
start = time.perf_counter()
for _ in range(N):
    current = scaffold_original
    for pos in ORDER:
        label = POSITION_LABELS[pos] - 1
        current = _attach_fragment(current, label, fragment_smiles[pos])
    result = Chem.MolToSmiles(current, canonical=True, isomericSmiles=True)
elapsed_current = time.perf_counter() - start

# Molzip implementation
start = time.perf_counter()
for _ in range(N):
    result = assemble_with_molzip(scaffold, labeled_fragments)
elapsed_molzip = time.perf_counter() - start

print(f"Current implementation: {elapsed_current/N*1e6:.1f} μs/call ({N/elapsed_current:.0f} calls/s)")
print(f"Molzip implementation:  {elapsed_molzip/N*1e6:.1f} μs/call ({N/elapsed_molzip:.0f} calls/s)")
print(f"Speedup: {elapsed_current/elapsed_molzip:.2f}x")
