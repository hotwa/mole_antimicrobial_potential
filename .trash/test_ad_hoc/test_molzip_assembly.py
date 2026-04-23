"""Test molzip-based assembly for 4 positions."""
from rdkit import Chem
from rdkit.Chem import rdmolops

# Test scaffold with 4 attachment points
scaffold_smiles = "[*:0]CC([*:1])C([*:2])C[*:3]"

# Test fragments (one per position) - must contain * for dummy atom
fragment_smiles = {
    "pos4": "*CC",
    "pos3": "*NN",
    "pos13": "*OO",
    "pos12": "*CC",  # Use CC instead of FF to avoid valence error
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
    """Replace * with [*:label] for atom map labeling."""
    # Need to use [*:label] format, not *:label
    labeled = fragment_smiles.replace("*", f"[*:{label}]", 1)
    mol = Chem.MolFromSmiles(labeled)
    if mol is None:
        raise ValueError(f"Bad fragment: {fragment_smiles} -> {labeled}")
    return mol

def assemble_with_molzip(scaffold: Chem.Mol, fragments_by_pos: dict[str, Chem.Mol]) -> str:
    """Assemble molecule using molzip."""
    cur = scaffold
    for pos in ORDER:
        frag = fragments_by_pos[pos]
        # Create molzip params with AtomMapNumber
        params = rdmolops.MolzipParams()
        params.label = rdmolops.MolzipLabel.AtomMapNumber
        # Zip current molecule with fragment
        cur = rdmolops.molzip(cur, frag, params)
    return Chem.MolToSmiles(cur, canonical=True, isomericSmiles=True)

# Test 1: Label fragments
print("=== Test 1: Label fragments ===")
labeled_fragments = {}
for pos, smiles in fragment_smiles.items():
    label = POSITION_LABELS[pos]
    mol = label_fragment(smiles, label)
    labeled_fragments[pos] = mol
    print(f"{pos}: {smiles} -> {Chem.MolToSmiles(mol)}")

# Test 2: Assemble molecule
print("\n=== Test 2: Assemble molecule ===")
scaffold = Chem.MolFromSmiles(scaffold_smiles)
print(f"Scaffold: {Chem.MolToSmiles(scaffold)}")

result_smiles = assemble_with_molzip(scaffold, labeled_fragments)
print(f"Result: {result_smiles}")

# Test 3: Compare with current implementation
print("\n=== Test 3: Compare with current implementation ===")
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

# Run current implementation
current = scaffold
for pos in ORDER:
    label = POSITION_LABELS[pos]
    current = _attach_fragment(current, label, fragment_smiles[pos])
current_smiles = Chem.MolToSmiles(current, canonical=True, isomericSmiles=True)

print(f"Current implementation: {current_smiles}")
print(f"Molzip implementation:  {result_smiles}")
print(f"Match: {current_smiles == result_smiles}")

# Test 4: Performance comparison
print("\n=== Test 4: Performance comparison ===")
import time

N = 1000

# Current implementation
start = time.perf_counter()
for _ in range(N):
    current = scaffold
    for pos in ORDER:
        label = POSITION_LABELS[pos]
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
