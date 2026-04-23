"""Profile RDKit materialization operations."""
import sys
sys.path.insert(0, '.')

import time
from rdkit import Chem
from src.stream_enumeration_screen import (
    _load_scaffolds, _load_fragment_smiles, IndexSpace, _materialize_batch,
    _attach_fragment, _assemble_molecule, DEFAULT_POSITION_TO_ATTACHMENT
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

# Profile individual operations
N = 1000

# 1. MolFromSmiles
start = time.perf_counter()
for _ in range(N):
    mol = Chem.MolFromSmiles("CCO")
elapsed = time.perf_counter() - start
print(f"\nChem.MolFromSmiles('CCO'): {elapsed/N*1e6:.1f} μs/call ({N/elapsed:.0f} calls/s)")

# 2. MolToSmiles
mol = Chem.MolFromSmiles("CCO")
start = time.perf_counter()
for _ in range(N):
    smi = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
elapsed = time.perf_counter() - start
print(f"Chem.MolToSmiles(canonical=True): {elapsed/N*1e6:.1f} μs/call ({N/elapsed:.0f} calls/s)")

# 3. MolToSmiles without canonical
start = time.perf_counter()
for _ in range(N):
    smi = Chem.MolToSmiles(mol, canonical=False, isomericSmiles=False)
elapsed = time.perf_counter() - start
print(f"Chem.MolToSmiles(canonical=False): {elapsed/N*1e6:.1f} μs/call ({N/elapsed:.0f} calls/s)")

# 4. _attach_fragment
scaffold_mol = Chem.MolFromSmiles(scaffolds[0].scaffold_smiles)
fragment_smiles = ordinary_fragments[0]
start = time.perf_counter()
for _ in range(N):
    result = _attach_fragment(scaffold_mol, DEFAULT_POSITION_TO_ATTACHMENT["pos3"], fragment_smiles)
elapsed = time.perf_counter() - start
print(f"_attach_fragment: {elapsed/N*1e6:.1f} μs/call ({N/elapsed:.0f} calls/s)")

# 5. Full _assemble_molecule
start = time.perf_counter()
for i in range(N):
    smi = _assemble_molecule(
        scaffolds[0].scaffold_smiles,
        {
            "pos3": ordinary_fragments[i % len(ordinary_fragments)],
            "pos4": ordinary_fragments[(i+1) % len(ordinary_fragments)],
            "pos12": ordinary_fragments[(i+2) % len(ordinary_fragments)],
            "pos13": pos13_fragments[i % len(pos13_fragments)],
        },
    )
elapsed = time.perf_counter() - start
print(f"\n_assemble_molecule (full): {elapsed/N*1e6:.1f} μs/call ({N/elapsed:.0f} calls/s)")

# 6. Full _materialize_batch
space = IndexSpace(
    scaffold_count=len(scaffolds),
    pos3_count=len(ordinary_fragments),
    pos4_count=len(ordinary_fragments),
    pos12_count=len(ordinary_fragments),
    pos13_count=len(pos13_fragments),
)
start = time.perf_counter()
rows = _materialize_batch(
    start_idx=0,
    end_idx=N,
    space=space,
    scaffolds=scaffolds,
    ordinary_fragments=ordinary_fragments,
    pos13_fragments=pos13_fragments,
)
elapsed = time.perf_counter() - start
print(f"_materialize_batch ({N} mol): {elapsed/N*1e6:.1f} μs/call ({N/elapsed:.0f} calls/s)")

# 7. Parallel test
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

def _worker_batch(args):
    start_idx, end_idx, space, scaffolds, ordinary_fragments, pos13_fragments = args
    return _materialize_batch(
        start_idx=start_idx,
        end_idx=end_idx,
        space=space,
        scaffolds=scaffolds,
        ordinary_fragments=ordinary_fragments,
        pos13_fragments=pos13_fragments,
    )

for workers in [1, 2, 4, 8]:
    batch_size = N // workers
    args_list = [
        (i * batch_size, (i+1) * batch_size, space, scaffolds, ordinary_fragments, pos13_fragments)
        for i in range(workers)
    ]
    start = time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as executor:
        results = list(executor.map(_worker_batch, args_list))
    elapsed = time.perf_counter() - start
    total_mol = sum(len(r) for r in results)
    print(f"ProcessPoolExecutor(workers={workers}): {elapsed:.3f}s for {total_mol} mol ({total_mol/elapsed:.0f} mol/s)")
