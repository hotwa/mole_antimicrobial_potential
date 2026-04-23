"""Test parallel optimization strategies - simple version."""
import sys
sys.path.insert(0, '.')

import time
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
from src.stream_enumeration_screen import (
    _load_scaffolds, _load_fragment_smiles, IndexSpace,
    _init_materialization_worker, _materialize_batch_in_worker
)

if __name__ == '__main__':
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

    # Create IndexSpace
    space = IndexSpace(
        scaffold_count=len(scaffolds),
        pos3_count=len(ordinary_fragments),
        pos4_count=len(ordinary_fragments),
        pos12_count=len(ordinary_fragments),
        pos13_count=len(pos13_fragments),
    )

    # Test different batch sizes
    N = 10000
    batch_sizes = [256, 512, 1000]
    workers_list = [1, 2, 4]

    print(f"\nTesting {N} molecules with different configurations...")

    for workers in workers_list:
        for batch_size in batch_sizes:
            if batch_size > N:
                continue

            # Create batches
            batches = []
            for start in range(0, N, batch_size):
                end = min(N, start + batch_size)
                batches.append((start, end))

            # Test with ProcessPoolExecutor
            start_time = time.perf_counter()
            with ProcessPoolExecutor(
                max_workers=workers,
                mp_context=mp.get_context("spawn"),
                initializer=_init_materialization_worker,
                initargs=(space, tuple(scaffolds), tuple(ordinary_fragments), tuple(pos13_fragments)),
            ) as executor:
                # Submit all batches
                futures = [executor.submit(_materialize_batch_in_worker, start, end) for start, end in batches]
                # Wait for all results
                results = [f.result() for f in futures]
                total_mol = sum(len(r) for r in results)
            elapsed = time.perf_counter() - start_time

            print(f"workers={workers}, batch_size={batch_size}: {elapsed:.3f}s for {total_mol} mol ({total_mol/elapsed:.0f} mol/s)")
