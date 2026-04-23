"""Benchmark materialization speed."""
import sys
sys.path.insert(0, '.')

import time
from src.stream_enumeration_screen import (
    _load_scaffolds, _load_fragment_smiles, IndexSpace, _materialize_batch
)

# Load scaffold and fragments
scaffolds = _load_scaffolds(
    scaffold_file=None,
    scaffold_dir=None,
    scaffold_catalog=None,
)
ordinary_fragments = _load_fragment_smiles('data/05.stream_tasks/thought2_stream_screen_2026-04-23/thought2_stream_screen_2026-04-23/validation_output/thought2_enumeration/input_libraries/shared_ordinary_library.csv')
pos13_fragments = _load_fragment_smiles('data/05.stream_tasks/thought2_stream_screen_2026-04-23/thought2_stream_screen_2026-04-23/validation_output/thought2_enumeration/input_libraries/pos13_sugar_library.csv')

space = IndexSpace(
    scaffold_count=len(scaffolds),
    pos3_count=len(ordinary_fragments),
    pos4_count=len(ordinary_fragments),
    pos12_count=len(ordinary_fragments),
    pos13_count=len(pos13_fragments),
)

print(f'Space: {space}')
print(f'Total combinations: {space.total_combinations:,}')

# Time materialization at different batch sizes
for batch_size in [100, 500, 1000]:
    start = time.perf_counter()
    rows = _materialize_batch(
        start_idx=0,
        end_idx=batch_size,
        space=space,
        scaffolds=scaffolds,
        ordinary_fragments=ordinary_fragments,
        pos13_fragments=pos13_fragments,
    )
    elapsed = time.perf_counter() - start
    print(f'Materialized {batch_size} molecules in {elapsed:.3f}s ({batch_size/elapsed:.0f} mol/s)')
