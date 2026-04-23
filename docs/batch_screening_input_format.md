# Screening Input Format Guide

This is the canonical handbook for high-throughput `mole screen` runs.

For the full multi-command CLI parameter reference, see
[docs/cli_reference.md](/home/lingyuzeng/project/mole_antimicrobial_potential/docs/cli_reference.md).

Use it for:

- supported screening input formats
- preprocessing large CSV/TSV inputs into Parquet shards
- `mole screen` parameter meanings and defaults
- deciding which knobs change result content vs. which only change
  performance, scheduling, or reproducibility

## Preferred Execution Format

`mole screen` supports these input types:

- CSV / TSV tables
- Parquet files
- Parquet shard directories
- `tar` / `tar.gz` archives with CSV members
- SQLite database files

For throughput-oriented runs, prefer Parquet files or Parquet shard
directories.

Raw `tar.gz` bundles and one huge CSV are transport/storage formats. They are
not the preferred execution format because they delay first-batch prediction and
increase repeated parsing cost.

## Why Preprocess Into Parquet Shards

Preprocessing helps because it:

- keeps only the columns used by screening: `smiles`, `chem_id`,
  `source_group`
- avoids repeated text parsing of giant CSV files
- lets the screening loop start consuming ready chunks earlier
- gives stable shard sizes and explicit Parquet row groups
- makes multi-process and multi-GPU sharding straightforward

If you want a lightweight comparison before preprocessing, run:

```bash
pixi run mole benchmark-screening-inputs \
  --input-path <raw_csv_or_tsv> \
  --input-path <prepared_parquet_or_dir>
```

This command is diagnostic only. It does not run MolE prediction, but it is
useful for checking candidate input paths before a full screening job.

Recommended shard properties:

- rows per shard: `50,000` to `200,000`
- row group size: `2,000` to `8,000`
- shard grouping: one natural `source_group` per scaffold, batch, or source

## `mole preprocess-screening-input`

Use this command to convert large CSV/TSV inputs into Parquet shards:

```bash
pixi run mole preprocess-screening-input \
  --input-path <raw_csv_or_tsv> \
  --output-dir <prepared_dir> \
  --smiles-colname <smiles_col> \
  --chem-id-colname <chem_id_col> \
  --source-group <group_name> \
  --rows-per-shard 100000 \
  --row-group-size 4000
```

### Parameters

| Parameter | Default | Effect type | Meaning |
| --- | --- | --- | --- |
| `--input-path` | required | content | Input CSV/TSV path |
| `--output-dir` | required | output only | Root directory for written shards |
| `--smiles-colname` | `smiles` | content | Source SMILES column name |
| `--chem-id-colname` | `chem_id` | content | Source chem_id column name; missing values are auto-filled |
| `--source-group` | `input` | content | Logical label written into shard rows |
| `--rows-per-shard` | `100000` | performance/layout | Rows written per Parquet shard |
| `--row-group-size` | `4000` | performance/layout | Parquet row-group size |
| `--output` | unset | output only | Optional manifest JSON path |

### Output Schema

Written Parquet shards contain:

- `smiles`
- `chem_id`
- `source_group`

If `chem_id` is missing, stable IDs are generated in the form
`<source_group>__<row_index>`.

## `mole screen`

Screen prepared Parquet shards or other supported inputs:

```bash
pixi run mole screen \
  --input-path <prepared_parquet_or_dir> \
  --output-dir <screen_run_dir> \
  --classifier-backend pickle \
  --grouping-mode auto \
  --cpu-workers auto \
  --input-chunk-size 4000 \
  --max-batch-size 16384 \
  --target-gpu-memory-fraction 0.8 \
  --prediction-row-budget 8192 \
  --num-graph-workers 0 \
  --graph-batch-size 1024 \
  --profiling
```

### Input Formats

`--input-path` may point to:

- one CSV or TSV file
- one Parquet file
- one directory containing Parquet shards
- one `tar` or `tar.gz` archive
- one SQLite database file

### Parameters That Change Result Content

These knobs can change which rows are loaded, how duplicates are handled, or
what output structure is emitted:

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--smiles-colname` | `smiles` | SMILES column for tabular inputs |
| `--chem-id-colname` | `chem_id` | chem_id column for tabular inputs |
| `--archive-pattern` | `*_scheme_b_unique_products.csv` | Which CSV members inside tar archives are read |
| `--archive-smiles-colname` | `product_smiles_canonical` | SMILES column inside archive CSV members |
| `--archive-chem-id-colname` | `example_combo_id` | chem_id column inside archive CSV members |
| `--sqlite-table` | unset | Explicit SQLite table selection |
| `--sqlite-query` | unset | Explicit SQLite query selection |
| `--no-dedupe-smiles` | dedupe enabled | Keep duplicate SMILES instead of collapsing to first occurrence |
| `--aggregate-scores` | enabled | Emit one aggregate row per molecule |
| `--per-strain` | disabled | Emit per-strain probabilities instead of aggregate rows |
| `--app-threshold` | `0.04374140128493309` | Growth inhibition threshold |
| `--min-nkill` | `10` | Broad-spectrum threshold |
| `--classifier-backend` | `auto` | Backend selector: `auto`, `pickle`, `timber` |

Notes:

- `auto` currently prefers the original pickle/XGBoost model when available.
- `pickle` is the closest path to the original MolE prediction workflow.
- `timber` is a portability/performance option and may not be numerically
  identical to pickle.

### Parameters That Primarily Affect Performance or Scheduling

These knobs are intended to change throughput, queueing, memory pressure, or
device utilization, not the semantic prediction path for a fixed model/backend:

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--grouping-mode` | `auto` | Work-unit planning mode: `auto`, `source`, `chunk`, `none` |
| `--cpu-workers` | `auto` | CPU preprocessing worker count |
| `--target-rows-per-group` | `auto` | Row threshold used by the planner |
| `--target-bytes-per-group` | `auto` | Byte threshold used by the planner |
| `--input-chunk-size` | `10000` | Rows read/emitted per normalized chunk |
| `--max-batch-size` | `16384` | Hard cap for adaptive GPU prediction batch size |
| `--target-gpu-memory-fraction` | `0.8` | Fraction of free VRAM the scheduler tries to fill |
| `--prefetch-queue-size` | `auto` | Ready-chunk queue depth between CPU and GPU stages |
| `--prediction-row-budget` | `8192` | Consecutive ready rows merged into one prediction call |
| `--num-graph-workers` | `auto` | CPU graph workers for the MolE pre-forward graph path |
| `--graph-batch-size` | `1024` | MolE graph / forward mini-batch size |
| `--prefetch-batches` | `2` | Prefetched graph mini-batches per graph worker |
| `--profiling` | disabled | Write stage timing fields into `manifest.json` |

### Reproducibility Knob

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--deterministic-representation` | disabled | Force deterministic CUDA MolE forward for reproducibility checks |

This should keep the same model semantics, but it can change low-level CUDA
execution behavior and reduce throughput.

## `num_graph_workers=0`

`num_graph_workers` controls the MolE CPU graph construction path before MolE
forward:

`SMILES -> RDKit -> graph features -> PyG DataLoader batch`

`--num-graph-workers 0` means:

- do not launch extra DataLoader graph workers
- build graphs synchronously in the main process
- keep the same prediction semantics for a fixed model/backend
- only change throughput and CPU contention characteristics

`--num-graph-workers auto|N>0` means:

- use additional DataLoader workers
- prefetch graph mini-batches ahead of MolE forward
- potentially improve overlap, but also increase CPU overhead and contention

On the current 2080 Ti host, `--num-graph-workers 0` benchmarked faster than
worker auto-selection because extra workers contended with the rest of the
pipeline. Re-benchmark this on larger CPU/GPU hosts.

## Runtime Environment Variables

These environment variables affect screening runs even when you do not pass
equivalent CLI flags:

| Variable | Meaning | Result impact |
| --- | --- | --- |
| `MOLE_CLASSIFIER_BACKEND` | Default backend preference: `auto`, `pickle`, `timber` | Can change result content |
| `MOLE_MOLE_MODEL_PATH` | Override MolE checkpoint directory | Can change result content |
| `MOLE_PICKLE_MODEL_PATH` | Override pickle/XGBoost model path | Can change result content |
| `MOLE_TIMBER_MODEL_DIR` | Override Timber compiled model directory | Can change result content |
| `CUDA_VISIBLE_DEVICES` | Bind the process to a visible GPU subset | Performance/device only |
| `MOLE_TORCH_VERSION` | CUDA wheel version for `pixi run install-cuda-torch` | Install/runtime only |
| `MOLE_TORCH_CUDA_TAG` | CUDA wheel tag for `pixi run install-cuda-torch` | Install/runtime only |
| `MOLE_TORCH_INDEX_URL` | Custom PyTorch wheel index URL | Install/runtime only |

## Practical Recommendations

- Use `--classifier-backend pickle` when you want the original MolE/XGBoost
  path.
- Preprocess large CSV/TSV sources into Parquet shards before screening.
- Use Parquet files or shard directories as the preferred execution format.
- Treat raw `tar.gz` and giant CSV files as transport/storage formats.
- Start with `--num-graph-workers 0` on the current host, then benchmark again
  on larger machines.
- Use `--profiling` when you need to see whether time is dominated by input
  preparation, graph construction, MolE forward, or classifier aggregation.
