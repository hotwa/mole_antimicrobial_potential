# New Prediction Batches

Use this directory for MolE screening inputs and screening outputs.

This file is the operator-facing guide for real batch runs.

## Recommended Layout

Keep raw inputs, prepared inputs, and screening runs separated:

```text
data/04.new_predictions/<batch_name>/
  raw/
  prepared/
  runs/<run_id>/
```

Typical `mole screen` outputs under `runs/<run_id>/`:

- `normalized_input.tsv`
- `predictions_all.tsv`
- `by_source/<source_group>/predictions.tsv`
- `manifest.json`

`manifest.json` is the first file to inspect when you need to confirm:

- backend selection
- device
- adaptive batch size
- graph worker count
- stage profiling

## Supported Input Types

`mole screen` accepts:

- CSV / TSV files
- Parquet files
- Parquet shard directories
- `tar` / `tar.gz` archives containing CSV members
- SQLite database files

For real throughput-oriented runs, prefer Parquet files or Parquet shard
directories.

## Preferred Workflow

1. Keep the original source under `raw/`.
2. Normalize it to `smiles`, `chem_id`, `source_group`.
3. Convert large CSV/TSV sources into Parquet shards under `prepared/`.
4. Point `mole screen` at the prepared Parquet file or directory.
5. Keep only the run artifacts you actually need under `runs/<run_id>/`.

## Preprocess Large CSV/TSV Inputs

```bash
pixi run mole preprocess-screening-input \
  --input-path data/04.new_predictions/<batch_name>/raw/input.csv \
  --output-dir data/04.new_predictions/<batch_name>/prepared \
  --smiles-colname smiles \
  --chem-id-colname chem_id \
  --source-group <batch_name> \
  --rows-per-shard 100000 \
  --row-group-size 4000
```

Use this preprocessing step when the source is a large CSV/TSV. It keeps only:

- `smiles`
- `chem_id`
- `source_group`

If `chem_id` is missing, stable IDs are generated automatically.

## Recommended Screening Command

```bash
pixi run mole screen \
  --input-path data/04.new_predictions/<batch_name>/prepared/<batch_name> \
  --output-dir data/04.new_predictions/<batch_name>/runs/<run_id> \
  --classifier-backend pickle \
  --aggregate-scores \
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

## Operational Notes

- `--classifier-backend pickle`
  - use this when you want the original MolE/XGBoost prediction path
- `--classifier-backend auto`
  - currently prefers pickle when the pickle model is present
- `--num-graph-workers 0`
  - disables extra graph workers
  - keeps graph construction in the main process
  - should not change prediction semantics for a fixed model/backend
  - only changes throughput and CPU contention
- `--profiling`
  - writes timing fields into `manifest.json`
  - use it when you need to see whether the run is bottlenecked by input
    parsing, graph construction, MolE forward, or aggregation
- `--aggregate-scores`
  - preferred for broad-spectrum screening summaries
- `--per-strain`
  - use only when you need full strain-level outputs
- `--no-dedupe-smiles`
  - keep duplicate SMILES if you need every original row preserved
  - leave dedupe enabled when you want a compact broad-spectrum screening pass

## Parameters That Change Result Content

For run planning, treat these as content-affecting knobs:

- `--classifier-backend`
- `--archive-pattern`
- `--archive-smiles-colname`
- `--archive-chem-id-colname`
- `--sqlite-table`
- `--sqlite-query`
- `--no-dedupe-smiles`
- `--aggregate-scores` / `--per-strain`
- `--app-threshold`
- `--min-nkill`

## Parameters That Primarily Affect Throughput

These are the main operator knobs for performance tuning:

- `--grouping-mode`
- `--cpu-workers`
- `--target-rows-per-group`
- `--target-bytes-per-group`
- `--input-chunk-size`
- `--max-batch-size`
- `--target-gpu-memory-fraction`
- `--prefetch-queue-size`
- `--prediction-row-budget`
- `--num-graph-workers`
- `--graph-batch-size`
- `--prefetch-batches`
- `--profiling`
- `--deterministic-representation`

## Runtime Environment Variables

These environment variables are relevant for screening jobs:

- `MOLE_CLASSIFIER_BACKEND`
- `MOLE_MOLE_MODEL_PATH`
- `MOLE_PICKLE_MODEL_PATH`
- `MOLE_TIMBER_MODEL_DIR`
- `CUDA_VISIBLE_DEVICES`
- `MOLE_TORCH_VERSION`
- `MOLE_TORCH_CUDA_TAG`
- `MOLE_TORCH_INDEX_URL`

Use `CUDA_VISIBLE_DEVICES` to pin one screening process to one GPU. Use
`MOLE_TORCH_*` only when installing the CUDA PyTorch wheel on a new machine.

## More Detail

For the full screening parameter handbook, read:

- [docs/batch_screening_input_format.md](/home/lingyuzeng/project/mole_antimicrobial_potential/docs/batch_screening_input_format.md)
- [docs/cli_reference.md](/home/lingyuzeng/project/mole_antimicrobial_potential/docs/cli_reference.md)

For file ownership and parameter landing points in code, read:

- [docs/repo_layout.md](/home/lingyuzeng/project/mole_antimicrobial_potential/docs/repo_layout.md)
