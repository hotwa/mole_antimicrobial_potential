# Process Screening

This document describes the operational model for:

```bash
pixi run mole screen --execution-mode process
```

## Purpose

Use process mode when you want one GPU-bound predictor process plus multiple
CPU producer processes for large multi-source screening runs on a single GPU.

The intended workflow is:

1. Repeat `--input-path` for multiple CSV/TSV-like source files.
2. Pre-shard those inputs into Parquet under the run output directory.
3. Run producer/predictor/writer multiprocessing over the prepared shards.
4. Persist only broad-spectrum hits and resumable batch metadata.

## Current Semantics

- One predictor process owns one GPU and one resident MolE model.
- Producer processes do CPU-heavy input parsing and batch assembly.
- Process mode is aggregate-only. It expects broad-spectrum summaries and uses
  those summaries to decide what to persist.
- Non-hit molecules are counted in manifests and run state, but their detailed
  rows are not written to disk.
- Batch-level resume safety comes from `batch_manifest.jsonl`, not only from
  whole-shard completion.

## Supported Inputs

Process mode currently supports:

- repeated CSV/TSV-like tabular `--input-path` values that can be pre-sharded
- prepared Parquet files
- prepared Parquet directories

Thread mode remains the path for archive and SQLite inputs.

## Output Layout

A process-mode run writes at least:

- `manifest.json`
- `prepared_manifest.json`
- `prepared_manifests/<source_group>.json`
- `prepared/<source_group>/shard_*.parquet`
- `batch_manifest.jsonl`
- `run_state.json`
- `hits/...`

The key difference from thread mode is that process mode does not write
`normalized_input.tsv` or `predictions_all.tsv` for the full population.

## Example

```bash
pixi run mole screen \
  --input-path data/04.new_predictions/raw/tylosin.csv \
  --input-path data/04.new_predictions/raw/tilmicosin.csv \
  --output-dir data/04.new_predictions/runs/process_batch \
  --execution-mode process \
  --classifier-backend pickle \
  --producer-processes auto \
  --predict-queue-max-batches auto \
  --result-queue-max-batches auto \
  --batch-checkpoint-size 2048 \
  --rows-per-shard 100000 \
  --row-group-size 4096
```

## Operational Notes

- Prefer one OS process per GPU via `CUDA_VISIBLE_DEVICES`.
- Do not start multiple predictor processes on the same GPU.
- If you need full row-level outputs for every molecule, use thread mode.
- If you rerun the same output directory, completed batches with existing hit
  parquet outputs are skipped.
