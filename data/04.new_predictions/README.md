# New Prediction Batches

Use this directory for molecule batches that need MolE/XGBoost screening.

## Recommended Layout

For a screening batch, keep the raw inputs and generated outputs under a
dedicated run directory:

```text
data/04.new_predictions/<batch_name>/runs/<run_id>/
```

Typical outputs written by `pixi run mole screen`:

- `normalized_input.tsv`
- `predictions_all.tsv`
- `by_source/<source_group>/predictions.tsv`
- `manifest.json`

## Supported Input Types

`mole screen` can ingest:

- `TSV` / `CSV` tables
- `tar.gz` archives containing nested CSV files
- `SQLite` databases with auto-detected tables, or an explicit table/query
  selection when multiple tables match

## Minimum Input Columns

Only `smiles` is required.

If `chem_id` is missing, the tool auto-generates stable IDs such as:

- `mol1`
- `mol2`
- `tylosin_scheme_b_unique_products__1`

All other columns are preserved in the output where possible.

## Usage

```bash
pixi run mole screen \
  --input-path data/04.new_predictions/2026-04-21_screening/macro_split_ring16_scheme_b_fix_pos13_per_scaffold_2026-04-21.tar.gz \
  --output-dir data/04.new_predictions/2026-04-21_screening/runs/demo
```

Prefer `mole screen` for reusable batch screening. Use the legacy batch
predictor scripts only for reproducing historical results.
