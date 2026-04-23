# CLI Reference

This document is the canonical parameter reference for the local `mole` CLI.
It reflects the current implementation in this repository, not the upstream
`rolayoalarcon/mole_antimicrobial_potential` repository.

Use this page when you need:

- the full command surface for `pixi run mole ...`
- exact default behavior for each command
- a quick answer to whether a parameter changes result content or only changes
  throughput, scheduling, or reproducibility

## Scope

This document covers the current canonical CLI commands:

- `mole doctor`
- `mole embed`
- `mole predict`
- `mole preprocess-screening-input`
- `mole benchmark-screening-inputs`
- `mole screen`
- `mole stream-enumeration-screen`
- `mole score`
- `mole optimize`

For production-scale streaming enumeration, also read
[docs/stream_enumeration_screen_runbook.md](/home/lingyuzeng/project/mole_antimicrobial_potential/docs/stream_enumeration_screen_runbook.md).

## Conventions

Effect types used below:

- `content`: can change which rows, molecules, scores, or output fields are produced
- `performance`: intended to change throughput, memory pressure, queueing, device use, or scheduling
- `reproducibility`: intended to keep the same model semantics while changing low-level execution choices
- `output`: only changes where or how results are written
- `validation`: only changes checking behavior or exit status

Important defaults:

- `mole predict` defaults to `per-strain` output
- `mole screen` defaults to `aggregate` output
- `mole screen` supports CSV/TSV, Parquet files, Parquet directories, `tar` / `tar.gz`, and SQLite
- `mole stream-enumeration-screen` writes only resumable hit shards and state files; it does not produce `predictions_all.tsv`
- `--classifier-backend auto` currently prefers `pickle` when the pickle model is available, then falls back to `timber`
- `--num-graph-workers 0` means no extra graph workers; graph construction stays in the main process and should not change prediction semantics for a fixed model/backend

## Runtime Environment Variables

These environment variables affect the CLI even when you do not pass equivalent
flags:

| Variable | Default | Meaning | Effect type |
| --- | --- | --- | --- |
| `MOLE_CLASSIFIER_BACKEND` | `auto` | Default classifier backend preference: `auto`, `pickle`, `timber` | content |
| `MOLE_MOLE_MODEL_PATH` | `pretrained_model/model_ginconcat_btwin_100k_d8000_l0.0001` | Override MolE checkpoint directory | content |
| `MOLE_PICKLE_MODEL_PATH` | `data/03.model_evaluation/MolE-XGBoost-08.03.2024_14.20.pkl` | Override pickle/XGBoost model path | content |
| `MOLE_TIMBER_MODEL_DIR` | `~/.timber/models/mole-antimicrobial` | Override Timber model directory | content |
| `CUDA_VISIBLE_DEVICES` | unset | Bind a process to one GPU or a GPU subset | performance |
| `MOLE_TORCH_VERSION` | unset | Override wheel version for `pixi run install-cuda-torch` | performance |
| `MOLE_TORCH_CUDA_TAG` | unset | Override CUDA wheel tag for `pixi run install-cuda-torch` | performance |
| `MOLE_TORCH_INDEX_URL` | unset | Override PyTorch wheel index URL for `pixi run install-cuda-torch` | performance |

## `mole doctor`

### Purpose

Validate the local runtime, model assets, classifier backends, and the main
REINVENT4 support files.

### Default behavior

- prints a human-readable environment report to stdout
- returns exit code `0` on success
- returns exit code `1` when required assets are missing
- if CUDA is unavailable, prints a warning by default instead of failing

### Parameters

| Parameter | Default | Meaning | Effect type |
| --- | --- | --- | --- |
| `--strict-gpu` | disabled | Fail when CUDA is unavailable instead of only warning | validation |
| `--env-file` | `workflows/reinvent4/configs/local.env.recommended.example` | Optional REINVENT4 env file to validate | validation |
| `--scaffold-file` | `workflows/reinvent4/inputs/scaffolds/mother_scaffold.template.smi` | Optional scaffold file to validate | validation |
| `--objective-file` | `workflows/reinvent4/inputs/objectives/pathogen_group_a.site_reward.prototype.json` | Optional objective file to validate | validation |

### Minimal example

```bash
pixi run mole doctor
```

## `mole embed`

### Purpose

Generate MolE embeddings from one or more SMILES strings.

### Default behavior

- requires one or more `--smiles`
- emits JSON records to stdout by default
- auto-generates `chem_id` values as `mol1`, `mol2`, ...
- runs MolE on `auto` device selection unless overridden

### Parameters

| Parameter | Default | Meaning | Effect type |
| --- | --- | --- | --- |
| `--smiles` | required | One or more SMILES strings | content |
| `--chem-id` | auto-generated | Optional chemical identifiers aligned with `--smiles` | content |
| `--mole-model` | `MOLE_MOLE_MODEL_PATH` or repo default MolE directory | MolE checkpoint directory | content |
| `--device` | `auto` | MolE device: `auto`, `cpu`, `cuda`, or `cuda:<index>` | performance |
| `--deterministic-representation` | disabled | Force deterministic CUDA MolE forward passes | reproducibility |
| `--output` | stdout | Optional output file path | output |
| `--format` | `json` | Output format: `json` or `tsv` | output |

### Minimal example

```bash
pixi run mole embed --smiles CCO
```

## `mole predict`

### Purpose

Predict antimicrobial activity for one or more molecules.

### Default behavior

- requires one or more `--smiles`
- defaults to `per-strain` output
- emits JSON to stdout by default
- auto-generates `chem_id` values as `mol1`, `mol2`, ...
- uses the shared scheduler and predictor singleton

### Parameters

| Parameter | Default | Meaning | Effect type |
| --- | --- | --- | --- |
| `--smiles` | required | One or more SMILES strings | content |
| `--chem-id` | auto-generated | Optional chemical identifiers aligned with `--smiles` | content |
| `--aggregate-scores` | disabled | Return aggregate per-molecule scores instead of per-strain rows | content |
| `--classifier-backend` | env or `auto` | Backend selector: `auto`, `pickle`, `timber` | content |
| `--app-threshold` | `0.04374140128493309` | Growth inhibition threshold | content |
| `--min-nkill` | `10` | Broad-spectrum threshold | content |
| `--num-graph-workers` | `auto` | CPU workers used for `SMILES -> RDKit -> graph -> PyG DataLoader batch` | performance |
| `--graph-batch-size` | `1024` | MolE graph construction and forward mini-batch size | performance |
| `--prefetch-batches` | `2` | Prefetched graph mini-batches per graph worker | performance |
| `--classifier-workers` | `auto` | CPU classifier workers used after MolE embedding generation | performance |
| `--classifier-inflight-batches` | `auto` | Max aggregate-classifier batches kept in flight before ordered merge | performance |
| `--profiling` | disabled | Emit MolE graph-build and prediction runtime profiling | performance |
| `--deterministic-representation` | disabled | Force deterministic CUDA MolE forward passes | reproducibility |
| `--output` | stdout | Optional output file path | output |

### Backend semantics

| Backend | Meaning |
| --- | --- |
| `auto` | Prefer the original pickle/XGBoost classifier when available; otherwise fall back to Timber |
| `pickle` | Use the original pickle/XGBoost classifier path; closest to the original MolE prediction workflow |
| `timber` | Use the compiled Timber backend; may differ numerically from `pickle` |

### `num_graph_workers=0`

This does not disable graph construction. It only disables extra graph worker
processes/threads so graph construction stays in the main process. For a fixed
model and backend, it should not change prediction semantics. It only changes
throughput and CPU resource contention.

### `classifier_workers=auto`

`classifier_workers` controls the CPU classifier/aggregation overlap that runs
after MolE embeddings are produced.

- `auto` is backend-aware
- for `pickle`, it resolves from the effective CPU quota / affinity using a
  conservative heuristic:
  `max(1, min(8, available_cpu_count // 4))`
- for `timber`, it resolves to `1` by default because the compiled backend uses
  shared native state and is treated conservatively until explicit concurrent
  safety is proven
- `classifier_inflight_batches=auto` stays coupled to the resolved worker count
  as `max(2, classifier_workers + 1)`

Use `pixi run python scripts/benchmark_classifier_workers.py ...` to sweep
candidate worker counts on a specific machine before overriding `auto`.

### Minimal example

```bash
pixi run mole predict --smiles CCO
```

## `mole preprocess-screening-input`

### Purpose

Convert large CSV/TSV screening inputs into Parquet shards that are faster and
more predictable for later `mole screen` runs.

### Default behavior

- reads one CSV/TSV input file
- writes Parquet shards under `--output-dir/<source_group>/`
- keeps only `smiles`, `chem_id`, and `source_group`
- auto-fills missing `chem_id` values
- emits a preprocessing manifest as JSON to stdout by default

### Parameters

| Parameter | Default | Meaning | Effect type |
| --- | --- | --- | --- |
| `--input-path` | required | Source CSV/TSV path | content |
| `--output-dir` | required | Root directory where Parquet shards are written | output |
| `--smiles-colname` | `smiles` | SMILES column name in the source table | content |
| `--chem-id-colname` | `chem_id` | chem_id column name in the source table | content |
| `--source-group` | `input` | Logical source label written into each shard row | content |
| `--rows-per-shard` | `100000` | Rows written per Parquet shard | performance |
| `--row-group-size` | `4000` | Parquet row-group size | performance |
| `--output` | stdout | Optional path for writing the preprocessing manifest JSON | output |

### Minimal example

```bash
pixi run mole preprocess-screening-input \
  --input-path data/04.new_predictions/raw/input.csv \
  --output-dir data/04.new_predictions/prepared \
  --source-group batch_a
```

## `mole benchmark-screening-inputs`

### Purpose

Inspect one or more screening input paths before a full run.

### Default behavior

- accepts one or more `--input-path`
- writes JSON to stdout by default
- current implementation is lightweight: it only collects basic path timing and
  file size information
- it does not run MolE prediction and is not a true end-to-end throughput benchmark

### Parameters

| Parameter | Default | Meaning | Effect type |
| --- | --- | --- | --- |
| `--input-path` | required, repeatable | One or more input paths to inspect | validation |
| `--output` | stdout | Optional path for writing benchmark JSON | output |

### Minimal example

```bash
pixi run mole benchmark-screening-inputs \
  --input-path data/04.new_predictions/raw/input.csv \
  --input-path data/04.new_predictions/prepared/batch_a
```

## `mole screen`

### Purpose

Batch-screen large molecule sets for antimicrobial activity.

### Default behavior

- requires `--input-path` and `--output-dir`
- defaults to aggregate output, not per-strain output
- defaults to `--execution-mode thread`
- thread mode writes run artifacts into the output directory:
  - `normalized_input.tsv`
  - `predictions_all.tsv`
  - `by_source/<source_group>/predictions.tsv`
  - `manifest.json`
- process mode first preprocesses repeated CSV/TSV-like `--input-path` values
  into Parquet shards, then writes:
  - `prepared_manifest.json`
  - `prepared_manifests/<source_group>.json`
  - `batch_manifest.jsonl`
  - `run_state.json`
  - `hits/...`
  - `manifest.json`
- process mode persists only broad-spectrum hit rows; it does not write
  `predictions_all.tsv` for non-hit molecules
- supports these input types:
  - thread mode: CSV / TSV files
  - Parquet files
  - Parquet shard directories
  - thread mode only: `tar` / `tar.gz` archives
  - thread mode only: SQLite database files

### Parameters That Can Change Result Content

| Parameter | Default | Meaning | Effect type |
| --- | --- | --- | --- |
| `--input-path` | required | Source CSV/TSV, Parquet file, Parquet directory, `tar` / `tar.gz`, or SQLite file | content |
| `--output-dir` | required | Output directory for run artifacts | output |
| `--execution-mode` | `thread` | `thread` for full normalized/predictions output; `process` for pre-sharded hit-only screening | content |
| `--smiles-colname` | `smiles` | SMILES column name for tabular inputs | content |
| `--chem-id-colname` | `chem_id` | chem_id column name for tabular inputs | content |
| `--archive-pattern` | `*_scheme_b_unique_products.csv` | Archive member filename pattern for tar inputs | content |
| `--archive-smiles-colname` | `product_smiles_canonical` | SMILES column name inside archive CSV members | content |
| `--archive-chem-id-colname` | `example_combo_id` | chem_id column name inside archive CSV members | content |
| `--sqlite-table` | unset | Explicit SQLite table selection | content |
| `--sqlite-query` | unset | Explicit SQLite query selection | content |
| `--classifier-backend` | env or `auto` | Backend selector: `auto`, `pickle`, `timber` | content |
| `--no-dedupe-smiles` | dedupe enabled | Keep duplicate SMILES instead of collapsing to first occurrence | content |
| `--aggregate-scores` | enabled | Return one aggregate row per molecule | content |
| `--per-strain` | disabled | Return per-strain probabilities instead of aggregate rows | content |
| `--app-threshold` | `0.04374140128493309` | Growth inhibition threshold | content |
| `--min-nkill` | `10` | Broad-spectrum threshold | content |

### Parameters That Primarily Affect Throughput or Scheduling

| Parameter | Default | Meaning | Effect type |
| --- | --- | --- | --- |
| `--grouping-mode` | `auto` | Work-unit planning mode: `auto`, `source`, `chunk`, `none` | performance |
| `--cpu-workers` | `auto` | CPU preprocessing worker count | performance |
| `--producer-processes` | `auto` | CPU producer process count in process mode | performance |
| `--predict-queue-max-batches` | `auto` | Bound on queued ready-to-predict batches in process mode | performance |
| `--result-queue-max-batches` | `auto` | Bound on queued predicted batches in process mode | performance |
| `--batch-checkpoint-size` | `2048` | Rows committed per resumable process-mode batch | performance |
| `--rows-per-shard` | `100000` | Rows written per prepared Parquet shard in process mode | performance |
| `--row-group-size` | `4096` | Parquet row-group size used during process-mode preprocessing | performance |
| `--target-rows-per-group` | `auto` | Planner row threshold | performance |
| `--target-bytes-per-group` | `auto` | Planner byte threshold | performance |
| `--input-chunk-size` / `--chunk-size` | `10000` | Rows read and emitted per normalized chunk | performance |
| `--max-batch-size` | `16384` | Hard cap for adaptive GPU prediction batch size | performance |
| `--target-gpu-memory-fraction` | `0.8` | Fraction of free VRAM the scheduler tries to fill | performance |
| `--prefetch-queue-size` | `auto` | Ready-chunk queue depth between CPU and GPU stages | performance |
| `--num-graph-workers` | `auto` | CPU graph workers for the MolE pre-forward graph path | performance |
| `--graph-batch-size` | `1024` | MolE graph / forward mini-batch size | performance |
| `--prefetch-batches` | `2` | Prefetched graph mini-batches per graph worker | performance |
| `--prediction-row-budget` | `8192` | Consecutive ready rows merged into a single prediction call up to this row budget | performance |
| `--profiling` | disabled | Write stage timing into `manifest.json` | performance |
| `--classifier-workers` | `auto` | CPU classifier workers used after MolE embedding generation | performance |
| `--classifier-inflight-batches` | `auto` | Max aggregate-classifier batches kept in flight before ordered merge | performance |

### Reproducibility Parameter

| Parameter | Default | Meaning | Effect type |
| --- | --- | --- | --- |
| `--deterministic-representation` | disabled | Force deterministic CUDA MolE forward passes | reproducibility |

### Backend semantics

| Backend | Meaning |
| --- | --- |
| `auto` | Prefer the original pickle/XGBoost classifier when available; otherwise fall back to Timber |
| `pickle` | Use the original pickle/XGBoost classifier path; closest to the original MolE prediction workflow |
| `timber` | Use the compiled Timber backend; may differ numerically from `pickle` |

### `num_graph_workers=0`

`num_graph_workers` controls the MolE CPU graph construction path before MolE
forward:

`SMILES -> RDKit -> graph features -> PyG DataLoader batch`

`--num-graph-workers 0` means:

- no extra graph workers are launched
- graph construction stays in the main process
- prediction semantics should stay the same for a fixed model/backend
- only throughput and CPU contention change

### `classifier_workers=auto`

For `mole screen`, `classifier_workers=auto` uses the same backend-aware rule
as `mole predict`:

- `pickle`: quota-aware `max(1, min(8, available_cpu_count // 4))`
- `timber`: `1`

This setting only changes throughput. It should not change prediction content
for a fixed model/backend. If you want to tune it explicitly, benchmark a
representative shard first with
`scripts/benchmark_classifier_workers.py --mode predictor` or `--mode stream`.

### Minimal example

```bash
pixi run mole screen \
  --input-path data/04.new_predictions/prepared/batch_a \
  --output-dir data/04.new_predictions/runs/batch_a \
  --classifier-backend pickle \
  --num-graph-workers 0 \
  --graph-batch-size 1024 \
  --prediction-row-budget 8192
```

### Process-mode multi-input example

```bash
pixi run mole screen \
  --input-path data/04.new_predictions/raw/tylosin.csv \
  --input-path data/04.new_predictions/raw/tilmicosin.csv \
  --output-dir data/04.new_predictions/runs/process_batch \
  --execution-mode process \
  --producer-processes auto \
  --predict-queue-max-batches auto \
  --result-queue-max-batches auto \
  --batch-checkpoint-size 2048
```

Process mode is the operational path for multi-source hit-only screening on one
GPU. It assumes aggregate output, writes batch-level resume metadata, and keeps
only broad-spectrum hits on disk. See
[docs/process_screening.md](/home/lingyuzeng/project/mole_antimicrobial_potential/docs/process_screening.md).

## `mole stream-enumeration-screen`

### Purpose

Enumerate scaffold/fragment combinations directly from scaffold and fragment
libraries, predict aggregate antimicrobial scores in bounded batches, and
persist only broad-spectrum hit shards plus resumable run state.

### Default behavior

- requires `--output-dir`, `--ordinary-library`, and `--pos13-library`
- defaults to `workflows/reinvent4/inputs/scaffolds/mother_scaffold.template.smi`
  when no explicit scaffold source is passed
- uses deterministic index mapping:
  `global_combination_index <-> (scaffold_idx, pos3_idx, pos4_idx, pos12_idx, pos13_idx)`
- writes only:
  - `run_state.json`
  - `shard_manifest.jsonl`
  - `hits/<shard_id>.parquet`
- does not write `normalized_input.tsv` or `predictions_all.tsv`
- skips completed shards on rerun and recomputes interrupted/failed shards from
  the shard boundary

### Parameters That Can Change Result Content

| Parameter | Default | Meaning | Effect type |
| --- | --- | --- | --- |
| `--output-dir` | required | Output directory for run state, manifest, and hit shards | output |
| `--scaffold-file` | `workflows/reinvent4/inputs/scaffolds/mother_scaffold.template.smi` | Single scaffold `.smi` file | content |
| `--scaffold-dir` | unset | Directory of scaffold `.smi` files for multi-scaffold runs | content |
| `--scaffold-catalog` | unset | CSV/TSV catalog with `scaffold_slug` and `scaffold_smiles` or `scaffold_file` | content |
| `--ordinary-library` | required | Shared ordinary fragment library CSV for positions `3/4/12` | content |
| `--pos13-library` | required | Position-13 fragment library CSV | content |
| `--run-state` | unset | Optional upstream `run_state.json` copied into provenance | output |
| `--chunk-manifest` | unset | Optional upstream chunk manifest preview copied into provenance | output |
| `--start-index` | `0` | Inclusive global combination start index | content |
| `--stop-index` | full resolved space | Exclusive global combination stop index | content |
| `--classifier-backend` | env or `auto` | Backend selector: `auto`, `pickle`, `timber` | content |
| `--app-threshold` | `0.04374140128493309` | Growth inhibition threshold | content |
| `--min-nkill` | `10` | Broad-spectrum threshold | content |

### Parameters That Primarily Affect Throughput or Scheduling

| Parameter | Default | Meaning | Effect type |
| --- | --- | --- | --- |
| `--shard-size` | `100000` | Combinations per resumable shard | performance |
| `--prediction-batch-size` | `1024` | Combinations per prediction call within a shard | performance |
| `--num-graph-workers` | `auto` | CPU graph workers for the MolE pre-forward graph path | performance |
| `--graph-batch-size` | `1024` | MolE graph / forward mini-batch size | performance |
| `--prefetch-batches` | `2` | Prefetched graph mini-batches per graph worker | performance |
| `--classifier-workers` | `auto` | CPU classifier workers used after MolE embedding generation | performance |
| `--classifier-inflight-batches` | `auto` | Max aggregate-classifier batches kept in flight before ordered merge | performance |
| `--profiling` | disabled | Enable prediction profiling during streaming enumeration | performance |

### Reproducibility Parameter

| Parameter | Default | Meaning | Effect type |
| --- | --- | --- | --- |
| `--deterministic-representation` | disabled | Force deterministic CUDA MolE forward passes | reproducibility |

### Resume semantics

- `run_state.json` stores the current top-level status and the immutable
  parameter snapshot for compatibility checks
- `shard_manifest.jsonl` is the per-shard source of truth
- hit shard writes are atomic: write temp parquet, rename, then update manifest
- only `completed` shards with an existing hit parquet are skipped
- `interrupted` / `failed` shards are recomputed from `start_idx`

### Minimal example

```bash
pixi run mole stream-enumeration-screen \
  --output-dir data/05.stream_tasks/smoke_runs/demo \
  --run-state data/05.stream_tasks/thought2_stream_screen_2026-04-23/thought2_stream_screen_2026-04-23/run_state.json \
  --ordinary-library data/05.stream_tasks/thought2_stream_screen_2026-04-23/thought2_stream_screen_2026-04-23/validation_output/thought2_enumeration/input_libraries/shared_ordinary_library.csv \
  --pos13-library data/05.stream_tasks/thought2_stream_screen_2026-04-23/thought2_stream_screen_2026-04-23/validation_output/thought2_enumeration/input_libraries/pos13_sugar_library.csv \
  --chunk-manifest data/05.stream_tasks/thought2_stream_screen_2026-04-23/thought2_stream_screen_2026-04-23/validation_output/thought2_enumeration/chunk_manifest.csv \
  --stop-index 100000 \
  --shard-size 20000 \
  --prediction-batch-size 512 \
  --classifier-backend pickle \
  --num-graph-workers 0
```

### Auxiliary benchmark script

Use the non-CLI benchmark helper when you need to compare candidate
`classifier_workers` values on a specific machine:

```bash
pixi run python scripts/benchmark_classifier_workers.py \
  --mode predictor \
  --ordinary-library data/05.stream_tasks/thought2_stream_screen_2026-04-23/thought2_stream_screen_2026-04-23/validation_output/thought2_enumeration/input_libraries/shared_ordinary_library.csv \
  --pos13-library data/05.stream_tasks/thought2_stream_screen_2026-04-23/thought2_stream_screen_2026-04-23/validation_output/thought2_enumeration/input_libraries/pos13_sugar_library.csv \
  --workers 1 2 4 6 8 12 \
  --sample-size 8192 \
  --repeats 3 \
  --classifier-backend pickle \
  --num-graph-workers 0 \
  --output data/05.stream_tasks/benchmarks/classifier_workers/predictor_sample8192.json
```

## `mole score`

### Purpose

Compute REINVENT4 reward values from MolE strain-level probabilities.

### Default behavior

- requires `--objective-file` and one or more `--smiles`
- always predicts per-strain probabilities internally, then computes reward
- writes JSON to stdout by default
- auto-generates `chem_id` values as `mol1`, `mol2`, ...
- if `site_reward.enabled=true` and the objective file does not include
  `site_reward.scaffold_smiles`, the CLI fills it from `--scaffold-file`

### Parameters

| Parameter | Default | Meaning | Effect type |
| --- | --- | --- | --- |
| `--objective-file` | required | Path to a normalized objective JSON file | content |
| `--smiles` | required | One or more SMILES strings | content |
| `--chem-id` | auto-generated | Optional chemical identifiers aligned with `--smiles` | content |
| `--scaffold-file` | `workflows/reinvent4/inputs/scaffolds/mother_scaffold.template.smi` | Used only when `site_reward.enabled=true` and `scaffold_smiles` is absent in the objective file | content |
| `--app-threshold` | `0.04374140128493309` | Shared inhibition threshold | content |
| `--min-nkill` | `10` | Broad-spectrum reference threshold reported in metadata | content |
| `--output` | stdout | Optional output file path | output |

### Objective file expectations

The objective JSON is validated against the current local schema:

| Field | Default | Meaning |
| --- | --- | --- |
| `mode` | required | `single_strain` or `broad_spectrum_soft` |
| `strain` | unset | Shortcut for one target strain in `single_strain` mode |
| `strains` | unset | Explicit target strain list |
| `weights` | unset | Optional weights for `single_strain`; not allowed in `broad_spectrum_soft` |
| `tau` | `0.02` | Soft-threshold temperature for `broad_spectrum_soft` |
| `site_reward.enabled` | `false` | Enable auxiliary site reward |
| `site_reward.scaffold_smiles` | unset | Required when `site_reward.enabled=true`, unless the CLI fills it from `--scaffold-file` |
| `site_reward.range_min` | `4` | Preferred lower heavy-atom bound per site |
| `site_reward.range_max` | `27` | Preferred upper heavy-atom bound per site |
| `site_reward.alpha` | `1.5` | Soft lower-bound transition width |
| `site_reward.beta` | `2.5` | Soft upper-bound transition width |
| `site_reward.coverage_weight` | `0.7` | Coverage term weight inside `site_reward` |
| `site_reward.balance_weight` | `0.3` | Balance term weight inside `site_reward` |
| `site_reward.lambda` | `0.85` | Final blend weight in `final_score = lambda * Mole_reward + (1 - lambda) * site_reward` |

### Minimal example

```bash
pixi run mole score \
  --objective-file workflows/reinvent4/inputs/objectives/pathogen_group_a.site_reward.prototype.json \
  --smiles CCO
```

## `mole optimize`

### Purpose

Launch the chunked REINVENT4 optimization workflow wrapper.

### Default behavior

- requires template, env, experiment spec, and scaffold file paths
- shells out to the configured workflow script
- writes child process output directly to the terminal
- returns the child process exit code
- by default uses the chunked long-run runner under `workflows/reinvent4/scripts/run_long_rl.py`

### Parameters

| Parameter | Default | Meaning | Effect type |
| --- | --- | --- | --- |
| `--template` | required | REINVENT4 TOML template path | content |
| `--env-file` | required | Local `.env`-style runtime file | content |
| `--experiment-spec` | required | Long-run RL spec JSON path | content |
| `--scaffold-file` | required | LibInvent scaffold `.smi` file | content |
| `--output-root` | `workflows/reinvent4/results/runs` | Root directory for run outputs | output |
| `--run-id` | empty string | Optional deterministic run identifier | output |
| `--dry-run` | disabled | Render configuration and stop before launching REINVENT4 | validation |
| `--script` | `workflows/reinvent4/scripts/run_long_rl.py` | Workflow script to invoke | content |

### Minimal example

```bash
pixi run mole optimize \
  --template workflows/reinvent4/configs/templates/multi_strain_rl_macrocycle.toml.tpl \
  --env-file workflows/reinvent4/configs/local.env.recommended.example \
  --experiment-spec workflows/reinvent4/configs/long_runs/pathogen_group_a.json \
  --scaffold-file workflows/reinvent4/inputs/scaffolds/mother_scaffold.template.smi
```

## Suggested Reading Order

For later agents, the shortest path is:

1. this file for the full CLI surface
2. [README.md](../README.md) for the quick-start and high-level guidance
3. [docs/batch_screening_input_format.md](batch_screening_input_format.md) for batch-screening input and throughput guidance
4. [docs/repo_layout.md](repo_layout.md) for parameter ownership in code
5. [README_API.md](../README_API.md) if the task moves from CLI to HTTP
