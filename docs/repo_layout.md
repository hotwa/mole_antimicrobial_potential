# Repository Layout

This repository now keeps only a few canonical entrypoints at the root. Most
implementation logic lives under `src/`, and older compatibility scripts live
under `scripts/legacy/`.

## Project Overview

This project is a MolE-based antimicrobial discovery system. It takes SMILES
strings as input, converts them into MolE embeddings, and then uses the original
pickle/XGBoost classifier or an optional Timber backend to estimate antimicrobial activity across many
bacterial strains. The same core predictor powers the CLI, the FastAPI server,
the MCP server, and the REINVENT4 scoring workflow.

The main capabilities are:

- MolE embedding generation for small molecules
- Per-strain antimicrobial probability prediction
- Aggregated antimicrobial scoring, including broad-spectrum summaries
- REINVENT4 reward generation for single-strain, multi-strain, and
  broad-spectrum optimization
- Chunked long-run RL workflows with scaffold-aware optimization

If you are a new agent, this is the page to read first when you need to
understand how the repository is organized and which file owns which feature.

## Usage Summary

The shortest path for common tasks is:

- Bootstrap a CUDA-capable runtime on a new NVIDIA machine:
  - `pixi install`
  - `pixi run install-cuda-torch`
  - optionally override the CUDA wheel with `MOLE_TORCH_CUDA_TAG=cu124`
    and `MOLE_TORCH_VERSION=2.5.1+cu124`
- Environment and model sanity check:
  - `pixi run mole doctor`
- MolE embedding:
  - `pixi run mole embed --smiles CCO`
- Antimicrobial prediction:
  - `pixi run mole predict --smiles CCO`
- Batch screening from CSV/TSV/archive inputs:
  - `pixi run mole screen --input-path data/04.new_predictions/2026-04-21_screening/macro_split_ring16_scheme_b_fix_pos13_per_scaffold_2026-04-21.tar.gz --output-dir data/04.new_predictions/2026-04-21_screening/runs/demo`
- Preferred high-throughput preprocessing:
  - `pixi run mole preprocess-screening-input --input-path <raw_csv> --output-dir <prepared_dir> --smiles-colname <smiles_col> --chem-id-colname <chem_id_col> --source-group <group_name>`
- REINVENT4 reward computation:
  - `pixi run mole score --objective-file workflows/reinvent4/inputs/objectives/pathogen_group_a.site_reward.prototype.json --smiles CCO`
- REINVENT4 optimization:
  - `pixi run mole optimize --template workflows/reinvent4/configs/templates/multi_strain_rl_macrocycle.toml.tpl --env-file workflows/reinvent4/configs/local.env.recommended.example --experiment-spec workflows/reinvent4/configs/long_runs/pathogen_group_a.json --scaffold-file workflows/reinvent4/inputs/scaffolds/mother_scaffold.template.smi`
- FastAPI service:
  - `python api_server.py`
- MCP service:
  - `python mcp_server_enhanced.py`

Use the canonical CLI/API/MCP entrypoints above for new work. Use the legacy
wrappers only when reproducing old behavior or comparing against historical
scripts.

For batch-screening input preparation guidance, read:

- [docs/cli_reference.md](/home/lingyuzeng/project/mole_antimicrobial_potential/docs/cli_reference.md)
- [docs/batch_screening_input_format.md](/home/lingyuzeng/project/mole_antimicrobial_potential/docs/batch_screening_input_format.md)

That document is the canonical explanation of why raw `tar.gz` bundles and one
huge CSV are transport/storage formats, not preferred screening formats.

## Root Directory Policy

The root should stay small and predictable:

- keep only canonical entrypoints and repository-level configuration at the root
- keep repository-level agent instruction files such as `AGENTS.md` and `CLAUDE.md` at the root
- keep implementation under `src/`
- keep workflow-specific assets under `workflow/` or `workflows/`
- keep one-off recovery artifacts under `.trash/`
- keep auxiliary tool subprojects under `tools/`

If a file is not a canonical entrypoint, project-level config, or a standard
repository document, it should usually not stay at the root.

## Canonical Entrypoints

Use these first for new work:

- [`mole_cli.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/mole_cli.py)
  - Unified CLI for MolE, Timber, and REINVENT4
  - Recommended usage: `pixi run mole ...`
- [`api_server.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/api_server.py)
  - FastAPI app exposing `/health`, `/predict`, and `/score`
  - Recommended usage: `python api_server.py`
- [`mcp_server_enhanced.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/mcp_server_enhanced.py)
  - Streamable HTTP FastMCP server
  - Recommended usage: `python mcp_server_enhanced.py`

Additional layout notes:

- [`docker/nginx/`](/home/lingyuzeng/project/mole_antimicrobial_potential/docker/nginx)
  - Nginx configuration used by the legacy container deployment path
- [`tools/timber/`](/home/lingyuzeng/project/mole_antimicrobial_potential/tools/timber)
  - Timber conversion/export helper subproject
- [`docs/references/papers/`](/home/lingyuzeng/project/mole_antimicrobial_potential/docs/references/papers)
  - Local reference PDFs and related literature assets
- [`docs/superpowers/plans/`](/home/lingyuzeng/project/mole_antimicrobial_potential/docs/superpowers/plans)
  - Retained implementation plans for substantial refactors and throughput work
  - Not canonical user-facing docs; use for historical execution context only
- [`.trash/`](/home/lingyuzeng/project/mole_antimicrobial_potential/.trash)
  - Quarantine area for ad hoc scripts, scratch inputs, and runtime leftovers
  - `.trash/test_ad_hoc/` is the preferred holding area for one-off test probes
    that are not part of the maintained `test/` suite

## Parameter Ownership Index

If you need to understand where a CLI knob actually lands in the code, use this
index first.

- [`mole_cli.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/mole_cli.py)
  - owns the public `mole` CLI surface
  - defines `mole predict`, `mole screen`, `mole preprocess-screening-input`,
    `mole benchmark-screening-inputs`, `mole doctor`, `mole score`,
    `mole optimize`
  - important parameter families:
    `--classifier-backend`, `--num-graph-workers`, `--graph-batch-size`,
    `--prefetch-batches`, `--profiling`,
    `--deterministic-representation`, `--grouping-mode`, `--cpu-workers`,
    `--input-chunk-size`, `--max-batch-size`,
    `--target-gpu-memory-fraction`, `--prediction-row-budget`
- [`src/prediction_scheduler.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/prediction_scheduler.py)
  - owns adaptive GPU micro-batching and runtime snapshots
  - consumes GPU-facing knobs such as `--max-batch-size`,
    `--target-gpu-memory-fraction`, device selection, and profiling snapshots
- [`src/predictor.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/predictor.py)
  - bridges CLI/API requests to MolE representation and classifier backends
  - forwards `--classifier-backend`, `--num-graph-workers`,
    `--graph-batch-size`, `--prefetch-batches`,
    `--deterministic-representation`, and profiling enablement
- [`workflow/dataset/dataset_representation.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/workflow/dataset/dataset_representation.py)
  - owns the MolE pre-forward CPU graph path:
    `SMILES -> RDKit -> graph -> PyG DataLoader`
  - this is where `num_graph_workers`, `graph_batch_size`,
    `prefetch_batches`, `pin_memory`, and deterministic CUDA controls become
    DataLoader / execution settings
  - `num_graph_workers=0` means no extra graph workers; graph construction stays
    in the main process and should not change prediction semantics
- [`src/screening_planner.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/screening_planner.py)
  - owns input work-unit planning for `mole screen`
  - consumes `--grouping-mode`, `--target-rows-per-group`,
    `--target-bytes-per-group`, `--cpu-workers`, and source-specific selectors
- [`src/screening_sources.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/screening_sources.py)
  - owns archive / SQLite / tabular / parquet normalization into screening
    chunks
  - consumes `--archive-pattern`, `--archive-smiles-colname`,
    `--archive-chem-id-colname`, `--sqlite-table`, `--sqlite-query`,
    `--smiles-colname`, `--chem-id-colname`, `--input-chunk-size`
- [`src/batch_screening.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/batch_screening.py)
  - owns producer/consumer flow for `mole screen`
  - consumes `--no-dedupe-smiles`, `--aggregate-scores` / `--per-strain`,
    `--prefetch-queue-size`, `--prediction-row-budget`, output writing, and
    manifest runtime fields
- [`src/preprocess_screening_input.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/preprocess_screening_input.py)
  - owns `mole preprocess-screening-input`
  - consumes `--input-path`, `--output-dir`, `--smiles-colname`,
    `--chem-id-colname`, `--source-group`, `--rows-per-shard`,
    `--row-group-size`
- [`src/classifier_backend.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/classifier_backend.py)
  - owns backend discovery and backend-specific model paths
  - reads `MOLE_CLASSIFIER_BACKEND`, `MOLE_PICKLE_MODEL_PATH`,
    `MOLE_TIMBER_MODEL_DIR`

## Runtime Environment Index

These environment variables matter most when moving between machines:

- `MOLE_CLASSIFIER_BACKEND`
  - default classifier backend preference
  - can change result content
- `MOLE_MOLE_MODEL_PATH`
  - override MolE checkpoint directory
  - can change result content
- `MOLE_PICKLE_MODEL_PATH`
  - override pickle/XGBoost model path
  - can change result content
- `MOLE_TIMBER_MODEL_DIR`
  - override Timber compiled model directory
  - can change result content
- `CUDA_VISIBLE_DEVICES`
  - bind a process to a specific GPU or GPU subset
  - affects device placement and throughput, not intended prediction semantics
- `MOLE_TORCH_VERSION`, `MOLE_TORCH_CUDA_TAG`, `MOLE_TORCH_INDEX_URL`
  - control `pixi run install-cuda-torch`
  - install/runtime only

## Core Library Modules

Import these from other code instead of duplicating behavior:

- [`src/mole_representation.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/mole_representation.py)
  - MolE SMILES-to-embedding pipeline
- [`src/predictor.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/predictor.py)
  - Shared antimicrobial predictor service
- [`src/prediction_scheduler.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/prediction_scheduler.py)
  - Adaptive single-GPU batch scheduler
- [`src/screening_sources.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/screening_sources.py)
  - Streaming archive, CSV, and SQLite multi-worker data loaders
- [`src/screening_planner.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/screening_planner.py)
  - Planners that dynamically split large sources into parallel manageable work-units
- [`src/batch_screening.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/batch_screening.py)
  - Batch screening loop managing ThreadPoolExecutor and Single GPU producer/consumer stream
- [`src/classifier_backend.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/classifier_backend.py)
  - Pickle/XGBoost and Timber backend selection; `auto` prefers pickle when available
- [`src/models.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/models.py)
  - Request schemas, objective normalization, and validation
- [`src/reinvent_scoring.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/reinvent_scoring.py)
  - REINVENT4 reward computation
- [`src/site_reward.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/site_reward.py)
  - Optional site-aware auxiliary reward prototype
- [`src/reinvent4_workflow.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/reinvent4_workflow.py)
  - REINVENT4 config rendering, chunk orchestration, and validation
- [`src/service.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/service.py)
  - Predictor singleton shared by API and MCP servers

## Batch Data Layout

Batch screening directories under
[`data/04.new_predictions/`](/home/lingyuzeng/project/mole_antimicrobial_potential/data/04.new_predictions)
should keep staging assets, benchmark artifacts, and final run outputs
separated:

- `raw/`
  - original source files kept for provenance
- `prepared/`
  - canonical screening-ready inputs such as Parquet shards
- `cache/`
  - transient or staging assets that support runs but are not final outputs
  - examples: extracted archive trees, copied `inputs/`, copied `prepared/`
- `benchmarks/`
  - throughput probes, patch comparison runs, and their logs
  - these should not be mixed with canonical long-running screening outputs
- `runs/<run_id>/`
  - canonical screening outputs intended for downstream consumption
  - examples: `manifest.json`, `normalized_input.tsv`,
    `predictions_all.tsv`, `by_source/...`

When reorganizing an old batch directory, prefer moving temporary extraction,
probe runs, and copied staging assets into `cache/` or `benchmarks/` instead of
leaving them beside canonical `runs/`.

## Legacy Compatibility Scripts

These scripts are kept for reproduction or reference. Prefer the canonical
entrypoints above for new work.

- [`scripts/legacy/mole_representation.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/scripts/legacy/mole_representation.py)
  - Legacy MolE embedding CLI
- [`scripts/legacy/mole_antimicrobial_prediction.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/scripts/legacy/mole_antimicrobial_prediction.py)
  - Legacy TSV batch predictor
- [`scripts/legacy/predict_api.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/scripts/legacy/predict_api.py)
  - Legacy FastAPI app
- [`scripts/legacy/mcp_server.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/scripts/legacy/mcp_server.py)
  - Legacy FastMCP app
- [`scripts/legacy/health_check.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/scripts/legacy/health_check.py)
  - Legacy health probe helper

## REINVENT4 Workflow Scripts

The optimization workflow lives under:

- [`workflows/reinvent4/scripts/`](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts)

Main utilities:

- [`run_reinvent4.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/run_reinvent4.py)
  - Render config and launch a single REINVENT4 run
- [`run_long_rl.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/run_long_rl.py)
  - Chunked long-run RL orchestration
- [`score_bridge.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/score_bridge.py)
  - Bridge REINVENT4 batches to local `/score`
- [`evaluate_results.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/evaluate_results.py)
  - Re-score and summarize generated molecules
- [`analyze_rgroup_sites.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/analyze_rgroup_sites.py)
  - Inspect per-site decoration sizes
- [`validate_scaffold.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/validate_scaffold.py)
  - Validate scaffold attachment points

## Legacy Manual Test Scripts

These are kept for historical/manual checks and are not part of the default
automated test suite.

- [`test/legacy/test_improvements.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/test/legacy/test_improvements.py)
  - Old integration script for `/health`, `/mcp`, and validation flows
- [`test/legacy/test_mcp_connection.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/test/legacy/test_mcp_connection.py)
  - Old SSE/MCP connection smoke test

## Quick Rules

- Use `pixi run mole doctor` before running anything on a new machine.
- If CUDA is missing after `pixi install`, run `pixi run install-cuda-torch`
  before the first `mole embed` or `mole predict` call.
- Use `pixi run mole predict` and `pixi run mole score` for day-to-day work.
- Use `pixi run mole screen` for reusable batch screening of broad-spectrum candidates.
- Use `python api_server.py` only when you specifically need the HTTP server.
- Use `python mcp_server_enhanced.py` only when you specifically need MCP.
- Use the scripts in `scripts/legacy/` only when reproducing older flows.
- Use `workflows/reinvent4/scripts/` for all REINVENT4-specific automation.
