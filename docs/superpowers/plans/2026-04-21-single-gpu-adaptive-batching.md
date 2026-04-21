# Single-GPU Adaptive Batching Prediction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `mole screen`, `mole predict`, the FastAPI `/predict` endpoint, and the MCP prediction tool use one single-GPU adaptive batching path that auto-tunes batch size from available VRAM and model footprint, so a single GPU is kept busy without requiring multiple Python processes on the same device.

**Architecture:** Split the current prediction path into three layers: input normalization, GPU batch scheduling, and prediction execution. The new scheduler runs one process per GPU, keeps exactly one MolE model resident per process, uses CPU-side prefetching to keep the GPU fed, autotunes batch size from `torch.cuda.mem_get_info()` and a warmup pass, and backs off on OOM. Multi-GPU remains manual: users launch one process per GPU and set `CUDA_VISIBLE_DEVICES` themselves.

**Tech Stack:** Python 3.10, Pixi, PyTorch CUDA wheels, pandas, asyncio/threading, sqlite3, tarfile, FastAPI, FastMCP, unittest.

---

### Task 1: Add a dedicated GPU batch scheduler

**Files:**
- Create: `src/prediction_scheduler.py`
- Modify: `src/service.py`
- Test: `test/test_prediction_scheduler.py`

- [ ] **Step 1: Write the failing tests**
  - Add a test that mocks `torch.cuda.mem_get_info()` for a 22 GB GPU and a 32 GB GPU and asserts that the scheduler picks different batch sizes for the two memory budgets.
  - Add a test that mocks the predictor to raise `torch.cuda.OutOfMemoryError` on an oversized batch and asserts that the scheduler halves the batch size and retries successfully.
  - Add a test that verifies the MolE checkpoint is loaded once per process and reused across multiple scheduled batches.
  - Add a test that verifies the scheduler records runtime metadata, including `gpu_name`, `device`, `torch_version`, `total_memory_bytes`, `free_memory_bytes`, `selected_batch_size`, and `used_cuda`.

- [ ] **Step 2: Run the tests and confirm they fail for the right reason**
  - Run: `pixi run test-prediction-scheduler`
  - Expected: fail because `PredictionScheduler` does not exist yet or because the batch-sizing and retry behavior is still missing.

- [ ] **Step 3: Implement the scheduler minimally**
  - Create `PredictionScheduler` with a single GPU worker per process.
  - Implement `select_batch_size(...)` from available GPU memory, a target memory fraction, and a hard max batch cap.
  - Implement a warmup path that estimates per-item memory cost once and caches the result for the process.
  - Implement OOM backoff by retrying with a smaller batch until the batch fits or the minimum batch size is reached.
  - Keep the MolE model resident and reuse it for all scheduled batches in the process.
  - Expose a runtime snapshot object from the scheduler so callers can include the selected batch size and device details in manifests and outputs.

- [ ] **Step 4: Re-run the tests until they pass**
  - Run: `pixi run test-prediction-scheduler`
  - Expected: PASS, with scheduler batch selection, OOM backoff, and model reuse all covered.

- [ ] **Step 5: Commit the scheduler layer**
  - Commit message: `feat: add single-gpu prediction scheduler`

### Task 2: Split ingestion from scoring and stream batches into the scheduler

**Files:**
- Create: `src/screening_sources.py`
- Modify: `src/batch_screening.py`
- Test: `test/test_screening_sources.py`, `test/test_batch_screening.py`

- [ ] **Step 1: Write the failing tests**
  - Add a test for CSV/TSV input that provides only `smiles` and asserts that `chem_id` is auto-generated as `mol1`, `mol2`, and so on.
  - Add a test for tar.gz input that asserts archive members under `*_scheme_b_unique_products.csv` are loaded with `source_group`, `source_file`, `source_member`, `source_row`, and auto-generated `chem_id` when missing.
  - Add a test for SQLite input that asserts the loader auto-discovers a unique table with a SMILES column, and that ambiguous table detection fails with a clear error unless `--sqlite-table` or `--sqlite-query` is provided.
  - Add a test that `screen_path()` no longer waits to materialize a single giant frame before the first prediction batch is sent to the scheduler.
  - Add a test that `manifest.json` records the normalized row count, selected batch size, GPU metadata, and per-source grouped outputs.

- [ ] **Step 2: Run the tests and confirm they fail for the right reason**
  - Run: `pixi run test-screening-sources`
  - Run: `pixi run test-batch-screening`
  - Expected: fail because the source adapters and streaming batch flow are not yet implemented.

- [ ] **Step 3: Implement the source adapters and streaming flow**
  - Move the input normalization logic into `src/screening_sources.py` so CSV/TSV, tar.gz, and SQLite all normalize to the same record shape.
  - Keep the rule that `smiles` is required and `chem_id` is optional; auto-generate IDs when missing.
  - Make `src/batch_screening.py` consume normalized chunks instead of forcing a full in-memory concat before prediction.
  - Feed each normalized chunk into the new scheduler so the GPU can start working while the archive is still being read.
  - Preserve the current output contract:
    - `normalized_input.tsv`
    - `predictions_all.tsv`
    - `by_source/<source_group>/predictions.tsv`
    - `manifest.json`

- [ ] **Step 4: Re-run the tests until they pass**
  - Run: `pixi run test-screening-sources`
  - Run: `pixi run test-batch-screening`
  - Expected: PASS, with auto-generated `chem_id` and streaming archive handling preserved.

- [ ] **Step 5: Commit the ingestion/screening refactor**
  - Commit message: `feat: stream batch screening inputs into the gpu scheduler`

### Task 3: Wire the scheduler into CLI, API, and MCP entrypoints

**Files:**
- Modify: `mole_cli.py`
- Modify: `api_server.py`
- Modify: `mcp_server_enhanced.py`
- Modify: `src/service.py`
- Create: `test/test_cli_scheduler.py`
- Create: `test/test_service_scheduler.py`

- [ ] **Step 1: Write the failing tests**
  - Add a CLI test that calls `mole screen` on a small input and asserts the command accepts the new tuning flags:
    - `--input-chunk-size`
    - `--max-batch-size`
    - `--target-gpu-memory-fraction`
    - `--prefetch-queue-size`
  - Add a test that `--chunk-size` still works as a legacy alias for the input chunk size so older scripts do not break.
  - Add a test that `mole predict` includes `prediction_runtime` fields showing the selected batch size and whether CUDA was used.
  - Add a test that API and MCP reuse the same predictor/scheduler singleton instead of creating separate MolE loads per request path.

- [ ] **Step 2: Run the tests and confirm they fail for the right reason**
  - Run: `pixi run test-cli`
  - Run: `pixi run test-score`
  - Run: `pixi run test-batch-screening`
  - Expected: fail because the new tuning flags and scheduler wiring are not yet exposed everywhere.

- [ ] **Step 3: Wire the entrypoints**
  - Update `mole_cli.py` so `screen` exposes the new tuning knobs and passes them to the scheduler.
  - Keep `predict` and `/predict` behavior stable, but route them through the same singleton predictor path so the MolE model is loaded once per process.
  - Ensure `api_server.py` and `mcp_server_enhanced.py` use `src/service.py` as the canonical place to obtain the shared predictor.
  - Ensure multi-GPU scaling remains manual: the docs should explicitly say one process per GPU and `CUDA_VISIBLE_DEVICES` is the supported way to spread across GPUs.
  - Add `prediction_runtime` to the returned payload or manifest so users can see the selected batch size, device, and CUDA status.

- [ ] **Step 4: Re-run the tests until they pass**
  - Run: `pixi run test-cli`
  - Run: `pixi run test-score`
  - Run: `pixi run test-batch-screening`
  - Run: `pixi run test-classifier-backend`
  - Expected: PASS, with the scheduler path shared across CLI/API/MCP and no duplicate MolE loads per request path.

- [ ] **Step 5: Commit the wiring**
  - Commit message: `feat: route prediction entrypoints through the adaptive gpu scheduler`

### Task 4: Update docs and runtime guidance for single-GPU and multi-GPU usage

**Files:**
- Modify: `README.md`
- Modify: `README_API.md`
- Modify: `docs/repo_layout.md`
- Modify: `data/04.new_predictions/README.md`
- Modify: `pixi.toml`

- [ ] **Step 1: Write the documentation checks**
  - Verify the docs explicitly say:
    - `pixi install`
    - `pixi run install-cuda-torch`
    - `mole doctor` should report CUDA availability
    - `mole screen` uses one process per GPU and autotunes batch size on that GPU
    - multi-GPU means “start multiple processes and pin each one with `CUDA_VISIBLE_DEVICES`”
  - Verify the docs explain the new tuning knobs and what they do:
    - `--input-chunk-size`
    - `--max-batch-size`
    - `--target-gpu-memory-fraction`
    - `--prefetch-queue-size`

- [ ] **Step 2: Update the docs**
  - Add a short performance section to `README.md` that explains why one process on one GPU is the default and why same-GPU multi-process is not the recommended first choice.
  - Update `README_API.md` so the API entrypoints are described as sharing the same batch scheduler as the CLI.
  - Update `docs/repo_layout.md` so new agents can find the scheduler and streaming input modules quickly.
  - Update `data/04.new_predictions/README.md` with the expected directory structure for long-running screening jobs and the runtime manifest fields.
  - Add `pixi` task names for the new tests and a smoke command for `mole screen` on a small sample input.

- [ ] **Step 3: Re-run the relevant tests and smoke checks**
  - Run: `pixi run test-cli`
  - Run: `pixi run test-batch-screening`
  - Run: `pixi run mole doctor`
  - Run: `pixi run mole screen --input-path <small.tsv> --output-dir <tmpdir>`
  - Expected: documentation matches the actual CLI shape and the smoke run writes the standard output files.

- [ ] **Step 4: Commit the docs and task metadata**
  - Commit message: `docs: document adaptive single-gpu screening and tuning knobs`

## Test Plan

- Unit tests must prove batch selection scales with GPU memory on both a 22 GB class card and a 32 GB class card.
- Unit tests must prove OOM backoff reduces batch size and retries instead of crashing the whole run.
- Unit tests must prove missing `chem_id` is auto-generated and archive / SQLite sources keep their provenance metadata.
- Integration tests must prove `mole screen` writes the total table first, then grouped outputs, and records the selected batch size in the manifest.
- CLI tests must prove the new tuning flags work and the legacy `--chunk-size` alias does not break existing scripts.
- API / MCP tests must prove both code paths reuse the same singleton predictor and do not duplicate model loads in the same process.

## Assumptions and Defaults

- v1 is **single GPU, single process per GPU** only.
- Same-GPU multi-process is intentionally out of scope for v1 because it duplicates model weights and tends to reduce throughput.
- Multi-GPU scaling is manual and documented: one process per GPU, each pinned by `CUDA_VISIBLE_DEVICES`.
- `smiles` remains the only required input field for screening; `chem_id` is optional and auto-generated when missing.
- `mole screen` is the primary throughput target, but `mole predict`, `/predict`, and MCP should all reuse the same predictor singleton and scheduler plumbing.
- The tuning knobs should remain opt-in overrides; the default path should autotune batch size automatically from the GPU memory budget and model footprint.
