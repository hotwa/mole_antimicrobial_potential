# Multiprocess Hit-Only Screening Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a single-GPU, multi-CPU-process screening pipeline that can read multiple per-scaffold CSV inputs, pre-shard them into efficient Parquet sources, feed one MolE predictor process per GPU, and persist only broad-spectrum hit molecules with resumable, idempotent checkpoints.

**Architecture:** Keep the existing invariant that one predictor process owns one GPU and one resident MolE model. Move CPU-heavy CSV parsing and batch assembly into a bounded multiprocessing producer pool, with Parquet shards as the stable intermediate format. Add batch-level manifests so restart safety does not depend on finishing a whole shard before progress is visible.

**Tech Stack:** Python 3.10, Pixi, pandas, pyarrow, multiprocessing `spawn`, asyncio, PyTorch CUDA, RDKit, MolE, XGBoost, unittest.

---

### Task 1: Add multi-input CLI and process-screening configuration

**Files:**
- Create: `src/screening_process_pipeline.py`
- Modify: `mole_cli.py`
- Modify: `src/batch_screening.py`
- Test: `test/test_cli.py`

- [ ] **Step 1: Write the failing CLI/config tests**

```python
def test_screen_parser_accepts_repeated_input_paths():
    from mole_cli import _build_parser

    parser = _build_parser()
    args = parser.parse_args(
        [
            "screen",
            "--input-path", "a.csv",
            "--input-path", "b.csv",
            "--output-dir", "out",
            "--execution-mode", "process",
        ]
    )

    assert args.input_path == ["a.csv", "b.csv"]
    assert args.execution_mode == "process"
```

```python
def test_screen_parser_accepts_process_tuning_flags():
    from mole_cli import _build_parser

    parser = _build_parser()
    args = parser.parse_args(
        [
            "screen",
            "--input-path", "a.csv",
            "--output-dir", "out",
            "--execution-mode", "process",
            "--producer-processes", "auto",
            "--predict-queue-max-batches", "auto",
            "--result-queue-max-batches", "8",
            "--batch-checkpoint-size", "2048",
        ]
    )

    assert args.producer_processes == "auto"
    assert args.predict_queue_max_batches == "auto"
    assert args.result_queue_max_batches == "8"
    assert args.batch_checkpoint_size == 2048
```

- [ ] **Step 2: Run tests to verify failure**

Run:

```bash
pixi run python -m unittest discover -s test -p 'test_cli.py' -v
```

Expected: FAIL because `screen` does not yet accept repeated `--input-path` or process-mode tuning flags.

- [ ] **Step 3: Add the configuration surface**

```python
@dataclass(frozen=True)
class ProcessScreenConfig:
    input_paths: list[str]
    output_dir: str
    execution_mode: str = "thread"
    producer_processes: int | str = "auto"
    predict_queue_max_batches: int | str = "auto"
    result_queue_max_batches: int | str = "auto"
    batch_checkpoint_size: int = 2048
    rows_per_shard: int = 100_000
    row_group_size: int = 4096
```

```python
screen.add_argument(
    "--input-path",
    action="append",
    required=True,
    help="One or more CSV/TSV/Parquet/SQLite/tar inputs. Repeat the flag for multiple sources.",
)
screen.add_argument(
    "--execution-mode",
    choices=["thread", "process"],
    default="thread",
    help="Use thread mode for current behavior or process mode for CPU-heavy producer pools.",
)
screen.add_argument("--producer-processes", default="auto")
screen.add_argument("--predict-queue-max-batches", default="auto")
screen.add_argument("--result-queue-max-batches", default="auto")
screen.add_argument("--batch-checkpoint-size", type=int, default=2048)
```

- [ ] **Step 4: Re-run the CLI tests**

Run:

```bash
pixi run python -m unittest discover -s test -p 'test_cli.py' -v
```

Expected: PASS for the new parser/config surface.

- [ ] **Step 5: Commit**

```bash
git add mole_cli.py src/batch_screening.py src/screening_process_pipeline.py test/test_cli.py
git commit -m "feat: add multiprocess screening cli configuration"
```

### Task 2: Pre-shard multiple CSV sources into Parquet files with explicit source groups

**Files:**
- Modify: `src/preprocess_screening_input.py`
- Modify: `src/screening_sources.py`
- Create: `test/test_preprocess_screening_input.py`

- [ ] **Step 1: Write the failing preprocessing tests**

```python
def test_preprocess_multiple_csvs_writes_one_source_group_per_input(tmp_path):
    from src.preprocess_screening_input import preprocess_many_to_parquet

    a = tmp_path / "tylosin.csv"
    b = tmp_path / "tilmicosin.csv"
    a.write_text("product_smiles_canonical,example_combo_id\nCCO,t1\n", encoding="utf-8")
    b.write_text("product_smiles_canonical,example_combo_id\nCCN,m1\n", encoding="utf-8")

    manifest = preprocess_many_to_parquet(
        input_paths=[a, b],
        output_dir=tmp_path / "prepared",
        smiles_colname="product_smiles_canonical",
        chem_id_colname="example_combo_id",
        rows_per_shard=1,
        row_group_size=1,
    )

    assert manifest["source_groups"] == ["tilmicosin", "tylosin"]
    assert (tmp_path / "prepared" / "tylosin" / "shard_0001.parquet").exists()
    assert (tmp_path / "prepared" / "tilmicosin" / "shard_0001.parquet").exists()
```

```python
def test_parquet_shard_keeps_only_required_columns(tmp_path):
    from src.preprocess_screening_input import preprocess_many_to_parquet
    import pandas as pd

    path = tmp_path / "tylosin.csv"
    path.write_text(
        "product_smiles_canonical,example_combo_id,unused\nCCO,t1,x\n",
        encoding="utf-8",
    )

    preprocess_many_to_parquet(
        input_paths=[path],
        output_dir=tmp_path / "prepared",
        smiles_colname="product_smiles_canonical",
        chem_id_colname="example_combo_id",
        rows_per_shard=10,
        row_group_size=10,
    )

    frame = pd.read_parquet(tmp_path / "prepared" / "tylosin" / "shard_0001.parquet")
    assert list(frame.columns) == ["smiles", "chem_id", "source_group", "source_file", "source_row"]
```

- [ ] **Step 2: Run tests to verify failure**

Run:

```bash
pixi run python -m unittest discover -s test -p 'test_preprocess_screening_input.py' -v
```

Expected: FAIL because multi-input preprocessing and source-group Parquet shards do not exist yet.

- [ ] **Step 3: Implement multi-input preprocessing**

```python
def preprocess_many_to_parquet(
    input_paths: list[str | Path],
    output_dir: str | Path,
    smiles_colname: str,
    chem_id_colname: str,
    rows_per_shard: int,
    row_group_size: int,
) -> dict[str, Any]:
    manifests = []
    for input_path in input_paths:
        source_group = Path(input_path).stem
        manifests.append(
            preprocess_to_parquet(
                input_path=input_path,
                output_dir=output_dir,
                smiles_colname=smiles_colname,
                chem_id_colname=chem_id_colname,
                source_group=source_group,
                rows_per_shard=rows_per_shard,
                row_group_size=row_group_size,
            )
        )
    return {
        "source_groups": sorted(item["source_group"] for item in manifests),
        "inputs": manifests,
    }
```

```python
shard = frame[[smiles_colname, chem_id_colname]].rename(
    columns={smiles_colname: "smiles", chem_id_colname: "chem_id"}
)
shard["source_group"] = source_group
shard["source_file"] = str(input_path)
shard["source_row"] = range(row_offset + 1, row_offset + 1 + len(shard))
shard.to_parquet(shard_path, index=False, row_group_size=row_group_size)
```

- [ ] **Step 4: Teach `process_work_unit()` to prefer Parquet shards**

```python
elif unit.source_type == "parquet":
    parquet_file = pq.ParquetFile(unit.source_path)
    for batch in parquet_file.iter_batches(batch_size=chunk_size, columns=["smiles", "chem_id", "source_group", "source_file", "source_row"]):
        frame = batch.to_pandas()
        yield frame.reset_index(drop=True)
```

- [ ] **Step 5: Re-run preprocessing tests**

Run:

```bash
pixi run python -m unittest discover -s test -p 'test_preprocess_screening_input.py' -v
pixi run python -m unittest discover -s test -p 'test_screening_sources.py' -v
```

Expected: PASS, with one Parquet shard tree per CSV source group.

- [ ] **Step 6: Commit**

```bash
git add src/preprocess_screening_input.py src/screening_sources.py test/test_preprocess_screening_input.py
git commit -m "feat: pre-shard multi-source screening inputs to parquet"
```

### Task 3: Build the multiprocessing producer pipeline around a single predictor process

**Files:**
- Create: `src/screening_process_pipeline.py`
- Modify: `src/batch_screening.py`
- Test: `test/test_screening_process_pipeline.py`

- [ ] **Step 1: Write the failing orchestration tests**

```python
def test_resolve_producer_processes_scales_from_cpu_count():
    from src.screening_process_pipeline import resolve_producer_processes

    assert resolve_producer_processes("auto", cpu_count=24, graph_workers=4) >= 6
    assert resolve_producer_processes("auto", cpu_count=96, graph_workers=12) >= 16
```

```python
def test_process_pipeline_uses_one_predictor_and_multiple_producers(tmp_path):
    from src.screening_process_pipeline import ProcessScreenConfig, plan_runtime

    config = ProcessScreenConfig(
        input_paths=["a.parquet", "b.parquet", "c.parquet", "d.parquet"],
        output_dir=str(tmp_path / "out"),
        execution_mode="process",
        producer_processes="auto",
        predict_queue_max_batches="auto",
        result_queue_max_batches="auto",
        batch_checkpoint_size=32,
    )

    runtime = plan_runtime(config, cpu_count=24, gpu_count=1)
    assert runtime.predictor_processes == 1
    assert runtime.producer_processes >= 6
```

- [ ] **Step 2: Run tests to verify failure**

Run:

```bash
pixi run python -m unittest discover -s test -p 'test_screening_process_pipeline.py' -v
```

Expected: FAIL because the process-screening runtime does not exist yet.

- [ ] **Step 3: Implement process topology helpers**

```python
def resolve_producer_processes(value: int | str, *, cpu_count: int, graph_workers: int) -> int:
    if value != "auto":
        return max(1, int(value))
    reserved = max(2, cpu_count // 12)
    available = max(1, cpu_count - reserved - graph_workers)
    return max(4, min(32, available // 2))
```

```python
def resolve_queue_depth(value: int | str, *, producer_processes: int) -> int:
    if value != "auto":
        return max(1, int(value))
    return max(8, producer_processes * 2)
```

```python
@dataclass(frozen=True)
class RuntimePlan:
    predictor_processes: int
    producer_processes: int
    predict_queue_max_batches: int
    result_queue_max_batches: int
```

- [ ] **Step 4: Implement the producer/predictor/writer topology**

```python
ctx = multiprocessing.get_context("spawn")
work_queue = ctx.JoinableQueue(maxsize=producer_processes * 2)
predict_queue = ctx.Queue(maxsize=predict_queue_max_batches)
result_queue = ctx.Queue(maxsize=result_queue_max_batches)
error_queue = ctx.Queue()
```

```python
def producer_main(work_queue, predict_queue, error_queue, batch_checkpoint_size):
    while True:
        item = work_queue.get()
        if item is None:
            break
        for frame in process_work_unit(...):
            for batch in split_frame(frame, batch_checkpoint_size):
                predict_queue.put(batch)
        work_queue.task_done()
```

```python
async def predictor_main(...):
    scheduler = get_scheduler(...)
    while True:
        batch = await next_batch_from_queue(...)
        if batch is None:
            break
        predicted = await scheduler.predict_molecules(...)
        result_queue.put(predicted_batch)
```

- [ ] **Step 5: Re-run the orchestration tests**

Run:

```bash
pixi run python -m unittest discover -s test -p 'test_screening_process_pipeline.py' -v
```

Expected: PASS, with one predictor plan per GPU and auto-scaled producer counts.

- [ ] **Step 6: Commit**

```bash
git add src/screening_process_pipeline.py src/batch_screening.py test/test_screening_process_pipeline.py
git commit -m "feat: add multiprocess producer pipeline for screening"
```

### Task 4: Add batch-level manifests and hit-only writer semantics

**Files:**
- Modify: `src/screening_process_pipeline.py`
- Modify: `src/batch_screening.py`
- Test: `test/test_screening_process_pipeline.py`

- [ ] **Step 1: Write the failing checkpoint tests**

```python
def test_writer_persists_only_broad_spectrum_hits(tmp_path):
    from src.screening_process_pipeline import commit_prediction_batch
    import pandas as pd

    batch = pd.DataFrame(
        [
            {"chem_id": "a", "smiles": "CCO", "broad_spectrum": 1, "ginhib_total": 11, "apscore_total": 0.1},
            {"chem_id": "b", "smiles": "CCN", "broad_spectrum": 0, "ginhib_total": 2, "apscore_total": 1.2},
        ]
    )

    manifest = {}
    commit_prediction_batch(tmp_path, "shard_1", "batch_1", batch, manifest)

    hits = pd.read_parquet(tmp_path / "hits" / "shard_1" / "batch_1.parquet")
    assert list(hits["chem_id"]) == ["a"]
```

```python
def test_rerun_skips_committed_batch(tmp_path):
    from src.screening_process_pipeline import should_skip_batch

    manifest = {"batch_1": {"status": "committed", "output_path": str(tmp_path / "hits" / "x.parquet")}}
    assert should_skip_batch("batch_1", manifest) is True
```

- [ ] **Step 2: Run tests to verify failure**

Run:

```bash
pixi run python -m unittest discover -s test -p 'test_screening_process_pipeline.py' -v
```

Expected: FAIL because the writer/manifest helpers do not exist yet.

- [ ] **Step 3: Implement atomic hit-only commit**

```python
def commit_prediction_batch(output_dir, shard_id, batch_id, predicted_frame, manifest):
    hits = predicted_frame[predicted_frame["broad_spectrum"] == 1].copy()
    target = Path(output_dir) / "hits" / shard_id / f"{batch_id}.parquet"
    tmp_target = target.with_suffix(".parquet.tmp")
    tmp_target.parent.mkdir(parents=True, exist_ok=True)
    hits.to_parquet(tmp_target, index=False)
    tmp_target.replace(target)
    manifest[batch_id] = {
        "status": "committed",
        "attempted_count": len(predicted_frame),
        "hit_count": len(hits),
        "output_path": str(target),
    }
```

```python
def should_skip_batch(batch_id, manifest):
    row = manifest.get(batch_id)
    return bool(row and row.get("status") == "committed" and Path(row["output_path"]).exists())
```

- [ ] **Step 4: Persist batch manifest and run state**

```python
record = {
    "batch_id": batch_id,
    "shard_id": shard_id,
    "source_group": source_group,
    "start_row": start_row,
    "end_row": end_row,
    "status": "committed",
    "attempted_count": attempted_count,
    "hit_count": hit_count,
    "output_path": str(output_path),
    "updated_at": _utc_now(),
    "error": None,
}
```

- [ ] **Step 5: Re-run the checkpoint tests**

Run:

```bash
pixi run python -m unittest discover -s test -p 'test_screening_process_pipeline.py' -v
pixi run python -m unittest discover -s test -p 'test_batch_screening.py' -v
```

Expected: PASS, with hit-only persistence and idempotent rerun semantics.

- [ ] **Step 6: Commit**

```bash
git add src/screening_process_pipeline.py src/batch_screening.py test/test_screening_process_pipeline.py
git commit -m "feat: add batch-level hit-only checkpointing for process screening"
```

### Task 5: Wire process mode into `mole screen`, benchmark it, and document the operational guidance

**Files:**
- Modify: `mole_cli.py`
- Modify: `README.md`
- Modify: `docs/cli_reference.md`
- Modify: `docs/repo_layout.md`
- Create: `docs/process_screening.md`

- [ ] **Step 1: Add a failing CLI integration test**

```python
def test_screen_command_routes_process_mode(tmp_path):
    from mole_cli import main

    rc = main(
        [
            "screen",
            "--input-path", str(tmp_path / "a.csv"),
            "--output-dir", str(tmp_path / "out"),
            "--execution-mode", "process",
        ]
    )
    assert rc == 0
```

- [ ] **Step 2: Run the integration tests to verify failure**

Run:

```bash
pixi run python -m unittest discover -s test -p 'test_cli.py' -v
```

Expected: FAIL because `mole screen` does not dispatch to the process pipeline yet.

- [ ] **Step 3: Dispatch `screen` by execution mode**

```python
if args.execution_mode == "process":
    summary = asyncio_run(
        screen_paths_multiprocess(
            input_paths=args.input_path,
            output_dir=args.output_dir,
            producer_processes=args.producer_processes,
            predict_queue_max_batches=args.predict_queue_max_batches,
            result_queue_max_batches=args.result_queue_max_batches,
            batch_checkpoint_size=args.batch_checkpoint_size,
            rows_per_shard=args.rows_per_shard,
            row_group_size=args.row_group_size,
        )
    )
else:
    summary = asyncio_run(screen_path(...))
```

- [ ] **Step 4: Document the operating model**

```markdown
- One predictor process per GPU.
- Use repeated `--input-path` for multiple per-scaffold CSVs.
- Process mode first pre-shards CSV inputs to Parquet, then launches CPU producer processes.
- Only broad-spectrum hits are written under `hits/`; non-hit molecules are counted but not persisted.
- Resume safety is batch-level, not only shard-level.
```

- [ ] **Step 5: Run the verification set**

Run:

```bash
pixi run python -m unittest discover -s test -p 'test_cli.py' -v
pixi run python -m unittest discover -s test -p 'test_preprocess_screening_input.py' -v
pixi run python -m unittest discover -s test -p 'test_screening_process_pipeline.py' -v
pixi run python -m unittest discover -s test -p 'test_batch_screening.py' -v
```

Expected: PASS, with process mode wired end-to-end.

- [ ] **Step 6: Run a small real benchmark on the four-source pattern**

Run:

```bash
CUDA_VISIBLE_DEVICES=0 pixi run mole screen \
  --input-path data/04.new_predictions/2026-04-21_screening/cache/extracted/2026-04-21_screening/tylosin/scheme_b_fix_pos13/tylosin_scheme_b_unique_products.csv \
  --input-path data/04.new_predictions/2026-04-21_screening/cache/extracted/2026-04-21_screening/tilmicosin/scheme_b_fix_pos13/tilmicosin_scheme_b_unique_products.csv \
  --input-path data/04.new_predictions/2026-04-21_screening/cache/extracted/2026-04-21_screening/tildipirosin/scheme_b_fix_pos13/tildipirosin_scheme_b_unique_products.csv \
  --input-path data/04.new_predictions/2026-04-21_screening/cache/extracted/2026-04-21_screening/tylvalosin/scheme_b_fix_pos13/tylvalosin_scheme_b_unique_products.csv \
  --output-dir data/04.new_predictions/2026-04-21_screening/runs/process_mode_smoke \
  --execution-mode process \
  --producer-processes auto \
  --predict-queue-max-batches auto \
  --result-queue-max-batches auto \
  --batch-checkpoint-size 2048 \
  --classifier-backend pickle \
  --num-graph-workers 4 \
  --graph-batch-size 1024
```

Expected: output directory contains Parquet preparation manifests, batch manifests, and only hit parquet files.

- [ ] **Step 7: Commit**

```bash
git add mole_cli.py README.md docs/cli_reference.md docs/repo_layout.md docs/process_screening.md
git commit -m "docs: document multiprocess hit-only screening workflow"
```

---

### Self-Review

- Spec coverage: the plan covers multi-input CSV ingestion, Parquet pre-sharding, auto-scaled CPU producer processes, single-predictor-per-GPU runtime, batch-level checkpointing, hit-only persistence, and documentation.
- Placeholder scan: no `TODO`, `TBD`, or undefined implementation steps remain.
- Type consistency: the same configuration names are used throughout the plan: `execution_mode`, `producer_processes`, `predict_queue_max_batches`, `result_queue_max_batches`, and `batch_checkpoint_size`.
