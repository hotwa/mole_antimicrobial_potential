# Stream Enumeration Screen Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a resumable `mole stream-enumeration-screen` CLI that enumerates molecules directly from scaffold and fragment libraries, predicts broad-spectrum activity in batches, and persists only hit shards with atomic commit semantics.

**Architecture:** Introduce a dedicated streaming enumeration module instead of extending `mole screen`. The new module owns deterministic global-index math, scaffold/fragment loading, shard checkpoint state, hit-only writing, and resume/idempotency rules, while reusing the existing scheduler/predictor stack for aggregate prediction.

**Tech Stack:** Python 3.10, pandas, pyarrow, RDKit, existing prediction scheduler, unittest, pixi tasks

---

### Task 1: Define CLI Surface And Test The Parser

**Files:**
- Modify: `mole_cli.py`
- Test: `test/test_cli.py`

- [ ] **Step 1: Write the failing parser test**

```python
    def test_parser_exposes_stream_enumeration_screen_subcommand(self) -> None:
        parser = mole_cli._build_parser()
        args = parser.parse_args(["stream-enumeration-screen", "--output-dir", "out"])
        self.assertEqual(args.command, "stream-enumeration-screen")
        self.assertEqual(args.output_dir, "out")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run python -m unittest test.test_cli.TestCli.test_parser_exposes_stream_enumeration_screen_subcommand -v`
Expected: FAIL because the parser does not define `stream-enumeration-screen`

- [ ] **Step 3: Add the parser entry and command dispatch**

Add a new subcommand in `mole_cli.py` with explicit scaffold input support:

```python
    stream = subparsers.add_parser(
        "stream-enumeration-screen",
        help="Enumerate scaffold/fragments on the fly and persist broad-spectrum hit shards",
    )
```

Include `--output-dir`, `--scaffold-file`, `--scaffold-dir`, `--scaffold-catalog`,
`--ordinary-library`, `--pos13-library`, `--run-state`, `--chunk-manifest`,
`--start-index`, `--stop-index`, `--shard-size`, `--prediction-batch-size`,
the existing prediction knobs, and defaults anchored to the task inputs.

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run python -m unittest test.test_cli.TestCli.test_parser_exposes_stream_enumeration_screen_subcommand -v`
Expected: PASS

### Task 2: Lock Deterministic Index Mapping With Red-Green Tests

**Files:**
- Create: `src/stream_enumeration_screen.py`
- Test: `test/test_stream_enumeration_screen.py`

- [ ] **Step 1: Write failing index mapping tests**

```python
    def test_global_index_round_trip(self) -> None:
        space = stream_module.IndexSpace(scaffold_count=2, pos3_count=3, pos4_count=5, pos12_count=7, pos13_count=11)
        for index in [0, 1, 10, 115, space.total_combinations - 1]:
            parts = stream_module.decode_global_index(index, space)
            self.assertEqual(stream_module.encode_global_index(parts, space), index)

    def test_global_index_uses_scaffold_major_order(self) -> None:
        space = stream_module.IndexSpace(scaffold_count=2, pos3_count=2, pos4_count=2, pos12_count=2, pos13_count=2)
        self.assertEqual(stream_module.decode_global_index(0, space), (0, 0, 0, 0, 0))
        self.assertEqual(stream_module.decode_global_index(1, space), (0, 0, 0, 0, 1))
        self.assertEqual(stream_module.decode_global_index(8, space), (1, 0, 0, 0, 0))
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pixi run python -m unittest test.test_stream_enumeration_screen.StreamEnumerationIndexTests -v`
Expected: FAIL because `src/stream_enumeration_screen.py` does not exist yet

- [ ] **Step 3: Implement the minimal index helpers**

Add an `IndexSpace` dataclass plus `encode_global_index()` and `decode_global_index()` with ordering:
`(scaffold_idx, pos3_idx, pos4_idx, pos12_idx, pos13_idx)`.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pixi run python -m unittest test.test_stream_enumeration_screen.StreamEnumerationIndexTests -v`
Expected: PASS

### Task 3: Add Shard Resume And Hit-Only Persistence Tests First

**Files:**
- Modify: `test/test_stream_enumeration_screen.py`
- Modify: `src/stream_enumeration_screen.py`

- [ ] **Step 1: Write failing behavior tests**

Cover these behaviors with a fake scheduler and tiny fragment libraries:

```python
    def test_hits_only_persistence(self) -> None:
        result = self.run_stream_screen(hit_indexes={1, 3}, stop_index=4)
        persisted = self.read_all_hits(result.output_dir)
        self.assertEqual(sorted(persisted["global_combination_index"].tolist()), [1, 3])

    def test_resume_skips_completed_shards(self) -> None:
        first = self.run_stream_screen(hit_indexes={2, 5}, stop_index=6, fail_after_shards=1)
        resumed = self.run_stream_screen(hit_indexes={2, 5}, stop_index=6)
        persisted = self.read_all_hits(resumed.output_dir)
        self.assertEqual(sorted(persisted["global_combination_index"].tolist()), [2, 5])

    def test_idempotent_rerun_does_not_duplicate_hits(self) -> None:
        self.run_stream_screen(hit_indexes={0, 2}, stop_index=4)
        rerun = self.run_stream_screen(hit_indexes={0, 2}, stop_index=4)
        persisted = self.read_all_hits(rerun.output_dir)
        self.assertEqual(len(persisted), 2)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pixi run python -m unittest test.test_stream_enumeration_screen.StreamEnumerationRunTests -v`
Expected: FAIL because shard execution, manifest writing, and rerun behavior are not implemented

- [ ] **Step 3: Implement minimal shard execution**

Implement:
- scaffold and library loading
- deterministic shard planning from `start_index`, `stop_index`, and `shard_size`
- hit-only parquet output per shard
- atomic temp write then rename
- `run_state.json` parameter snapshot
- `shard_manifest.jsonl` rewriting with one record per shard
- completed-shard skipping on rerun

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pixi run python -m unittest test.test_stream_enumeration_screen.StreamEnumerationRunTests -v`
Expected: PASS

### Task 4: Implement Chemistry Assembly And Scheduler Integration

**Files:**
- Modify: `src/stream_enumeration_screen.py`
- Modify: `mole_cli.py`

- [ ] **Step 1: Write a failing assembly/integration test**

```python
    def test_enumerated_rows_include_required_output_fields(self) -> None:
        result = self.run_stream_screen(hit_indexes={0}, stop_index=1)
        row = self.read_all_hits(result.output_dir).iloc[0].to_dict()
        self.assertIn("chem_id", row)
        self.assertIn("smiles", row)
        self.assertIn("broad_spectrum", row)
        self.assertIn("ginhib_total", row)
        self.assertIn("apscore_total", row)
        self.assertIn("global_combination_index", row)
        self.assertIn("shard_id", row)
        self.assertIn("scaffold_slug", row)
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `pixi run python -m unittest test.test_stream_enumeration_screen.StreamEnumerationRunTests.test_enumerated_rows_include_required_output_fields -v`
Expected: FAIL until the output schema is complete

- [ ] **Step 3: Implement the chemistry and prediction path**

Use RDKit to:
- read the scaffold SMILES
- find labeled scaffold dummy atoms
- attach each single-dummy fragment in the fixed position order
- sanitize and canonicalize the final molecule SMILES

Then call the scheduler in aggregate mode and filter rows with:

```python
(broad_spectrum == 1) | (ginhib_total >= min_nkill)
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `pixi run python -m unittest test.test_stream_enumeration_screen.StreamEnumerationRunTests.test_enumerated_rows_include_required_output_fields -v`
Expected: PASS

### Task 5: Run Focused Verification And A Pixi Smoke Run

**Files:**
- Modify: `pixi.toml`
- Optionally create: `test/test_stream_enumeration_screen.py`

- [ ] **Step 1: Add a pixi task for the new test file**

Add:

```toml
test-stream-enumeration-screen = { cmd = "python -m unittest discover -s test -p 'test_stream_enumeration_screen.py' -v", depends-on = ["py-compile-score"] }
```

- [ ] **Step 2: Run the focused test suite**

Run: `pixi run test-stream-enumeration-screen`
Expected: PASS

- [ ] **Step 3: Run CLI parser regression**

Run: `pixi run test-cli`
Expected: PASS

- [ ] **Step 4: Run a smoke command with an interrupt and resume**

Run a bounded test command against the provided libraries and default scaffold, for example:

```bash
pixi run mole stream-enumeration-screen \
  --output-dir data/05.stream_tasks/smoke_runs/stream_screen_demo \
  --run-state data/05.stream_tasks/thought2_stream_screen_2026-04-23/thought2_stream_screen_2026-04-23/run_state.json \
  --ordinary-library data/05.stream_tasks/thought2_stream_screen_2026-04-23/thought2_stream_screen_2026-04-23/validation_output/thought2_enumeration/input_libraries/shared_ordinary_library.csv \
  --pos13-library data/05.stream_tasks/thought2_stream_screen_2026-04-23/thought2_stream_screen_2026-04-23/validation_output/thought2_enumeration/input_libraries/pos13_sugar_library.csv \
  --chunk-manifest data/05.stream_tasks/thought2_stream_screen_2026-04-23/thought2_stream_screen_2026-04-23/validation_output/thought2_enumeration/chunk_manifest.csv \
  --stop-index 100000 \
  --shard-size 20000 \
  --prediction-batch-size 512 \
  --classifier-backend pickle
```

Interrupt once during execution, rerun the same command, and verify:
- completed shards are skipped
- the final hit row count is stable across reruns
- no duplicate `global_combination_index` values exist in persisted hits

