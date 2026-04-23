# Stream Enumeration Screening Runbook

This runbook is for `mole stream-enumeration-screen`, which streams
scaffold/fragment combinations directly into MolE/XGBoost prediction and writes
only broad-spectrum hit shards. It is intended for very large combinatorial
spaces where writing full molecule tables or `predictions_all.tsv` is not
acceptable.

## Scope

The current Thought2 snapshot has one scaffold-sized combination space:

- `pos3`: 1123 fragments
- `pos4`: 1123 fragments
- `pos12`: 1123 fragments
- `pos13`: 674 fragments
- per-scaffold total: `954,551,062,358`

The deterministic index order is:

```text
global_combination_index <-> (scaffold_idx, pos3_idx, pos4_idx, pos12_idx, pos13_idx)
```

`scaffold_idx` is the slowest-changing dimension. `pos13_idx` is the
fastest-changing dimension.

## Input Paths

Use these snapshot inputs for the 2026-04-23 Thought2 run:

```bash
TASK_ROOT=data/05.stream_tasks/thought2_stream_screen_2026-04-23/thought2_stream_screen_2026-04-23
RUN_STATE="$TASK_ROOT/run_state.json"
ORDINARY_LIBRARY="$TASK_ROOT/validation_output/thought2_enumeration/input_libraries/shared_ordinary_library.csv"
POS13_LIBRARY="$TASK_ROOT/validation_output/thought2_enumeration/input_libraries/pos13_sugar_library.csv"
CHUNK_MANIFEST="$TASK_ROOT/validation_output/thought2_enumeration/chunk_manifest.csv"
```

The default scaffold is:

```bash
SCAFFOLD_FILE=workflows/reinvent4/inputs/scaffolds/mother_scaffold.template.smi
```

If a future run has multiple scaffolds, use either `--scaffold-dir` or
`--scaffold-catalog`; do not rely on `run_state.json` `reference_slugs` as the
source of structure.

## Recommended Parameters

For the current RTX 2080 Ti 22GB machine:

| Parameter | Recommended value | Reason |
| --- | --- | --- |
| `--classifier-backend` | `pickle` | Closest to original MolE/XGBoost path and faster than Timber in local benchmarks |
| `--num-graph-workers` | `0` | Best local steady-state benchmark on this machine; avoids worker overhead |
| `--graph-batch-size` | `1024` | Stable default for MolE graph/forward batching |
| `--prediction-batch-size` | `512` | Keeps stream enumeration memory bounded while giving the scheduler enough rows |
| `--shard-size` | `1000000` | Recovery granularity of roughly one million combinations per shard |
| checkpoint frequency | one shard | State commits after every completed shard |

For RTX 5090 32GB, start with:

| Parameter | Starting value |
| --- | --- |
| `--classifier-backend` | `pickle` |
| `--num-graph-workers` | `0` |
| `--graph-batch-size` | `2048` |
| `--prediction-batch-size` | `1024` |
| `--shard-size` | `1000000` to `5000000` |

Re-benchmark `--num-graph-workers 0/2/4/8` on the new machine before production.
This parameter changes throughput and CPU contention, not prediction semantics
for a fixed model/backend.

## Dry-Run With Manual Interrupt And Resume

Start from a clean demo output directory:

```bash
rm -rf data/05.stream_tasks/smoke_runs/stream_screen_demo_resume
```

Launch a bounded smoke run:

```bash
pixi run mole stream-enumeration-screen \
  --output-dir data/05.stream_tasks/smoke_runs/stream_screen_demo_resume \
  --run-state "$RUN_STATE" \
  --ordinary-library "$ORDINARY_LIBRARY" \
  --pos13-library "$POS13_LIBRARY" \
  --chunk-manifest "$CHUNK_MANIFEST" \
  --stop-index 4000 \
  --shard-size 1000 \
  --prediction-batch-size 256 \
  --classifier-backend pickle \
  --num-graph-workers 0
```

Manually interrupt with `Ctrl-C` after at least one shard commits. Confirm the
interrupted state:

```bash
pixi run python - <<'PY'
from pathlib import Path
import json

base = Path("data/05.stream_tasks/smoke_runs/stream_screen_demo_resume")
print(json.loads((base / "run_state.json").read_text())["status"])
for line in (base / "shard_manifest.jsonl").read_text().splitlines():
    row = json.loads(line)
    print(row["shard_id"], row["status"], row["last_committed_idx"], row["attempted_count"], row["hit_count"])
PY
```

Resume by repeating the same command:

```bash
pixi run mole stream-enumeration-screen \
  --output-dir data/05.stream_tasks/smoke_runs/stream_screen_demo_resume \
  --run-state "$RUN_STATE" \
  --ordinary-library "$ORDINARY_LIBRARY" \
  --pos13-library "$POS13_LIBRARY" \
  --chunk-manifest "$CHUNK_MANIFEST" \
  --stop-index 4000 \
  --shard-size 1000 \
  --prediction-batch-size 256 \
  --classifier-backend pickle \
  --num-graph-workers 0
```

Check completion and duplicate safety:

```bash
pixi run python - <<'PY'
from pathlib import Path
import json
import pandas as pd

base = Path("data/05.stream_tasks/smoke_runs/stream_screen_demo_resume")
state = json.loads((base / "run_state.json").read_text())
manifest = [json.loads(line) for line in (base / "shard_manifest.jsonl").read_text().splitlines() if line.strip()]
files = sorted((base / "hits").glob("*.parquet"))
hits = pd.concat([pd.read_parquet(path) for path in files], ignore_index=True) if files else pd.DataFrame()

print("status", state["status"])
print("completed_shards", state["progress"]["completed_shards"])
print("attempted_count", state["progress"]["attempted_count"])
print("hit_count", state["progress"]["hit_count"])
print("manifest_status_counts", {s: sum(row["status"] == s for row in manifest) for s in sorted({row["status"] for row in manifest})})
print("hit_rows", len(hits))
print("unique_global_combination_index", hits["global_combination_index"].nunique() if not hits.empty else 0)
print("duplicate_global_combination_index", 0 if hits.empty else len(hits) - hits["global_combination_index"].nunique())
PY
```

Observed local smoke result on 2026-04-23:

```text
completed_shards: 4
attempted_count: 4000
hit_count: 3954
duplicate_global_combination_index: 0
resume elapsed for remaining 3000 combinations: 35.35 s
```

## Production Index Splits

Total combinations:

```text
954551062358
```

Use one independent output directory per worker. Do not let multiple workers
write the same `output-dir`.

### Four GPUs

```bash
CUDA_VISIBLE_DEVICES=0 pixi run mole stream-enumeration-screen \
  --output-dir data/05.stream_tasks/production_runs/thought2_gpu0 \
  --run-state "$RUN_STATE" \
  --ordinary-library "$ORDINARY_LIBRARY" \
  --pos13-library "$POS13_LIBRARY" \
  --chunk-manifest "$CHUNK_MANIFEST" \
  --start-index 0 \
  --stop-index 238637765589 \
  --shard-size 1000000 \
  --prediction-batch-size 512 \
  --classifier-backend pickle \
  --num-graph-workers 0

CUDA_VISIBLE_DEVICES=1 pixi run mole stream-enumeration-screen \
  --output-dir data/05.stream_tasks/production_runs/thought2_gpu1 \
  --run-state "$RUN_STATE" \
  --ordinary-library "$ORDINARY_LIBRARY" \
  --pos13-library "$POS13_LIBRARY" \
  --chunk-manifest "$CHUNK_MANIFEST" \
  --start-index 238637765589 \
  --stop-index 477275531179 \
  --shard-size 1000000 \
  --prediction-batch-size 512 \
  --classifier-backend pickle \
  --num-graph-workers 0

CUDA_VISIBLE_DEVICES=2 pixi run mole stream-enumeration-screen \
  --output-dir data/05.stream_tasks/production_runs/thought2_gpu2 \
  --run-state "$RUN_STATE" \
  --ordinary-library "$ORDINARY_LIBRARY" \
  --pos13-library "$POS13_LIBRARY" \
  --chunk-manifest "$CHUNK_MANIFEST" \
  --start-index 477275531179 \
  --stop-index 715913296768 \
  --shard-size 1000000 \
  --prediction-batch-size 512 \
  --classifier-backend pickle \
  --num-graph-workers 0

CUDA_VISIBLE_DEVICES=3 pixi run mole stream-enumeration-screen \
  --output-dir data/05.stream_tasks/production_runs/thought2_gpu3 \
  --run-state "$RUN_STATE" \
  --ordinary-library "$ORDINARY_LIBRARY" \
  --pos13-library "$POS13_LIBRARY" \
  --chunk-manifest "$CHUNK_MANIFEST" \
  --start-index 715913296768 \
  --stop-index 954551062358 \
  --shard-size 1000000 \
  --prediction-batch-size 512 \
  --classifier-backend pickle \
  --num-graph-workers 0
```

### Eight GPUs

```text
gpu0: start=0            stop=119318882794
gpu1: start=119318882794 stop=238637765589
gpu2: start=238637765589 stop=357956648384
gpu3: start=357956648384 stop=477275531179
gpu4: start=477275531179 stop=596594413973
gpu5: start=596594413973 stop=715913296768
gpu6: start=715913296768 stop=835232179563
gpu7: start=835232179563 stop=954551062358
```

Use the same command shape as the four-GPU section and change
`CUDA_VISIBLE_DEVICES`, `--output-dir`, `--start-index`, and `--stop-index`.

## Real-Time Progress

Run this command while workers are active:

```bash
pixi run python - <<'PY'
from pathlib import Path
import json
import time

TOTAL = 954_551_062_358
roots = sorted(Path("data/05.stream_tasks/production_runs").glob("thought2_gpu*"))
now = time.time()

completed = 0
hits = 0
attempted = 0
last_update = None

for root in roots:
    manifest = root / "shard_manifest.jsonl"
    if not manifest.exists():
        continue
    for line in manifest.read_text().splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        if row["status"] == "completed":
            completed += int(row["end_idx"]) - int(row["start_idx"])
            attempted += int(row.get("attempted_count", 0))
            hits += int(row.get("hit_count", 0))
        updated = row.get("updated_at")
        if updated:
            last_update = max(last_update or updated, updated)

rate_hint_per_minute = 5092
remaining = max(TOTAL - completed, 0)
eta_minutes_single_equivalent = remaining / rate_hint_per_minute if rate_hint_per_minute else None

print(f"workers_seen={len(roots)}")
print(f"completed_combinations={completed}")
print(f"completion_fraction={completed / TOTAL:.12%}")
print(f"attempted_count={attempted}")
print(f"hit_count={hits}")
print(f"hit_rate={hits / attempted:.6%}" if attempted else "hit_rate=NA")
print(f"remaining_combinations={remaining}")
print(f"eta_days_at_5092_per_minute_single_worker={eta_minutes_single_equivalent / 1440:.2f}" if eta_minutes_single_equivalent else "eta_days=NA")
print(f"last_manifest_update={last_update}")
PY
```

For accurate ETA on the target machine, replace `rate_hint_per_minute` with the
observed per-worker throughput from the first production hour.

## Failure Recovery SOP

### OOM

1. Stop the affected worker.
2. Reduce `--prediction-batch-size` by half.
3. Keep the same `--output-dir`, `--start-index`, and `--stop-index`.
4. Rerun the same worker command.

Completed shards are skipped. The interrupted or failed shard is recomputed from
its `start_idx`.

### Worker crash

1. Inspect `run_state.json` and `shard_manifest.jsonl` in that worker output dir.
2. If the run is `failed` or `interrupted`, rerun the identical worker command.
3. Do not delete hit parquet files unless you intentionally want to rerun
   completed shards.

### Machine reboot

1. Confirm all output directories still exist under
   `data/05.stream_tasks/production_runs/`.
2. Relaunch each worker with the same `CUDA_VISIBLE_DEVICES`, `--output-dir`,
   `--start-index`, and `--stop-index`.
3. Use the progress command to verify completed shard counts resume increasing.

### Parameter mismatch

The CLI refuses to resume if the new command's parameter snapshot does not match
the existing `run_state.json`.

If the mismatch is accidental, rerun the original command. If the mismatch is
intentional, create a new `--output-dir` instead of mutating the existing run.

## Merging Hits

The production command writes only hit parquet shards. Merge all worker hits
after completion:

```bash
pixi run python - <<'PY'
from pathlib import Path
import pandas as pd

roots = sorted(Path("data/05.stream_tasks/production_runs").glob("thought2_gpu*"))
files = []
for root in roots:
    files.extend(sorted((root / "hits").glob("*.parquet")))

out = Path("data/05.stream_tasks/production_runs/thought2_broad_spectrum_hits.parquet")
out.parent.mkdir(parents=True, exist_ok=True)

if files:
    frame = pd.concat([pd.read_parquet(path) for path in files], ignore_index=True)
    frame = frame.drop_duplicates(subset=["global_combination_index"]).sort_values(
        ["ginhib_total", "apscore_total", "global_combination_index"],
        ascending=[False, True, True],
    )
    frame.to_parquet(out, index=False)
    print({"hit_rows": len(frame), "output": str(out)})
else:
    print({"hit_rows": 0, "output": str(out), "note": "no hit shards found"})
PY
```

