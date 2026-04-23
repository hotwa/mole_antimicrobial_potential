from __future__ import annotations

import json
import tempfile
import threading
import time
import unittest
from concurrent.futures import Future
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd

import mole_cli
import src.stream_enumeration_screen as stream_module


class StreamEnumerationIndexTests(unittest.TestCase):
    def test_global_index_round_trip(self) -> None:
        space = stream_module.IndexSpace(scaffold_count=2, pos3_count=3, pos4_count=5, pos12_count=7, pos13_count=11)
        for index in [0, 1, 10, 115, space.total_combinations - 1]:
            parts = stream_module.decode_global_index(index, space)
            self.assertEqual(stream_module.encode_global_index(parts, space), index)

    def test_global_index_uses_scaffold_major_order(self) -> None:
        space = stream_module.IndexSpace(scaffold_count=2, pos3_count=2, pos4_count=2, pos12_count=2, pos13_count=2)
        self.assertEqual(stream_module.decode_global_index(0, space), (0, 0, 0, 0, 0))
        self.assertEqual(stream_module.decode_global_index(1, space), (0, 0, 0, 0, 1))
        self.assertEqual(stream_module.decode_global_index(8, space), (0, 1, 0, 0, 0))
        self.assertEqual(stream_module.decode_global_index(16, space), (1, 0, 0, 0, 0))


class _FakeScheduler:
    def __init__(self, hit_indexes: set[int]) -> None:
        self.hit_indexes = set(hit_indexes)
        self.calls: list[list[str]] = []

    async def predict_molecules(self, molecules, **kwargs):
        self.calls.append([m.chem_id for m in molecules])
        items = []
        for molecule in molecules:
            global_index = int(str(molecule.chem_id).rsplit("__g", 1)[1])
            is_hit = global_index in self.hit_indexes
            items.append(
                {
                    "chem_id": molecule.chem_id,
                    "broad_spectrum": 1 if is_hit else 0,
                    "ginhib_total": 10 if is_hit else 0,
                    "apscore_total": -1.0 if is_hit else 0.0,
                }
            )
        return items

    def runtime_snapshot(self):
        return {"backend": "fake"}


class _InterruptingScheduler:
    async def predict_molecules(self, molecules, **kwargs):
        raise KeyboardInterrupt()

    def runtime_snapshot(self):
        return {"backend": "fake"}


class _FakeMaterializationExecutor:
    def __init__(
        self,
        *,
        max_workers=None,
        space,
        scaffolds,
        ordinary_fragments,
        pos13_fragments,
        delay_by_start: dict[int, float] | None = None,
        on_batch_complete=None,
    ) -> None:
        self._space = space
        self._scaffolds = scaffolds
        self._ordinary_fragments = ordinary_fragments
        self._pos13_fragments = pos13_fragments
        self._delay_by_start = delay_by_start or {}
        self._on_batch_complete = on_batch_complete
        self._threads: list[threading.Thread] = []

    def submit_materialize_batch(self, *, batch_start: int, batch_end: int) -> Future:
        future: Future = Future()

        def _run() -> None:
            try:
                time.sleep(self._delay_by_start.get(batch_start, 0.0))
                rows = stream_module._materialize_batch(
                    start_idx=batch_start,
                    end_idx=batch_end,
                    space=self._space,
                    scaffolds=self._scaffolds,
                    ordinary_fragments=self._ordinary_fragments,
                    pos13_fragments=self._pos13_fragments,
                )
                if self._on_batch_complete is not None:
                    self._on_batch_complete(batch_start=batch_start, batch_end=batch_end, rows=rows)
                future.set_result(rows)
            except BaseException as exc:  # pragma: no cover
                future.set_exception(exc)

        thread = threading.Thread(target=_run, daemon=True)
        thread.start()
        self._threads.append(thread)
        return future

    def shutdown(self, *, wait: bool, cancel_futures: bool) -> None:
        if wait:
            for thread in self._threads:
                thread.join(timeout=1.0)


class StreamEnumerationRunTests(unittest.IsolatedAsyncioTestCase):
    def _write_inputs(self, root: Path) -> dict[str, Path]:
        scaffold = root / "scaffold.smi"
        scaffold.write_text("[*:0]CC([*:1])C([*:2])C[*:3]\n", encoding="utf-8")

        ordinary = root / "ordinary.csv"
        ordinary.write_text(
            "\n".join(
                [
                    "fragment_role,fragment_smiles_plain,fragment_smiles_canonical,heavy_atoms,row_count,fragment_smiles_labeled_example,source_cleavage_positions,source_parent_ml_ids,source_fragment_ids,source_bond_types",
                    "ordinary,*C,*C,1,1,[3*]C,3,ML1,frag1,SINGLE",
                    "ordinary,*N,*N,1,1,[3*]N,3,ML2,frag2,SINGLE",
                ]
            )
            + "\n",
            encoding="utf-8",
        )

        pos13 = root / "pos13.csv"
        pos13.write_text(
            "\n".join(
                [
                    "library_position,fragment_role,fragment_smiles_plain,fragment_smiles_canonical,fragment_smiles_labeled_example,heavy_atoms,sugar_structure_id,source_row_count,source_parent_ml_ids,source_fragment_ids,source_bond_types",
                    "13,sugar,*O,*O,*O,1,POS13_1,1,ML3,frag3,SINGLE",
                    "13,sugar,*F,*F,*F,1,POS13_2,1,ML4,frag4,SINGLE",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        return {"scaffold": scaffold, "ordinary": ordinary, "pos13": pos13}

    def _read_all_hits(self, output_dir: Path) -> pd.DataFrame:
        hits_dir = output_dir / "hits"
        files = sorted(hits_dir.glob("*.parquet"))
        if not files:
            return pd.DataFrame()
        return pd.concat([pd.read_parquet(path) for path in files], ignore_index=True)

    def _read_run_state(self, output_dir: Path) -> dict:
        return json.loads((output_dir / "run_state.json").read_text(encoding="utf-8"))

    async def _run_stream_screen(
        self,
        *,
        output_dir: Path,
        hit_indexes: set[int],
        stop_index: int,
        fail_after_shards: int | None = None,
        scheduler_override=None,
        enumeration_workers: int | str = 2,
        enumeration_prefetch_batches: int | str = 2,
        prediction_batch_size: int = 2,
        shard_size: int = 2,
    ):
        scheduler = scheduler_override or _FakeScheduler(hit_indexes)
        inputs = self._write_inputs(output_dir.parent)
        return await stream_module.stream_enumeration_screen(
            output_dir=output_dir,
            scaffold_file=inputs["scaffold"],
            scaffold_dir=None,
            scaffold_catalog=None,
            ordinary_library=inputs["ordinary"],
            pos13_library=inputs["pos13"],
            run_state_source=None,
            chunk_manifest_source=None,
            start_index=0,
            stop_index=stop_index,
            shard_size=shard_size,
            prediction_batch_size=prediction_batch_size,
            scheduler=scheduler,
            fail_after_shards=fail_after_shards,
            classifier_backend="pickle",
            enumeration_workers=enumeration_workers,
            enumeration_prefetch_batches=enumeration_prefetch_batches,
        )

    async def test_hits_only_persistence(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "run"
            await self._run_stream_screen(output_dir=output_dir, hit_indexes={1, 3}, stop_index=4)
            persisted = self._read_all_hits(output_dir)
            self.assertEqual(sorted(persisted["global_combination_index"].tolist()), [1, 3])

    async def test_resume_skips_completed_shards(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "run"
            with self.assertRaises(RuntimeError):
                await self._run_stream_screen(
                    output_dir=output_dir,
                    hit_indexes={2, 5},
                    stop_index=6,
                    fail_after_shards=1,
                )

            await self._run_stream_screen(output_dir=output_dir, hit_indexes={2, 5}, stop_index=6)
            persisted = self._read_all_hits(output_dir)
            self.assertEqual(sorted(persisted["global_combination_index"].tolist()), [2, 5])

            manifest_lines = (output_dir / "shard_manifest.jsonl").read_text(encoding="utf-8").splitlines()
            completed = [json.loads(line) for line in manifest_lines if json.loads(line)["status"] == "completed"]
            self.assertEqual(len(completed), 3)

    async def test_idempotent_rerun_does_not_duplicate_hits(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "run"
            await self._run_stream_screen(output_dir=output_dir, hit_indexes={0, 2}, stop_index=4)
            await self._run_stream_screen(output_dir=output_dir, hit_indexes={0, 2}, stop_index=4)
            persisted = self._read_all_hits(output_dir)
            self.assertEqual(len(persisted), 2)
            self.assertEqual(len(set(persisted["global_combination_index"].tolist())), 2)

    async def test_keyboard_interrupt_marks_run_interrupted_and_resume_completes(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "run"
            with self.assertRaises(KeyboardInterrupt):
                await self._run_stream_screen(
                    output_dir=output_dir,
                    hit_indexes=set(),
                    stop_index=4,
                    scheduler_override=_InterruptingScheduler(),
                )

            run_state = json.loads((output_dir / "run_state.json").read_text(encoding="utf-8"))
            self.assertEqual(run_state["status"], "interrupted")
            manifest_rows = [
                json.loads(line)
                for line in (output_dir / "shard_manifest.jsonl").read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]
            self.assertEqual(manifest_rows[0]["status"], "interrupted")

            await self._run_stream_screen(output_dir=output_dir, hit_indexes={1}, stop_index=4)
            persisted = self._read_all_hits(output_dir)
            self.assertEqual(persisted["global_combination_index"].tolist(), [1])

    async def test_enumerated_rows_include_required_output_fields(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "run"
            await self._run_stream_screen(output_dir=output_dir, hit_indexes={0}, stop_index=1)
            row = self._read_all_hits(output_dir).iloc[0].to_dict()
            self.assertIn("chem_id", row)
            self.assertIn("smiles", row)
            self.assertIn("broad_spectrum", row)
            self.assertIn("ginhib_total", row)
            self.assertIn("apscore_total", row)
            self.assertIn("global_combination_index", row)
            self.assertIn("shard_id", row)
            self.assertIn("scaffold_slug", row)

    async def test_enumeration_prefetches_multiple_batches_while_prediction_waits(self) -> None:
        class DelayedScheduler:
            def __init__(self) -> None:
                self.predict_started = threading.Event()
                self.release_prediction = threading.Event()
                self.materialized_batches = 0
                self.max_materialized_while_predicting = 0

            async def predict_molecules(self, molecules, **kwargs):
                self.predict_started.set()
                deadline = time.monotonic() + 1.5
                while time.monotonic() < deadline:
                    self.max_materialized_while_predicting = max(
                        self.max_materialized_while_predicting,
                        self.materialized_batches,
                    )
                    if self.materialized_batches >= 2:
                        break
                    await stream_module.asyncio.sleep(0.01)
                self.release_prediction.set()
                return [
                    {
                        "chem_id": molecule.chem_id,
                        "broad_spectrum": 0,
                        "ginhib_total": 0,
                        "apscore_total": 0.0,
                    }
                    for molecule in molecules
                ]

            def runtime_snapshot(self):
                return {"backend": "fake"}

        scheduler = DelayedScheduler()

        def fake_executor_factory(**kwargs):
            return _FakeMaterializationExecutor(
                **kwargs,
                delay_by_start={0: 0.02, 1: 0.02, 2: 0.02, 3: 0.02},
                on_batch_complete=lambda **_: _record_materialized_batch(),
            )

        def _record_materialized_batch() -> None:
            scheduler.materialized_batches += 1
            if scheduler.predict_started.is_set() and not scheduler.release_prediction.is_set():
                scheduler.max_materialized_while_predicting = max(
                    scheduler.max_materialized_while_predicting,
                    scheduler.materialized_batches,
                )

        with tempfile.TemporaryDirectory() as tmpdir, patch.object(
            stream_module,
            "_create_materialization_executor",
            side_effect=fake_executor_factory,
        ):
            output_dir = Path(tmpdir) / "run"
            inputs = self._write_inputs(Path(tmpdir))
            await stream_module.stream_enumeration_screen(
                output_dir=output_dir,
                scaffold_file=inputs["scaffold"],
                scaffold_dir=None,
                scaffold_catalog=None,
                ordinary_library=inputs["ordinary"],
                pos13_library=inputs["pos13"],
                start_index=0,
                stop_index=4,
                shard_size=4,
                prediction_batch_size=1,
                scheduler=scheduler,
                classifier_backend="pickle",
                enumeration_workers=2,
                enumeration_prefetch_batches=2,
            )

        self.assertGreaterEqual(scheduler.max_materialized_while_predicting, 2)

    async def test_enumeration_prefetch_respects_bound_and_preserves_order(self) -> None:
        scheduler = _FakeScheduler(set())
        def fake_executor_factory(**kwargs):
            return _FakeMaterializationExecutor(
                **kwargs,
                delay_by_start={0: 0.05, 1: 0.01, 2: 0.01, 3: 0.01},
            )

        with tempfile.TemporaryDirectory() as tmpdir, patch.object(
            stream_module,
            "_create_materialization_executor",
            side_effect=fake_executor_factory,
        ):
            output_dir = Path(tmpdir) / "run"
            inputs = self._write_inputs(Path(tmpdir))
            await stream_module.stream_enumeration_screen(
                output_dir=output_dir,
                scaffold_file=inputs["scaffold"],
                scaffold_dir=None,
                scaffold_catalog=None,
                ordinary_library=inputs["ordinary"],
                pos13_library=inputs["pos13"],
                start_index=0,
                stop_index=4,
                shard_size=4,
                prediction_batch_size=1,
                scheduler=scheduler,
                classifier_backend="pickle",
                enumeration_workers=2,
                enumeration_prefetch_batches=2,
            )

            run_state = self._read_run_state(output_dir)
            self.assertLessEqual(run_state["runtime"]["max_pending_materialized_batches"], 2)
            self.assertEqual(
                [call[0] for call in scheduler.calls],
                ["scaffold__g0", "scaffold__g1", "scaffold__g2", "scaffold__g3"],
            )

    async def test_run_state_tracks_runtime_progress_metrics(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "run"
            await self._run_stream_screen(output_dir=output_dir, hit_indexes={1, 3}, stop_index=4)
            run_state = self._read_run_state(output_dir)
            runtime = run_state["runtime"]
            self.assertEqual(runtime["attempted_count"], 4)
            self.assertEqual(runtime["hit_count"], 2)
            self.assertIn("pending_materialized_batches", runtime)
            self.assertIn("prefetch_depth", runtime)
            self.assertIn("scheduler", runtime)

    async def test_auto_enumeration_settings_use_process_materializer(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "run"
            await self._run_stream_screen(
                output_dir=output_dir,
                hit_indexes={0},
                stop_index=2,
                enumeration_workers="auto",
                enumeration_prefetch_batches="auto",
                prediction_batch_size=1,
                shard_size=2,
            )
            run_state = self._read_run_state(output_dir)
            runtime = run_state["runtime"]
            self.assertEqual(runtime["materialization_mode"], "process")
            self.assertGreaterEqual(runtime["resolved_enumeration_workers"], 1)
            self.assertGreaterEqual(runtime["prefetch_depth"]["configured"], 1)


class StreamEnumerationCliTests(unittest.TestCase):
    def test_parser_exposes_stream_enumeration_screen_subcommand(self) -> None:
        parser = mole_cli._build_parser()
        args = parser.parse_args(
            [
                "stream-enumeration-screen",
                "--output-dir",
                "out",
                "--ordinary-library",
                "ordinary.csv",
                "--pos13-library",
                "pos13.csv",
            ]
        )
        self.assertEqual(args.command, "stream-enumeration-screen")
        self.assertEqual(args.output_dir, "out")

    def test_parser_accepts_enumeration_tuning_flags(self) -> None:
        parser = mole_cli._build_parser()
        args = parser.parse_args(
            [
                "stream-enumeration-screen",
                "--output-dir",
                "out",
                "--ordinary-library",
                "ordinary.csv",
                "--pos13-library",
                "pos13.csv",
                "--enumeration-workers",
                "3",
                "--enumeration-prefetch-batches",
                "4",
                "--classifier-workers",
                "2",
                "--classifier-inflight-batches",
                "5",
            ]
        )
        self.assertEqual(args.enumeration_workers, "3")
        self.assertEqual(args.enumeration_prefetch_batches, "4")
        self.assertEqual(args.classifier_workers, "2")
        self.assertEqual(args.classifier_inflight_batches, "5")


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
