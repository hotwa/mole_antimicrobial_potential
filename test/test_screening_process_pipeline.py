from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from src.screening_planner import WorkUnit


class ScreeningProcessPipelineTestCase(unittest.TestCase):
    def test_screen_paths_multiprocess_prepares_inputs_and_routes_work_units(self) -> None:
        with self._temp_directory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            input_a = tmp_path / "tylosin.csv"
            input_b = tmp_path / "tilmicosin.csv"
            input_a.write_text("product_smiles_canonical,example_combo_id\nCCO,t1\n", encoding="utf-8")
            input_b.write_text("product_smiles_canonical,example_combo_id\nCCN,m1\n", encoding="utf-8")
            output_dir = tmp_path / "process_run"

            async def _run():
                from src.screening_process_pipeline import screen_paths_multiprocess

                def fake_runner(**kwargs):
                    work_items = kwargs["work_items"]
                    runtime = kwargs["runtime"]
                    self.assertEqual(len(work_items), 2)
                    self.assertTrue(all(item.work_unit.source_type == "parquet" for item in work_items))
                    self.assertTrue(all(str(output_dir / "prepared") in item.work_unit.source_path for item in work_items))
                    self.assertEqual(runtime.predictor_processes, 1)
                    self.assertGreaterEqual(runtime.producer_processes, 1)
                    self.assertEqual(kwargs["classifier_workers"], "2")
                    self.assertEqual(kwargs["classifier_inflight_batches"], "5")

                    batch_manifest_path = output_dir / "batch_manifest.jsonl"
                    run_state_path = output_dir / "run_state.json"
                    hits_path = output_dir / "hits" / "tilmicosin" / "batch_000000.parquet"
                    hits_path.parent.mkdir(parents=True, exist_ok=True)
                    pd.DataFrame(
                        [{"chem_id": "m1", "smiles": "CCN", "broad_spectrum": 1}]
                    ).to_parquet(hits_path, index=False)
                    batch_manifest_path.write_text(
                        json.dumps(
                            {
                                "batch_id": "tilmicosin:batch_000000",
                                "shard_id": "tilmicosin",
                                "source_group": "tilmicosin",
                                "start_row": 1,
                                "end_row": 1,
                                "status": "committed",
                                "attempted_count": 1,
                                "hit_count": 1,
                                "output_path": str(hits_path),
                                "updated_at": "2026-04-23T00:00:00+00:00",
                                "error": None,
                            }
                        )
                        + "\n",
                        encoding="utf-8",
                    )
                    run_state_path.write_text(
                        json.dumps(
                            {
                                "status": "completed",
                                "updated_at": "2026-04-23T00:00:00+00:00",
                                "attempted_count": 1,
                                "hit_count": 1,
                                "output_path": str(hits_path),
                                "progress": {
                                    "committed_batches": 1,
                                    "attempted_count": 1,
                                    "hit_count": 1,
                                },
                                "error": None,
                            }
                        )
                        + "\n",
                        encoding="utf-8",
                    )
                    return {
                        "runtime": runtime,
                        "batch_manifest_path": str(batch_manifest_path),
                        "run_state_path": str(run_state_path),
                        "hits_dir": str(output_dir / "hits"),
                    }

                with mock.patch("src.screening_process_pipeline.run_process_pipeline", side_effect=fake_runner):
                    return await screen_paths_multiprocess(
                        input_paths=[input_a, input_b],
                        output_dir=output_dir,
                        smiles_colname="product_smiles_canonical",
                        chem_id_colname="example_combo_id",
                        execution_mode="process",
                        batch_checkpoint_size=64,
                        classifier_workers="2",
                        classifier_inflight_batches="5",
                    )

            summary = self._run_async(_run())

            prepared_manifest_path = output_dir / "prepared_manifest.json"
            per_input_manifest_dir = output_dir / "prepared_manifests"

            self.assertEqual(summary["execution_mode"], "process")
            self.assertEqual(summary["source_groups"], ["tilmicosin", "tylosin"])
            self.assertEqual(summary["work_unit_count"], 2)
            self.assertEqual(summary["prepared_manifest_path"], str(prepared_manifest_path))
            self.assertEqual(len(summary["prepared_manifest_paths"]), 2)
            self.assertTrue(prepared_manifest_path.exists())
            self.assertTrue((per_input_manifest_dir / "tilmicosin.json").exists())
            self.assertTrue((per_input_manifest_dir / "tylosin.json").exists())
            self.assertTrue((output_dir / "prepared" / "tilmicosin" / "shard_0001.parquet").exists())
            self.assertTrue((output_dir / "prepared" / "tylosin" / "shard_0001.parquet").exists())
            self.assertTrue((output_dir / "batch_manifest.jsonl").exists())
            self.assertTrue((output_dir / "run_state.json").exists())
            self.assertTrue((output_dir / "hits" / "tilmicosin" / "batch_000000.parquet").exists())

    def test_commit_prediction_batch_persists_only_broad_spectrum_hits(self) -> None:
        from src.screening_process_pipeline import (
            commit_prediction_batch,
            load_batch_manifest,
        )

        with self.subTest("hit_only_commit"):
            with unittest.mock.patch("src.screening_process_pipeline._utc_now", return_value="2026-04-23T00:00:00+00:00"):
                with self._temp_directory() as tmp_dir:
                    tmp_path = Path(tmp_dir)
                    manifest_path = tmp_path / "batch_manifest.jsonl"
                    manifest = {}
                    predicted_frame = pd.DataFrame(
                        [
                            {
                                "chem_id": "hit-1",
                                "smiles": "CCO",
                                "broad_spectrum": 1,
                                "ginhib_total": 12,
                                "apscore_total": 0.1,
                            },
                            {
                                "chem_id": "miss-1",
                                "smiles": "CCN",
                                "broad_spectrum": 0,
                                "ginhib_total": 3,
                                "apscore_total": 1.5,
                            },
                        ]
                    )

                    record = commit_prediction_batch(
                        output_dir=tmp_path,
                        shard_id="source_a",
                        batch_id="batch_000000",
                        predicted_frame=predicted_frame,
                        manifest_records=manifest,
                        manifest_path=manifest_path,
                    )

                    hits_path = tmp_path / "hits" / "source_a" / "batch_000000.parquet"
                    hits = pd.read_parquet(hits_path)

                    self.assertEqual(list(hits["chem_id"]), ["hit-1"])
                    self.assertEqual(record["attempted_count"], 2)
                    self.assertEqual(record["hit_count"], 1)
                    self.assertEqual(record["output_path"], str(hits_path))
                    self.assertEqual(record["status"], "committed")
                    self.assertEqual(record["updated_at"], "2026-04-23T00:00:00+00:00")
                    self.assertIsNone(record["error"])
                    self.assertEqual(load_batch_manifest(manifest_path)["batch_000000"]["hit_count"], 1)

    def test_commit_prediction_batch_records_zero_hit_batch_with_empty_output(self) -> None:
        from src.screening_process_pipeline import (
            commit_prediction_batch,
            should_skip_batch,
        )

        with unittest.mock.patch("src.screening_process_pipeline._utc_now", return_value="2026-04-23T00:00:00+00:00"):
            with self._temp_directory() as tmp_dir:
                tmp_path = Path(tmp_dir)
                manifest = {}
                predicted_frame = pd.DataFrame(
                    [
                        {"chem_id": "miss-1", "smiles": "CCO", "broad_spectrum": 0},
                        {"chem_id": "miss-2", "smiles": "CCN", "broad_spectrum": 0},
                    ]
                )

                record = commit_prediction_batch(
                    output_dir=tmp_path,
                    shard_id="source_a",
                    batch_id="batch_000001",
                    predicted_frame=predicted_frame,
                    manifest_records=manifest,
                )

                output_path = tmp_path / "hits" / "source_a" / "batch_000001.parquet"
                persisted = pd.read_parquet(output_path)
                self.assertTrue(output_path.exists())
                self.assertEqual(len(persisted), 0)
                self.assertEqual(record["attempted_count"], 2)
                self.assertEqual(record["hit_count"], 0)
                self.assertEqual(record["output_path"], str(output_path))
                self.assertTrue(should_skip_batch("batch_000001", manifest))

    def test_should_skip_batch_requires_committed_status_and_existing_output(self) -> None:
        from src.screening_process_pipeline import should_skip_batch

        with self._temp_directory() as tmp_path:
            tmp_path = Path(tmp_path)
            output_path = tmp_path / "hits" / "source_a" / "batch_000000.parquet"
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_bytes(b"parquet")

            manifest = {
                "batch_000000": {
                    "status": "committed",
                    "attempted_count": 2,
                    "hit_count": 1,
                    "output_path": str(output_path),
                    "updated_at": "2026-04-23T00:00:00+00:00",
                    "error": None,
                }
            }

            self.assertTrue(should_skip_batch("batch_000000", manifest))
            output_path.unlink()
            self.assertFalse(should_skip_batch("batch_000000", manifest))
            manifest["batch_000000"]["output_path"] = ""
            self.assertFalse(should_skip_batch("batch_000000", manifest))

    def test_writer_main_materializes_manifest_and_skips_committed_batch_on_rerun(self) -> None:
        from src.screening_process_pipeline import (
            PredictedBatch,
            load_batch_manifest,
            writer_main,
            writer_run_state_path,
        )

        input_frame = pd.DataFrame([{"smiles": "CCO", "chem_id": "mol1"}])
        first_prediction = pd.DataFrame(
            [{"chem_id": "mol1", "smiles": "CCO", "apscore_total": 0.1, "broad_spectrum": 1}]
        )
        rerun_prediction = pd.DataFrame(
            [{"chem_id": "mol1", "smiles": "CCO", "apscore_total": 9.9, "broad_spectrum": 1}]
        )
        batch = PredictedBatch(
            producer_id=0,
            work_unit=WorkUnit(
                source_type="parquet",
                source_path="input.parquet",
                group_id="source_a",
                source_group="source_a",
                start_row=10,
            ),
            batch_index=0,
            input_frame=input_frame,
            predicted_frame=first_prediction,
        )
        rerun_batch = PredictedBatch(
            producer_id=0,
            work_unit=batch.work_unit,
            batch_index=0,
            input_frame=input_frame,
            predicted_frame=rerun_prediction,
        )

        with unittest.mock.patch("src.screening_process_pipeline._utc_now", return_value="2026-04-23T00:00:00+00:00"):
            with self._temp_directory() as tmp_dir:
                tmp_path = Path(tmp_dir)
                first_queue = self._queue_with_items(batch, None)
                writer_main(result_queue=first_queue, error_queue=self._queue_with_items(), output_dir=tmp_path)

                manifest_path = tmp_path / "batch_manifest.jsonl"
                manifest = load_batch_manifest(manifest_path)
                record = manifest["source_a:batch_000000"]
                output_path = Path(record["output_path"])
                first_hits = pd.read_parquet(output_path)
                self.assertEqual(record["status"], "committed")
                self.assertEqual(record["attempted_count"], 1)
                self.assertEqual(record["hit_count"], 1)
                self.assertEqual(record["error"], None)
                self.assertEqual(record["updated_at"], "2026-04-23T00:00:00+00:00")
                self.assertEqual(record["source_group"], "source_a")
                self.assertEqual(record["start_row"], 10)
                self.assertEqual(record["end_row"], 10)

                second_queue = self._queue_with_items(rerun_batch, None)
                writer_main(result_queue=second_queue, error_queue=self._queue_with_items(), output_dir=tmp_path)

                rerun_manifest = load_batch_manifest(manifest_path)
                rerun_hits = pd.read_parquet(output_path)
                run_state = json.loads(writer_run_state_path(tmp_path).read_text(encoding="utf-8"))

                self.assertEqual(rerun_manifest, manifest)
                self.assertEqual(first_hits.to_dict(orient="records"), rerun_hits.to_dict(orient="records"))
                self.assertEqual(run_state["status"], "completed")
                self.assertEqual(run_state["attempted_count"], 1)
                self.assertEqual(run_state["hit_count"], 1)
                self.assertEqual(run_state["output_path"], str(output_path))
                self.assertEqual(run_state["updated_at"], "2026-04-23T00:00:00+00:00")
                self.assertIsNone(run_state["error"])
                self.assertEqual(run_state["progress"]["committed_batches"], 1)
                self.assertEqual(run_state["progress"]["attempted_count"], 1)
                self.assertEqual(run_state["progress"]["hit_count"], 1)

    def test_writer_main_keeps_committed_manifest_when_sink_fails_after_commit(self) -> None:
        from src.screening_process_pipeline import (
            PredictedBatch,
            load_batch_manifest,
            should_skip_batch,
            writer_main,
            writer_run_state_path,
        )

        batch = PredictedBatch(
            producer_id=0,
            work_unit=WorkUnit(
                source_type="parquet",
                source_path="input.parquet",
                group_id="source_a",
                source_group="source_a",
                start_row=20,
            ),
            batch_index=1,
            input_frame=pd.DataFrame([{"smiles": "CCO", "chem_id": "mol2"}]),
            predicted_frame=pd.DataFrame(
                [{"chem_id": "mol2", "smiles": "CCO", "apscore_total": 0.2, "broad_spectrum": 1}]
            ),
        )

        with unittest.mock.patch("src.screening_process_pipeline._utc_now", return_value="2026-04-23T00:00:00+00:00"):
            with self._temp_directory() as tmp_dir:
                tmp_path = Path(tmp_dir)
                error_queue = self._queue_with_items()
                writer_main(
                    result_queue=self._queue_with_items(batch, None),
                    error_queue=error_queue,
                    output_dir=tmp_path,
                    write_batch=mock.Mock(side_effect=RuntimeError("sink failed")),
                )

                manifest = load_batch_manifest(tmp_path / "batch_manifest.jsonl")
                batch_id = "source_a:batch_000001"
                record = manifest[batch_id]
                run_state = json.loads(writer_run_state_path(tmp_path).read_text(encoding="utf-8"))
                error_record = error_queue.get_nowait()

                self.assertEqual(record["status"], "committed")
                self.assertEqual(record["hit_count"], 1)
                self.assertTrue(Path(record["output_path"]).exists())
                self.assertTrue(should_skip_batch(batch_id, manifest))
                self.assertEqual(run_state["status"], "failed")
                self.assertEqual(run_state["error"], "sink failed")
                self.assertEqual(error_record.stage, "writer")
                self.assertEqual(error_record.message, "sink failed")

    def test_resolve_graph_workers_accepts_zero(self) -> None:
        from src.screening_process_pipeline import resolve_graph_workers

        self.assertEqual(resolve_graph_workers(0, cpu_count=24), 0)

    def test_resolve_producer_processes_scales_from_cpu_count(self) -> None:
        from src.screening_process_pipeline import resolve_producer_processes

        self.assertGreaterEqual(resolve_producer_processes("auto", cpu_count=24, graph_workers=4), 6)
        self.assertGreaterEqual(resolve_producer_processes("auto", cpu_count=96, graph_workers=12), 16)

    def test_process_pipeline_uses_one_predictor_and_multiple_producers(self) -> None:
        from src.screening_process_pipeline import ProcessScreenConfig, plan_runtime

        config = ProcessScreenConfig(
            input_paths=["a.parquet", "b.parquet", "c.parquet", "d.parquet"],
            output_dir=str(Path("out")),
            execution_mode="process",
            producer_processes="auto",
            predict_queue_max_batches="auto",
            result_queue_max_batches="auto",
            batch_checkpoint_size=32,
        )

        runtime = plan_runtime(config, cpu_count=24, gpu_count=1)
        self.assertEqual(runtime.predictor_processes, 1)
        self.assertGreaterEqual(runtime.producer_processes, 6)
        self.assertEqual(runtime.work_queue_max_batches, runtime.producer_processes * 2)
        self.assertGreaterEqual(runtime.predict_queue_max_batches, 8)
        self.assertGreaterEqual(runtime.result_queue_max_batches, 8)

    def test_create_process_topology_uses_spawn_queues(self) -> None:
        from src.screening_process_pipeline import RuntimePlan, create_process_topology

        runtime = RuntimePlan(
            predictor_processes=1,
            producer_processes=6,
            work_queue_max_batches=12,
            predict_queue_max_batches=10,
            result_queue_max_batches=14,
            batch_checkpoint_size=32,
            graph_workers=2,
        )
        topology = create_process_topology(runtime)
        try:
            self.assertEqual(topology.context.get_start_method(), "spawn")
            self.assertEqual(topology.work_queue._maxsize, 12)
            self.assertEqual(topology.predict_queue._maxsize, 10)
            self.assertEqual(topology.result_queue._maxsize, 14)
        finally:
            topology.close()

    def _queue_with_items(self, *items: object) -> object:
        import queue

        q: queue.Queue[object] = queue.Queue()
        for item in items:
            q.put(item)
        return q

    def _temp_directory(self):
        return tempfile.TemporaryDirectory()

    def _run_async(self, coro):
        import asyncio

        return asyncio.run(coro)

    def test_producer_main_splits_frames_into_checkpoint_batches(self) -> None:
        from src.screening_process_pipeline import (
            ProcessWorkItem,
            RuntimePlan,
            create_process_topology,
            producer_main,
        )

        runtime = RuntimePlan(
            predictor_processes=1,
            producer_processes=1,
            work_queue_max_batches=2,
            predict_queue_max_batches=8,
            result_queue_max_batches=8,
            batch_checkpoint_size=2,
            graph_workers=0,
        )
        topology = create_process_topology(runtime)
        try:
            topology.work_queue.put(
                ProcessWorkItem(
                    work_unit=WorkUnit(
                        source_type="parquet",
                        source_path="input.parquet",
                        group_id="source_a",
                        source_group="source_a",
                    ),
                    smiles_colname="smiles",
                    chem_id_colname="chem_id",
                )
            )
            topology.work_queue.put(None)

            frame = pd.DataFrame(
                [
                    {"smiles": "CCO", "chem_id": "mol1"},
                    {"smiles": "CCN", "chem_id": "mol2"},
                    {"smiles": "CCC", "chem_id": "mol3"},
                    {"smiles": "CCF", "chem_id": "mol4"},
                    {"smiles": "CCl", "chem_id": "mol5"},
                ]
            )

            with mock.patch("src.screening_process_pipeline.process_work_unit", return_value=iter([frame])):
                producer_main(
                    producer_id=0,
                    work_queue=topology.work_queue,
                    predict_queue=topology.predict_queue,
                    error_queue=topology.error_queue,
                    batch_checkpoint_size=2,
                )

            first = topology.predict_queue.get(timeout=1)
            second = topology.predict_queue.get(timeout=1)
            third = topology.predict_queue.get(timeout=1)
            sentinel = topology.predict_queue.get(timeout=1)

            self.assertEqual([len(first.frame), len(second.frame), len(third.frame)], [2, 2, 1])
            self.assertEqual([first.batch_index, second.batch_index, third.batch_index], [0, 1, 2])
            self.assertIsNone(sentinel)
        finally:
            topology.close()

    def test_predictor_and_writer_topology_preserve_screen_frame_semantics(self) -> None:
        from src.screening_process_pipeline import (
            PredictedBatch,
            ProducerBatch,
            RuntimePlan,
            create_process_topology,
            predictor_main,
            writer_main,
        )

        runtime = RuntimePlan(
            predictor_processes=1,
            producer_processes=1,
            work_queue_max_batches=2,
            predict_queue_max_batches=8,
            result_queue_max_batches=8,
            batch_checkpoint_size=2,
            graph_workers=0,
        )
        topology = create_process_topology(runtime)
        written: list[PredictedBatch] = []
        try:
            input_frame = pd.DataFrame([{"smiles": "CCO", "chem_id": "mol1"}])
            predicted_frame = pd.DataFrame(
                [{"chem_id": "mol1", "apscore_total": 0.1, "broad_spectrum": 1}]
            )

            topology.predict_queue.put(
                ProducerBatch(
                    producer_id=0,
                    work_unit=WorkUnit(
                        source_type="parquet",
                        source_path="input.parquet",
                        group_id="source_a",
                        source_group="source_a",
                    ),
                    batch_index=0,
                    frame=input_frame,
                )
            )
            topology.predict_queue.put(None)

            with mock.patch("src.screening_process_pipeline.create_scheduler", return_value=mock.sentinel.scheduler), mock.patch(
                "src.screening_process_pipeline.screen_frame_sync",
                return_value=predicted_frame,
            ) as screen_mock:
                predictor_main(
                    predict_queue=topology.predict_queue,
                    result_queue=topology.result_queue,
                    error_queue=topology.error_queue,
                    producer_processes=1,
                    aggregate_scores=True,
                    app_threshold=0.1,
                    min_nkill=5,
                    classifier_workers=2,
                    classifier_inflight_batches=4,
                )

            result = topology.result_queue.get(timeout=1)
            sentinel = topology.result_queue.get(timeout=1)
            self.assertEqual(list(result.predicted_frame["chem_id"]), ["mol1"])
            self.assertIsNone(sentinel)

            topology.result_queue.put(result)
            topology.result_queue.put(None)
            writer_main(
                result_queue=topology.result_queue,
                error_queue=topology.error_queue,
                write_batch=written.append,
            )

            self.assertEqual(len(written), 1)
            self.assertEqual(written[0].batch_index, 0)
            screen_mock.assert_called_once()
            _, kwargs = screen_mock.call_args
            self.assertEqual(kwargs["classifier_workers"], 2)
            self.assertEqual(kwargs["classifier_inflight_batches"], 4)
        finally:
            topology.close()

    def test_predictor_main_waits_for_all_producer_sentinels(self) -> None:
        from src.screening_process_pipeline import (
            ProducerBatch,
            RuntimePlan,
            create_process_topology,
            predictor_main,
        )

        runtime = RuntimePlan(
            predictor_processes=1,
            producer_processes=2,
            work_queue_max_batches=4,
            predict_queue_max_batches=8,
            result_queue_max_batches=8,
            batch_checkpoint_size=2,
            graph_workers=0,
        )
        topology = create_process_topology(runtime)
        try:
            input_frame = pd.DataFrame([{"smiles": "CCO", "chem_id": "mol1"}])
            predicted_frame = pd.DataFrame(
                [{"chem_id": "mol1", "apscore_total": 0.1, "broad_spectrum": 1}]
            )

            topology.predict_queue.put(
                ProducerBatch(
                    producer_id=0,
                    work_unit=WorkUnit(
                        source_type="parquet",
                        source_path="input.parquet",
                        group_id="source_a",
                        source_group="source_a",
                    ),
                    batch_index=0,
                    frame=input_frame,
                )
            )
            topology.predict_queue.put(None)
            topology.predict_queue.put(None)

            with mock.patch("src.screening_process_pipeline.create_scheduler", return_value=mock.sentinel.scheduler), mock.patch(
                "src.screening_process_pipeline.screen_frame_sync",
                return_value=predicted_frame,
            ) as screen_mock:
                predictor_main(
                    predict_queue=topology.predict_queue,
                    result_queue=topology.result_queue,
                    error_queue=topology.error_queue,
                    producer_processes=2,
                    aggregate_scores=True,
                    app_threshold=0.1,
                    min_nkill=5,
                    classifier_workers=3,
                    classifier_inflight_batches=6,
                )

            result = topology.result_queue.get(timeout=1)
            sentinel = topology.result_queue.get(timeout=1)
            self.assertEqual(list(result.predicted_frame["chem_id"]), ["mol1"])
            self.assertIsNone(sentinel)
            screen_mock.assert_called_once()
            _, kwargs = screen_mock.call_args
            self.assertEqual(kwargs["classifier_workers"], 3)
            self.assertEqual(kwargs["classifier_inflight_batches"], 6)
        finally:
            topology.close()
