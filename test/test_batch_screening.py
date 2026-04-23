from __future__ import annotations

import tempfile
import unittest
import os
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pandas as pd

import mole_cli
from src.batch_screening import screen_path, ScreeningSummary

class BatchScreeningTestCase(unittest.IsolatedAsyncioTestCase):

    async def test_screen_path_streams_chunks(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.tsv"
            output_dir = Path(tmpdir) / "out"
            input_path.write_text("smiles\nCCO\nCCN\nCCC\nCCF\n", encoding="utf-8")

            # Predictor mockup
            async def fake_predict_molecules(molecules, **kwargs):
                return [{"chem_id": m.chem_id, "score": 0.5, "pred_id": f"{m.chem_id}:strain1", "source_group": "input"} for m in molecules]

            mock_scheduler = mock.MagicMock()
            mock_scheduler.predict_molecules = mock.AsyncMock(side_effect=fake_predict_molecules)

            from src.screening_planner import WorkUnit
            mock_plan = mock.MagicMock(return_value=[
                WorkUnit(source_type="tabular", source_path=str(input_path), group_id="input_part1", start_row=0, max_rows=2),
                WorkUnit(source_type="tabular", source_path=str(input_path), group_id="input_part2", start_row=2, max_rows=2),
            ])

            with mock.patch("src.batch_screening.get_scheduler", return_value=mock_scheduler), \
                 mock.patch("src.batch_screening.ScreeningPlanner.plan", mock_plan):
                summary = await screen_path(
                    input_path=input_path,
                    output_dir=output_dir,
                    chunk_size=2,
                    target_rows_per_group=2,
                    grouping_mode="chunk",
                    aggregate_scores=True
                )

            # Since chunk_size=2 and we have 4 inputs, predict_molecules should be called twice
            self.assertEqual(mock_scheduler.predict_molecules.call_count, 2)
            self.assertEqual(summary.normalized_rows, 4)
            self.assertEqual(summary.predicted_rows, 4)
            self.assertTrue((output_dir / "normalized_input.tsv").exists())
            self.assertTrue((output_dir / "predictions_all.tsv").exists())
            self.assertTrue((output_dir / "by_source" / "input" / "predictions.tsv").exists())

    async def test_screen_path_prefetches_chunks_concurrently(self) -> None:
        import asyncio
        import time
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.tsv"
            output_dir = Path(tmpdir) / "out"
            # 5 items, chunk size 1 = 5 chunks
            input_path.write_text("smiles\nC\nCC\nCCC\nCCCC\nCCCCC\n", encoding="utf-8")

            produced_times = []
            consumed_times = []

            original_process = mole_cli.screen_path.__globals__.get("process_work_unit")
            if not original_process:
                from src.screening_sources import process_work_unit as original_process

            def slow_process(*args, **kwargs):
                produced_times.append(time.time())
                return original_process(*args, **kwargs)

            async def slow_predict_molecules(molecules, **kwargs):
                await asyncio.sleep(0.05)
                consumed_times.append(time.time())
                return [{"chem_id": m.chem_id, "score": 0.5, "pred_id": f"{m.chem_id}:strain1"} for m in molecules]

            mock_scheduler = mock.MagicMock()
            mock_scheduler.predict_molecules = mock.AsyncMock(side_effect=slow_predict_molecules)

            from src.screening_planner import WorkUnit
            mock_plan = mock.MagicMock(return_value=[
                WorkUnit(source_type="tabular", source_path=str(input_path), group_id=f"input_part{i+1}", start_row=i, max_rows=1)
                for i in range(5)
            ])

            with mock.patch("src.batch_screening.get_scheduler", return_value=mock_scheduler), \
                 mock.patch("src.batch_screening.process_work_unit", side_effect=slow_process), \
                 mock.patch("src.batch_screening.ScreeningPlanner.plan", mock_plan):
                # prefetch 3 chunks
                summary = await screen_path(
                    input_path=input_path,
                    output_dir=output_dir,
                    chunk_size=1,
                    target_rows_per_group=1,
                    grouping_mode="chunk",
                    cpu_workers=5,
                    prefetch_queue_size=3,
                    aggregate_scores=True
                )

            t0 = produced_times[0]
            first_consumer_end = consumed_times[0]

            produced_before_first_consume = [t for t in produced_times if t < first_consumer_end]
            # Must be at least 3 chunks prefetched concurrently
            self.assertGreaterEqual(len(produced_before_first_consume), 3)

    async def test_screen_path_incremental_writing(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.tsv"
            output_dir = Path(tmpdir) / "out"
            input_path.write_text("smiles\nCCO\nCCN\nCCC\nCCF\n", encoding="utf-8")

            write_counts = [0]
            async def fake_predict_molecules(molecules, **kwargs):
                # check if file written yet
                pred_path = output_dir / "predictions_all.tsv"
                if write_counts[0] == 0:
                    self.assertFalse(pred_path.exists())
                else:
                    self.assertTrue(pred_path.exists())
                write_counts[0] += 1
                return [{"chem_id": m.chem_id, "score": 0.5, "pred_id": f"{m.chem_id}:strain1"} for m in molecules]

            mock_scheduler = mock.MagicMock()
            mock_scheduler.predict_molecules = mock.AsyncMock(side_effect=fake_predict_molecules)

            with mock.patch("src.batch_screening.get_scheduler", return_value=mock_scheduler):
                summary = await screen_path(
                    input_path=input_path,
                    output_dir=output_dir,
                    chunk_size=2,
                    aggregate_scores=True
                )

            # verify 1 group created
            self.assertEqual(len(summary.grouped_outputs), 1)
            self.assertEqual(summary.grouped_outputs[0]["source_group"], "input")

            # verify predictions_all
            df = pd.read_csv(output_dir / "predictions_all.tsv", sep="\t")
            self.assertEqual(len(df), 4)

    async def test_screen_path_stops_producer_on_consumer_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.tsv"
            output_dir = Path(tmpdir) / "out"
            # Provide many items so the producer normally has to work for a while
            smiles_list = ["C" * (i+1) for i in range(1, 101)]
            content = "smiles\tchem_id\n" + "\n".join([f"{s}\t{s}" for s in smiles_list])
            input_path.write_text(content, encoding="utf-8")

            # Use a mock stream that tracks how many items are yielded
            original_process_fn = mole_cli.screen_path.__globals__.get("process_work_unit")
            if not original_process_fn:
                from src.screening_sources import process_work_unit as original_process_fn

            yield_count = 0
            async def fake_predict_molecules(molecules, **kwargs):
                nonlocal yield_count
                yield_count += 1
                if yield_count == 2:
                    # Throw exception on the 2nd chunk
                    raise RuntimeError("Simulated consumer failure")
                return [{"chem_id": m.chem_id, "score": 0.5, "pred_id": f"{m.chem_id}:strain1"} for m in molecules]

            mock_scheduler = mock.MagicMock()
            mock_scheduler.predict_molecules = mock.AsyncMock(side_effect=fake_predict_molecules)

            producer_yielded = 0
            producer_finished = False

            def tracking_process(*args, **kwargs):
                nonlocal producer_yielded, producer_finished
                try:
                    res = original_process_fn(*args, **kwargs)
                    for chunk in res:
                        producer_yielded += 1
                        yield chunk
                finally:
                    producer_finished = True

            from src.screening_planner import WorkUnit
            mock_plan = mock.MagicMock(return_value=[
                WorkUnit(source_type="tabular", source_path=str(input_path), group_id=f"input_part{i+1}", start_row=i, max_rows=1)
                for i in range(100)
            ])
            with mock.patch("src.batch_screening.get_scheduler", return_value=mock_scheduler), \
                 mock.patch("src.batch_screening.process_work_unit", side_effect=tracking_process), \
                 mock.patch("src.batch_screening.ScreeningPlanner.plan", mock_plan):

                with self.assertRaises(RuntimeError) as context:
                    await screen_path(
                        input_path=input_path,
                        output_dir=output_dir,
                        chunk_size=1, # 1 item per chunk to easily control iterations
                        target_rows_per_group=1,
                        grouping_mode="chunk",
                        cpu_workers=2,
                        prefetch_queue_size=2,
                        aggregate_scores=True
                    )

                self.assertIn("Simulated consumer failure", str(context.exception))

                # Check that producer task is actually done and cleanup happened
                self.assertTrue(producer_finished)

                # Producer might have prefetched slightly more than the 2 consumed chunks
                # but it should definitely not have processed all 100 chunks.
                self.assertLess(producer_yielded, 100)

    async def test_screen_path_worker_processes_chunks_and_preserves_work_unit_order(self) -> None:
        import asyncio
        import time
        import threading
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.tsv"
            output_dir = Path(tmpdir) / "out"
            input_path.write_text("smiles\nCCO\n", encoding="utf-8")

            main_thread_id = threading.get_ident()
            yield_threads = set()
            yield_sequence = []

            def mock_process_work_unit(unit, *args, **kwargs):
                if unit.group_id == "large":
                    yield pd.DataFrame([{"smiles": "CCO", "chem_id": "A1"}])
                    yield_threads.add(threading.get_ident())
                    time.sleep(0.2)
                    yield pd.DataFrame([{"smiles": "CCN", "chem_id": "A2"}])
                    yield_threads.add(threading.get_ident())
                elif unit.group_id == "small":
                    time.sleep(0.05)
                    yield pd.DataFrame([{"smiles": "CCC", "chem_id": "B1"}])
                    yield_threads.add(threading.get_ident())

            async def mock_predict_molecules(molecules, **kwargs):
                for m in molecules:
                    yield_sequence.append(m.chem_id)
                return [{"chem_id": m.chem_id, "score": 0.5, "pred_id": f"{m.chem_id}:strain1", "source_group": "input"} for m in molecules]

            mock_scheduler = mock.MagicMock()
            mock_scheduler.predict_molecules = mock.AsyncMock(side_effect=mock_predict_molecules)

            from src.screening_planner import WorkUnit
            mock_plan = mock.MagicMock(return_value=[
                WorkUnit(source_type="tabular", source_path="dummy", group_id="large", start_row=0, max_rows=2),
                WorkUnit(source_type="tabular", source_path="dummy", group_id="small", start_row=0, max_rows=1),
            ])

            with mock.patch("src.batch_screening.get_scheduler", return_value=mock_scheduler), \
                 mock.patch("src.batch_screening.process_work_unit", side_effect=mock_process_work_unit), \
                 mock.patch("src.batch_screening.ScreeningPlanner.plan", mock_plan):

                await screen_path(
                    input_path=input_path,
                    output_dir=output_dir,
                    chunk_size=1,
                    cpu_workers=2,
                    prefetch_queue_size=4,
                    aggregate_scores=True
                )

            self.assertNotIn(main_thread_id, yield_threads)
            self.assertTrue(len(yield_threads) > 0)
            self.assertEqual(yield_sequence, ["A1", "A2", "B1"])

    async def test_screen_path_dedupe_keeps_first_source_occurrence(self) -> None:
        import time
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.tsv"
            output_dir = Path(tmpdir) / "out"
            input_path.write_text("smiles\nCCO\n", encoding="utf-8")

            def mock_process_work_unit(unit, *args, **kwargs):
                if unit.group_id == "first":
                    time.sleep(0.15)
                    yield pd.DataFrame([{"smiles": "CCO", "chem_id": "first_hit"}])
                elif unit.group_id == "second":
                    yield pd.DataFrame([{"smiles": "CCO", "chem_id": "second_hit"}])

            async def mock_predict_molecules(molecules, **kwargs):
                return [
                    {
                        "chem_id": molecules[0].chem_id,
                        "score": 0.5,
                        "pred_id": f"{molecules[0].chem_id}:strain1",
                        "source_group": "input",
                    }
                ]

            mock_scheduler = mock.MagicMock()
            mock_scheduler.predict_molecules = mock.AsyncMock(side_effect=mock_predict_molecules)

            from src.screening_planner import WorkUnit
            mock_plan = mock.MagicMock(return_value=[
                WorkUnit(source_type="tabular", source_path="dummy", group_id="first", start_row=0, max_rows=1),
                WorkUnit(source_type="tabular", source_path="dummy", group_id="second", start_row=1, max_rows=1),
            ])

            with mock.patch("src.batch_screening.get_scheduler", return_value=mock_scheduler), \
                 mock.patch("src.batch_screening.process_work_unit", side_effect=mock_process_work_unit), \
                 mock.patch("src.batch_screening.ScreeningPlanner.plan", mock_plan):
                await screen_path(
                    input_path=input_path,
                    output_dir=output_dir,
                    chunk_size=1,
                    cpu_workers=2,
                    prefetch_queue_size=4,
                    dedupe_smiles=True,
                    aggregate_scores=True,
                )

            predicted = pd.read_csv(output_dir / "predictions_all.tsv", sep="\t")
            self.assertEqual(predicted["chem_id"].tolist(), ["first_hit"])

    async def test_screen_path_globally_sorts_aggregate_predictions_across_chunks(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.tsv"
            output_dir = Path(tmpdir) / "out"
            input_path.write_text("smiles\tchem_id\nCCO\tmol1\nCCN\tmol2\n", encoding="utf-8")

            async def fake_predict_molecules(molecules, **kwargs):
                chem_id = molecules[0].chem_id
                if chem_id == "mol1":
                    return [{
                        "chem_id": "mol1",
                        "broad_spectrum": 0,
                        "ginhib_total": 1,
                        "apscore_total": 0.9,
                    }]
                return [{
                    "chem_id": "mol2",
                    "broad_spectrum": 1,
                    "ginhib_total": 10,
                    "apscore_total": 0.1,
                }]

            mock_scheduler = mock.MagicMock()
            mock_scheduler.predict_molecules = mock.AsyncMock(side_effect=fake_predict_molecules)

            from src.screening_planner import WorkUnit
            mock_plan = mock.MagicMock(return_value=[
                WorkUnit(source_type="tabular", source_path=str(input_path), group_id="input_part1", source_group="input", start_row=0, max_rows=1),
                WorkUnit(source_type="tabular", source_path=str(input_path), group_id="input_part2", source_group="input", start_row=1, max_rows=1),
            ])

            with mock.patch("src.batch_screening.get_scheduler", return_value=mock_scheduler), \
                 mock.patch("src.batch_screening.ScreeningPlanner.plan", mock_plan):
                await screen_path(
                    input_path=input_path,
                    output_dir=output_dir,
                    chunk_size=1,
                    aggregate_scores=True,
                )

            predicted = pd.read_csv(output_dir / "predictions_all.tsv", sep="\t")
            by_source_predicted = pd.read_csv(
                output_dir / "by_source" / "input" / "predictions.tsv",
                sep="\t",
            )
            self.assertEqual(predicted["chem_id"].tolist(), ["mol2", "mol1"])
            self.assertEqual(by_source_predicted["chem_id"].tolist(), ["mol2", "mol1"])

    async def test_screen_path_passes_chunk_size_to_worker(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.tsv"
            output_dir = Path(tmpdir) / "out"
            input_path.write_text("smiles\nCCO\n", encoding="utf-8")

            # Predictor mockup
            async def fake_predict_molecules(molecules, **kwargs):
                return [{"chem_id": m.chem_id, "score": 0.5, "pred_id": f"{m.chem_id}:strain1", "source_group": "input"} for m in molecules]

            mock_scheduler = mock.MagicMock()
            mock_scheduler.predict_molecules = mock.AsyncMock(side_effect=fake_predict_molecules)

            from src.screening_planner import WorkUnit
            mock_plan = mock.MagicMock(return_value=[
                WorkUnit(source_type="tabular", source_path=str(input_path), group_id="input_part1", start_row=0, max_rows=1),
            ])

            import src.batch_screening
            original_process = src.batch_screening.process_work_unit
            mock_process = mock.Mock(side_effect=original_process)

            with mock.patch("src.batch_screening.get_scheduler", return_value=mock_scheduler), \
                 mock.patch("src.batch_screening.process_work_unit", mock_process), \
                 mock.patch("src.batch_screening.ScreeningPlanner.plan", mock_plan):

                custom_chunk_size = 42

                await screen_path(
                    input_path=input_path,
                    output_dir=output_dir,
                    chunk_size=custom_chunk_size,
                    aggregate_scores=True
                )

            mock_process.assert_called()
            # Check if mock_process was called with chunk_size=42
            found = False
            for call in mock_process.call_args_list:
                args, kwargs = call
                if kwargs.get("chunk_size") == custom_chunk_size:
                    found = True
                    break
            self.assertTrue(found, f"process_work_unit not called with chunk_size={custom_chunk_size}")

    async def test_screen_path_explicitly_forwards_graph_build_settings(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.tsv"
            output_dir = Path(tmpdir) / "out"
            input_path.write_text("smiles\tchem_id\nCCO\tmol1\n", encoding="utf-8")

            async def fake_predict_molecules(molecules, **kwargs):
                return [
                    {
                        "chem_id": molecules[0].chem_id,
                        "score": 0.5,
                        "pred_id": f"{molecules[0].chem_id}:strain1",
                        "source_group": "input",
                    }
                ]

            mock_scheduler = mock.MagicMock()
            mock_scheduler.predict_molecules = mock.AsyncMock(side_effect=fake_predict_molecules)

            with mock.patch("src.batch_screening.get_scheduler", return_value=mock_scheduler):
                await screen_path(
                    input_path=input_path,
                    output_dir=output_dir,
                    chunk_size=1,
                    aggregate_scores=True,
                    num_graph_workers=3,
                    graph_batch_size=128,
                    prefetch_batches=4,
                )

            _, kwargs = mock_scheduler.predict_molecules.call_args
            self.assertEqual(kwargs["num_graph_workers"], 3)
            self.assertEqual(kwargs["graph_batch_size"], 128)
            self.assertEqual(kwargs["prefetch_batches"], 4)

    async def test_screen_path_uses_explicit_scheduler_instead_of_singleton(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.tsv"
            output_dir = Path(tmpdir) / "out"
            input_path.write_text("smiles\tchem_id\nCCO\tmol1\n", encoding="utf-8")

            explicit_scheduler = mock.MagicMock()
            explicit_scheduler.predict_molecules = mock.AsyncMock(
                return_value=[
                    {
                        "chem_id": "mol1",
                        "score": 0.5,
                        "pred_id": "mol1:strain1",
                        "source_group": "input",
                    }
                ]
            )

            with mock.patch("src.batch_screening.get_scheduler", side_effect=AssertionError("singleton should not be used")):
                await screen_path(
                    input_path=input_path,
                    output_dir=output_dir,
                    chunk_size=1,
                    aggregate_scores=True,
                    scheduler=explicit_scheduler,
                )

            explicit_scheduler.predict_molecules.assert_awaited_once()


    def test_parser_exposes_screen_subcommand(self) -> None:
        parser = mole_cli._build_parser()
        args = parser.parse_args(["screen", "--input-path", "x.tsv", "--output-dir", "out"])
        self.assertEqual(args.command, "screen")

    def test_screen_command_invokes_batch_screening_writer(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.tsv"
            output_dir = Path(tmpdir) / "out"
            input_path.write_text("smiles\nCCO\n", encoding="utf-8")
            args = SimpleNamespace(
                input_path=str(input_path),
                output_dir=str(output_dir),
                smiles_colname="smiles",
                chem_id_colname="chem_id",
                archive_pattern="*_scheme_b_unique_products.csv",
                archive_smiles_colname="product_smiles_canonical",
                archive_chem_id_colname="example_combo_id",
                dedupe_smiles=True,
                aggregate_scores=True,
                app_threshold=0.04374140128493309,
                min_nkill=10,
                input_chunk_size=256,
                max_batch_size=2048,
                target_gpu_memory_fraction=0.8,
                num_graph_workers=3,
                graph_batch_size=128,
                prefetch_batches=4,
                output_format="tsv",
            )

            summary = ScreeningSummary(
                normalized_rows=1,
                predicted_rows=1,
                normalized_input_path=output_dir / "normalized_input.tsv",
                predictions_all_path=output_dir / "predictions_all.tsv",
                grouped_outputs=[],
                grouping_mode="auto",
                cpu_workers_selected=1,
                prefetch_queue_size_selected=4,
                work_unit_count=1,
                target_rows_per_group=100000,
                target_bytes_per_group=50000000
            )

            mock_scheduler = mock.MagicMock()
            mock_scheduler.runtime_snapshot.return_value = {"used_cuda": False, "selected_batch_size": 256}

            with mock.patch.object(mole_cli, "screen_path", return_value=summary) as screen_path_mock, \
                 mock.patch.object(mole_cli, "create_scheduler", return_value=mock_scheduler):
                exit_code = mole_cli._command_screen(args)

            self.assertEqual(exit_code, 0)
            self.assertTrue(output_dir.exists())
            self.assertTrue((output_dir / "manifest.json").is_file())
            _, kwargs = screen_path_mock.call_args
            self.assertIs(kwargs["scheduler"], mock_scheduler)
            self.assertEqual(kwargs["num_graph_workers"], 3)
            self.assertEqual(kwargs["graph_batch_size"], 128)
            self.assertEqual(kwargs["prefetch_batches"], 4)

    async def test_screen_path_accepts_parquet_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            shard = Path(tmpdir) / "shard_0001.parquet"
            output_dir = Path(tmpdir) / "out"
            pd.DataFrame(
                {
                    "smiles": ["CCO", "CCN"],
                    "chem_id": ["m1", "m2"],
                    "source_group": ["demo", "demo"],
                }
            ).to_parquet(shard, index=False, engine="pyarrow", row_group_size=1)

            async def fake_predict_molecules(molecules, **kwargs):
                return [
                    {"chem_id": m.chem_id, "score": 0.5, "pred_id": f"{m.chem_id}:strain1"}
                    for m in molecules
                ]

            mock_scheduler = mock.MagicMock()
            mock_scheduler.predict_molecules = mock.AsyncMock(side_effect=fake_predict_molecules)

            with mock.patch("src.batch_screening.get_scheduler", return_value=mock_scheduler):
                summary = await screen_path(
                    input_path=shard,
                    output_dir=output_dir,
                    chunk_size=1,
                    aggregate_scores=True,
                )

            self.assertEqual(summary.normalized_rows, 2)
            self.assertEqual(summary.predicted_rows, 2)
            self.assertTrue((output_dir / "predictions_all.tsv").exists())
            self.assertTrue((output_dir / "by_source" / "demo" / "predictions.tsv").exists())

    async def test_screen_path_collects_profiling_summary_when_enabled(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.tsv"
            output_dir = Path(tmpdir) / "out"
            input_path.write_text("smiles\nCCO\nCCN\n", encoding="utf-8")

            async def fake_predict_molecules(molecules, **kwargs):
                return [{"chem_id": m.chem_id, "score": 0.5, "pred_id": f"{m.chem_id}:strain1"} for m in molecules]

            mock_scheduler = mock.MagicMock()
            mock_scheduler.predict_molecules = mock.AsyncMock(side_effect=fake_predict_molecules)
            mock_scheduler.runtime_snapshot = mock.Mock(
                return_value={
                    "device": "cpu",
                    "selected_batch_size": 32,
                    "last_profile": {
                        "representation_seconds": 1.0,
                        "prediction_frame_seconds": 0.25,
                        "growth_inhibition_seconds": 0.1,
                        "aggregate_scores_seconds": 0.2,
                        "xgboost_seconds": 0.5,
                        "graph_build": {
                            "graph_total_seconds": 0.8,
                            "graph_items": 2,
                        },
                    },
                }
            )

            summary = await screen_path(
                input_path=input_path,
                output_dir=output_dir,
                aggregate_scores=True,
                scheduler=mock_scheduler,
                enable_profiling=True,
            )

            self.assertIsInstance(summary.profiling, dict)
            self.assertIn("queue_wait_seconds", summary.profiling)
            self.assertIn("screen_frame_seconds", summary.profiling)
            self.assertEqual(summary.profiling["chunks_predicted"], 1)
            self.assertAlmostEqual(summary.profiling["prediction_seconds"], 2.05)
            self.assertEqual(summary.profiling["graph_build_seconds"], 0.8)
            self.assertAlmostEqual(summary.profiling["aggregate_scores_seconds"], 0.2)

    async def test_screen_path_coalesces_ready_frames_up_to_prediction_row_budget(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.tsv"
            output_dir = Path(tmpdir) / "out"
            input_path.write_text("smiles\nC\nCC\nCCC\nCCCC\n", encoding="utf-8")

            call_sizes = []

            async def fake_predict_molecules(molecules, **kwargs):
                call_sizes.append(len(molecules))
                return [{"chem_id": m.chem_id, "score": 0.5, "pred_id": f"{m.chem_id}:strain1"} for m in molecules]

            mock_scheduler = mock.MagicMock()
            mock_scheduler.predict_molecules = mock.AsyncMock(side_effect=fake_predict_molecules)

            from src.screening_planner import WorkUnit
            mock_plan = mock.MagicMock(return_value=[
                WorkUnit(source_type="tabular", source_path=str(input_path), group_id=f"input_part{i+1}", start_row=i, max_rows=1)
                for i in range(4)
            ])

            with mock.patch("src.batch_screening.get_scheduler", return_value=mock_scheduler), \
                 mock.patch("src.batch_screening.ScreeningPlanner.plan", mock_plan):
                summary = await screen_path(
                    input_path=input_path,
                    output_dir=output_dir,
                    chunk_size=1,
                    target_rows_per_group=1,
                    grouping_mode="chunk",
                    aggregate_scores=True,
                    prediction_row_budget=2,
                )

            self.assertEqual(summary.normalized_rows, 4)
            self.assertEqual(summary.predicted_rows, 4)
            self.assertLess(len(call_sizes), 4)
            self.assertIn(2, call_sizes)

    async def test_screen_path_drains_queue_while_prediction_is_in_flight(self) -> None:
        import asyncio

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.tsv"
            output_dir = Path(tmpdir) / "out"
            input_path.write_text("smiles\nC\nCC\nCCC\nCCCC\nCCCCC\n", encoding="utf-8")

            producer_yielded = 0
            prediction_started = asyncio.Event()
            release_prediction = asyncio.Event()

            def tracking_process(unit, *args, **kwargs):
                nonlocal producer_yielded
                yield pd.DataFrame(
                    [
                        {
                            "smiles": "C" * (unit.start_row + 1),
                            "chem_id": f"mol{unit.start_row + 1}",
                        }
                    ]
                )
                producer_yielded += 1

            async def blocked_predict_molecules(molecules, **kwargs):
                prediction_started.set()
                await release_prediction.wait()
                return [
                    {"chem_id": m.chem_id, "score": 0.5, "pred_id": f"{m.chem_id}:strain1"}
                    for m in molecules
                ]

            mock_scheduler = mock.MagicMock()
            mock_scheduler.predict_molecules = mock.AsyncMock(side_effect=blocked_predict_molecules)

            from src.screening_planner import WorkUnit
            mock_plan = mock.MagicMock(
                return_value=[
                    WorkUnit(
                        source_type="tabular",
                        source_path=str(input_path),
                        group_id=f"input_part{i+1}",
                        start_row=i,
                        max_rows=1,
                    )
                    for i in range(5)
                ]
            )

            with mock.patch("src.batch_screening.get_scheduler", return_value=mock_scheduler), \
                 mock.patch("src.batch_screening.process_work_unit", side_effect=tracking_process), \
                 mock.patch("src.batch_screening.ScreeningPlanner.plan", mock_plan):
                run_task = asyncio.create_task(
                    screen_path(
                        input_path=input_path,
                        output_dir=output_dir,
                        chunk_size=1,
                        target_rows_per_group=1,
                        grouping_mode="chunk",
                        cpu_workers=5,
                        prefetch_queue_size=1,
                        aggregate_scores=True,
                    )
                )

                await asyncio.wait_for(prediction_started.wait(), timeout=1.0)
                await asyncio.sleep(0.1)
                self.assertGreaterEqual(
                    producer_yielded,
                    4,
                    "consumer should keep draining the queue while prediction is running so producer threads do not stall at queue capacity",
                )

                release_prediction.set()
                summary = await asyncio.wait_for(run_task, timeout=2.0)

            self.assertEqual(summary.normalized_rows, 5)
            self.assertEqual(summary.predicted_rows, 5)


if __name__ == "__main__":
    unittest.main()
