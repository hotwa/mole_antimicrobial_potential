from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import mole_cli


class MoleCliTestCase(unittest.TestCase):
    def test_parser_exposes_expected_subcommands(self):
        parser = mole_cli._build_parser()
        self.assertEqual(parser.parse_args(["doctor"]).command, "doctor")
        self.assertEqual(
            parser.parse_args(["predict", "--smiles", "CCO"]).command,
            "predict",
        )
        self.assertEqual(
            parser.parse_args(["score", "--objective-file", "x.json", "--smiles", "CCO"]).command,
            "score",
        )
        self.assertEqual(
            parser.parse_args(["benchmark-screening-inputs", "--input-path", "input.csv"]).command,
            "benchmark-screening-inputs",
        )

    def test_parser_accepts_graph_construction_tuning_flags(self):
        parser = mole_cli._build_parser()
        args = parser.parse_args(
            [
                "predict",
                "--smiles",
                "CCO",
                "--num-graph-workers",
                "3",
                "--graph-batch-size",
                "128",
                "--prefetch-batches",
                "4",
                "--classifier-workers",
                "2",
                "--classifier-inflight-batches",
                "5",
                "--deterministic-representation",
                "--classifier-backend",
                "pickle",
            ]
        )
        self.assertEqual(args.num_graph_workers, "3")
        self.assertEqual(args.graph_batch_size, 128)
        self.assertEqual(args.prefetch_batches, 4)
        self.assertEqual(args.classifier_workers, "2")
        self.assertEqual(args.classifier_inflight_batches, "5")
        self.assertTrue(args.deterministic_representation)
        self.assertEqual(args.classifier_backend, "pickle")

    def test_parser_accepts_profiling_flags(self):
        parser = mole_cli._build_parser()
        predict_args = parser.parse_args(["predict", "--smiles", "CCO", "--profiling"])
        screen_args = parser.parse_args(["screen", "--input-path", "x.tsv", "--output-dir", "out", "--profiling"])
        self.assertTrue(predict_args.profiling)
        self.assertTrue(screen_args.profiling)

    def test_parser_accepts_prediction_row_budget(self):
        parser = mole_cli._build_parser()
        args = parser.parse_args(
            ["screen", "--input-path", "x.tsv", "--output-dir", "out", "--prediction-row-budget", "8192"]
        )
        self.assertEqual(args.prediction_row_budget, 8192)

    def test_screen_parser_accepts_repeated_input_paths(self):
        parser = mole_cli._build_parser()
        args = parser.parse_args(
            [
                "screen",
                "--input-path",
                "a.csv",
                "--input-path",
                "b.csv",
                "--output-dir",
                "out",
                "--execution-mode",
                "process",
            ]
        )

        self.assertEqual(args.input_path, ["a.csv", "b.csv"])
        self.assertEqual(args.execution_mode, "process")

    def test_screen_parser_accepts_process_tuning_flags(self):
        parser = mole_cli._build_parser()
        args = parser.parse_args(
            [
                "screen",
                "--input-path",
                "a.csv",
                "--output-dir",
                "out",
                "--execution-mode",
                "process",
                "--producer-processes",
                "auto",
                "--predict-queue-max-batches",
                "auto",
                "--result-queue-max-batches",
                "8",
                "--batch-checkpoint-size",
                "2048",
            ]
        )

        self.assertEqual(args.producer_processes, "auto")
        self.assertEqual(args.predict_queue_max_batches, "auto")
        self.assertEqual(args.result_queue_max_batches, "8")
        self.assertEqual(args.batch_checkpoint_size, 2048)

    def test_screen_process_config_tracks_multi_input_and_process_flags(self):
        from src.screening_process_pipeline import process_screen_config_from_args

        parser = mole_cli._build_parser()
        args = parser.parse_args(
            [
                "screen",
                "--input-path",
                "a.csv",
                "--input-path",
                "b.csv",
                "--output-dir",
                "out",
                "--execution-mode",
                "process",
                "--producer-processes",
                "3",
                "--predict-queue-max-batches",
                "auto",
                "--result-queue-max-batches",
                "8",
                "--batch-checkpoint-size",
                "4096",
            ]
        )

        config = process_screen_config_from_args(args)

        self.assertEqual(config.input_paths, ["a.csv", "b.csv"])
        self.assertEqual(config.output_dir, "out")
        self.assertEqual(config.execution_mode, "process")
        self.assertEqual(config.producer_processes, "3")
        self.assertEqual(config.predict_queue_max_batches, "auto")
        self.assertEqual(config.result_queue_max_batches, "8")
        self.assertEqual(config.batch_checkpoint_size, 4096)

    def test_screen_process_config_wraps_scalar_string_input_path(self):
        from src.screening_process_pipeline import process_screen_config_from_args

        args = SimpleNamespace(
            input_path="/tmp/input.csv",
            output_dir="out",
            execution_mode="process",
        )

        config = process_screen_config_from_args(args)

        self.assertEqual(config.input_paths, ["/tmp/input.csv"])

    def test_command_screen_rejects_multi_input_in_thread_mode(self):
        parser = mole_cli._build_parser()
        args = parser.parse_args(
            [
                "screen",
                "--input-path",
                "a.csv",
                "--input-path",
                "b.csv",
                "--output-dir",
                "out",
            ]
        )

        with self.assertRaisesRegex(
            ValueError,
            "screen thread mode currently supports exactly one --input-path; use --execution-mode process for multiple inputs",
        ):
            mole_cli._command_screen(args)

    def test_command_screen_thread_mode_does_not_forward_process_config_kw(self):
        parser = mole_cli._build_parser()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_path = tmp / "input.tsv"
            output_dir = tmp / "out"
            input_path.write_text("smiles\tchem_id\nCCO\tmol1\n", encoding="utf-8")
            args = parser.parse_args(
                [
                    "screen",
                    "--input-path",
                    str(input_path),
                    "--output-dir",
                    str(output_dir),
                    "--classifier-workers",
                    "2",
                    "--classifier-inflight-batches",
                    "5",
                ]
            )
            summary = SimpleNamespace(
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
                target_bytes_per_group=50000000,
                profiling=None,
            )
            fake_scheduler = SimpleNamespace(runtime_snapshot=mock.Mock(return_value={"used_cuda": False}))

            async def fake_screen_path(**kwargs):
                self.assertNotIn("process_config", kwargs)
                self.assertEqual(kwargs["input_path"], str(input_path))
                self.assertEqual(kwargs["classifier_workers"], "2")
                self.assertEqual(kwargs["classifier_inflight_batches"], "5")
                return summary

            with mock.patch.object(mole_cli, "create_scheduler", return_value=fake_scheduler) as create_scheduler_mock, mock.patch.object(
                mole_cli, "screen_path", side_effect=fake_screen_path
            ), mock.patch.object(mole_cli, "_dump_json"):
                self.assertEqual(mole_cli._command_screen(args), 0)
            _, scheduler_kwargs = create_scheduler_mock.call_args
            self.assertEqual(scheduler_kwargs["classifier_workers"], "2")
            self.assertEqual(scheduler_kwargs["classifier_inflight_batches"], "5")

    def test_screen_command_routes_process_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_path = tmp / "a.csv"
            output_dir = tmp / "out"
            input_path.write_text("smiles,chem_id\nCCO,m1\n", encoding="utf-8")

            process_summary = {
                "execution_mode": "process",
                "prepared_manifest_path": str(output_dir / "prepared_manifest.json"),
                "prepared_manifest_paths": [],
                "batch_manifest_path": str(output_dir / "batch_manifest.jsonl"),
                "run_state_path": str(output_dir / "run_state.json"),
                "hits_dir": str(output_dir / "hits"),
                "input_paths": [str(input_path.resolve())],
                "prepared_input_paths": [],
                "source_groups": ["a"],
                "work_unit_count": 1,
                "runtime": {
                    "predictor_processes": 1,
                    "producer_processes": 2,
                    "predict_queue_max_batches": 8,
                    "result_queue_max_batches": 8,
                    "batch_checkpoint_size": 2048,
                    "graph_workers": 0,
                },
                "prediction_runtime": None,
            }

            async def fake_screen_paths_multiprocess(**kwargs):
                self.assertEqual(kwargs["input_paths"], [str(input_path)])
                self.assertEqual(kwargs["output_dir"], str(output_dir))
                self.assertEqual(kwargs["execution_mode"], "process")
                self.assertEqual(kwargs["classifier_workers"], "3")
                self.assertEqual(kwargs["classifier_inflight_batches"], "6")
                return process_summary

            with mock.patch.object(
                mole_cli,
                "screen_paths_multiprocess",
                side_effect=fake_screen_paths_multiprocess,
                create=True,
            ) as process_mock, mock.patch.object(
                mole_cli,
                "create_scheduler",
                side_effect=AssertionError("thread scheduler should not be created in process mode"),
            ), mock.patch.object(mole_cli, "_dump_json"):
                rc = mole_cli.main(
                    [
                        "screen",
                        "--input-path",
                        str(input_path),
                        "--output-dir",
                        str(output_dir),
                        "--execution-mode",
                        "process",
                        "--classifier-workers",
                        "3",
                        "--classifier-inflight-batches",
                        "6",
                    ]
                )

            self.assertEqual(rc, 0)
            process_mock.assert_called_once()

    def test_doctor_reports_ok_when_fake_probe_is_available(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            mole_model = tmp / "model"
            mole_model.mkdir()
            (mole_model / "config.yaml").write_text("model: {}", encoding="utf-8")
            (mole_model / "model.pth").write_bytes(b"stub")
            args = SimpleNamespace(
                strict_gpu=False,
                env_file=str(Path("workflows/reinvent4/configs/local.env.recommended.example")),
                scaffold_file=str(Path("workflows/reinvent4/inputs/scaffolds/mother_scaffold.template.smi")),
                objective_file=str(Path("workflows/reinvent4/inputs/objectives/pathogen_group_a.site_reward.prototype.json")),
            )
            fake_probe = SimpleNamespace(
                pickle_available=True,
                pickle_path=Path("/tmp/mole.pkl"),
                timber_available=True,
                timber_library_path=Path("/tmp/libtimber_model.so"),
                timber_model_dir=Path("/tmp/timber"),
                selected_backend="timber",
            )
            with mock.patch.object(mole_cli, "inspect_classifier_backends", return_value=fake_probe), mock.patch.object(
                mole_cli.torch.cuda, "is_available", return_value=True
            ), mock.patch.dict(mole_cli.os.environ, {"MOLE_MOLE_MODEL_PATH": str(mole_model)}):
                self.assertEqual(mole_cli._command_doctor(args), 0)

    def test_optimize_translates_to_chunked_workflow_command(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            template = tmp / "template.toml"
            env_file = tmp / "local.env"
            spec = tmp / "spec.json"
            scaffold = tmp / "scaffold.smi"
            for path in (template, env_file, spec, scaffold):
                path.write_text("placeholder", encoding="utf-8")

            args = SimpleNamespace(
                script=str(tmp / "run_long_rl.py"),
                template=str(template),
                env_file=str(env_file),
                experiment_spec=str(spec),
                scaffold_file=str(scaffold),
                output_root=str(tmp / "runs"),
                run_id="demo",
                dry_run=True,
            )
            completed = SimpleNamespace(returncode=0)
            with mock.patch.object(mole_cli.subprocess, "run", return_value=completed) as run_mock:
                self.assertEqual(mole_cli._command_optimize(args), 0)
            invoked = run_mock.call_args.args[0]
            self.assertIn(str(template.resolve()), invoked)
            self.assertIn(str(env_file.resolve()), invoked)
            self.assertIn(str(spec.resolve()), invoked)
            self.assertIn(str(scaffold.resolve()), invoked)
            self.assertIn("--dry-run", invoked)

    def test_command_score_fills_site_reward_scaffold_from_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            scaffold = tmp / "scaffold.smi"
            scaffold.write_text("[*:0]c1c([*:1])ccc([*:2])c1\n", encoding="utf-8")
            args = SimpleNamespace(
                objective_file=str(tmp / "objective.json"),
                scaffold_file=str(scaffold),
                smiles=["CCO"],
                chem_ids=None,
                aggregate_scores=False,
                app_threshold=0.04374140128493309,
                min_nkill=10,
                output=None,
            )
            objective = {
                "mode": "single_strain",
                "strain": "Akkermansia muciniphila (NT5021)",
                "app_threshold": 0.04374140128493309,
                "min_nkill": 10,
                "tau": 0.02,
                "site_reward": {
                    "enabled": True,
                    "range_min": 4,
                    "range_max": 27,
                    "lambda": 0.85,
                },
                "label": "single_akkermansia",
            }
            predictor = SimpleNamespace(
                ensure_loaded=mock.AsyncMock(),
                predict=mock.AsyncMock(return_value=[{"chem_id": "mol1", "score": 0.5}]),
            )
            captured = {}

            def fake_score_reinvent_predictions(raw_items, request):
                captured["scaffold_smiles"] = request.objective.site_reward.scaffold_smiles
                captured["mode"] = request.objective.mode
                return [{"chem_id": "mol1", "score": 0.5}]

            with mock.patch.object(mole_cli, "load_objective_spec", return_value=objective), mock.patch.object(
                mole_cli, "get_predictor", return_value=predictor
            ), mock.patch.object(mole_cli, "score_reinvent_predictions", side_effect=fake_score_reinvent_predictions), mock.patch.object(
                mole_cli, "_dump_json"
            ):
                self.assertEqual(mole_cli._command_score(args), 0)

            self.assertEqual(captured["scaffold_smiles"], scaffold.read_text(encoding="utf-8").strip())
            self.assertEqual(captured["mode"], "single_strain")

    def test_command_benchmark_screening_inputs_dumps_benchmark_json(self):
        args = SimpleNamespace(
            input_path=["/tmp/input_a.csv", "/tmp/input_b.parquet"],
            output=None,
        )
        expected = {
            "benchmarks": [
                {
                    "input_path": "/tmp/input_a.csv",
                    "elapsed_seconds": 0.01,
                    "size_bytes": 128,
                }
            ]
        }

        with mock.patch.object(mole_cli, "benchmark_paths", return_value=expected) as benchmark_mock, mock.patch.object(
            mole_cli, "_dump_json"
        ) as dump_mock:
            exit_code = mole_cli._command_benchmark_screening_inputs(args)

        self.assertEqual(exit_code, 0)
        benchmark_mock.assert_called_once_with(args.input_path)
        dump_mock.assert_called_once_with(expected, None)

    def test_predict_async_threads_graph_construction_settings(self):
        fake_scheduler = SimpleNamespace(
            predict_molecules=mock.AsyncMock(
                return_value=[{"chem_id": "mol1", "apscore_total": -1.2}]
            ),
            runtime_snapshot=mock.Mock(return_value={"device": "cpu"}),
        )
        args = SimpleNamespace(
            smiles=["CCO"],
            chem_ids=["mol1"],
            aggregate_scores=True,
            app_threshold=0.04374140128493309,
            min_nkill=10,
            num_graph_workers="3",
            graph_batch_size=128,
            prefetch_batches=4,
            classifier_workers="2",
            classifier_inflight_batches="5",
            deterministic_representation=True,
            classifier_backend="pickle",
            profiling=False,
        )

        with mock.patch.object(mole_cli, "get_scheduler", return_value=fake_scheduler) as get_scheduler_mock, mock.patch.dict(
            mole_cli.os.environ, {}, clear=False
        ):
            payload = mole_cli.asyncio_run(mole_cli._predict_async(args))
            self.assertEqual(mole_cli.os.environ["MOLE_CLASSIFIER_BACKEND"], "pickle")

        _, scheduler_kwargs = get_scheduler_mock.call_args
        self.assertEqual(scheduler_kwargs["classifier_workers"], "2")
        self.assertEqual(scheduler_kwargs["classifier_inflight_batches"], "5")
        _, kwargs = fake_scheduler.predict_molecules.call_args
        self.assertEqual(kwargs["num_graph_workers"], "3")
        self.assertEqual(kwargs["graph_batch_size"], 128)
        self.assertEqual(kwargs["prefetch_batches"], 4)
        self.assertEqual(kwargs["classifier_workers"], "2")
        self.assertEqual(kwargs["classifier_inflight_batches"], "5")
        self.assertTrue(kwargs["deterministic_representation"])
        self.assertEqual(payload["items"], [{"chem_id": "mol1", "apscore_total": -1.2}])

    def test_command_stream_enumeration_screen_threads_classifier_stage_settings(self):
        args = SimpleNamespace(
            output_dir="/tmp/out",
            scaffold_file="/tmp/scaffold.smi",
            scaffold_dir=None,
            scaffold_catalog=None,
            ordinary_library="/tmp/ordinary.csv",
            pos13_library="/tmp/pos13.csv",
            run_state_source=None,
            chunk_manifest_source=None,
            start_index=0,
            stop_index=10,
            shard_size=5,
            prediction_batch_size=2,
            classifier_backend="pickle",
            app_threshold=0.04374140128493309,
            min_nkill=10,
            num_graph_workers="3",
            graph_batch_size=128,
            prefetch_batches=4,
            classifier_workers="2",
            classifier_inflight_batches="5",
            enumeration_workers="auto",
            enumeration_prefetch_batches="auto",
            deterministic_representation=True,
            profiling=False,
        )
        fake_scheduler = SimpleNamespace(runtime_snapshot=mock.Mock(return_value={"device": "cuda:0"}))
        fake_summary = SimpleNamespace(
            output_dir=Path("/tmp/out"),
            run_state_path=Path("/tmp/out/run_state.json"),
            shard_manifest_path=Path("/tmp/out/shard_manifest.jsonl"),
            attempted_count=10,
            hit_count=3,
            completed_shards=2,
            start_index=0,
            stop_index=10,
            total_combinations=100,
        )

        async def fake_stream_enumeration_screen(**kwargs):
            self.assertEqual(kwargs["classifier_workers"], "2")
            self.assertEqual(kwargs["classifier_inflight_batches"], "5")
            return fake_summary

        with mock.patch.object(mole_cli, "create_scheduler", return_value=fake_scheduler) as create_scheduler_mock, mock.patch.object(
            mole_cli, "stream_enumeration_screen", side_effect=fake_stream_enumeration_screen
        ), mock.patch.object(mole_cli, "_dump_json"):
            self.assertEqual(mole_cli._command_stream_enumeration_screen(args), 0)

        _, scheduler_kwargs = create_scheduler_mock.call_args
        self.assertEqual(scheduler_kwargs["classifier_workers"], "2")
        self.assertEqual(scheduler_kwargs["classifier_inflight_batches"], "5")


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
