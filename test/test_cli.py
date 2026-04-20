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


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
