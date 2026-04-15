#!/usr/bin/env python3
"""Unit tests for REINVENT4 workflow helpers."""

import tempfile
import unittest
from importlib import util
from pathlib import Path
from unittest import mock

from src.reinvent4_workflow import (
    build_reinvent_command,
    build_reinvent_context,
    build_score_request,
    chunk_reward_summary,
    create_run_directory,
    external_process_payload,
    extract_top_unique_smiles,
    load_long_run_spec,
    load_objective_spec,
    load_runtime_settings,
    render_template_file,
    unique_murcko_count,
    validate_scaffold_file,
    validate_scaffold_smiles,
)


class Reinvent4WorkflowTestCase(unittest.TestCase):
    def test_validate_scaffold_supports_four_attachment_points(self) -> None:
        scaffold = "[*:0]c1c([*:1])c([*:2])cc([*:3])c1"
        labels = validate_scaffold_smiles(scaffold)
        self.assertEqual(labels, [0, 1, 2, 3])

    def test_validate_scaffold_rejects_non_contiguous_labels(self) -> None:
        with self.assertRaisesRegex(ValueError, "contiguous set starting at 0"):
            validate_scaffold_smiles("[*:0]c1c([*:2])ccc([*:3])c1")

    def test_load_objective_spec_normalizes_defaults(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "objective.json"
            path.write_text('{"mode":"single_strain","strain":"Akkermansia muciniphila (NT5021)"}', encoding="utf-8")
            objective = load_objective_spec(path)
        self.assertEqual(objective["mode"], "single_strain")
        self.assertEqual(objective["strain"], "Akkermansia muciniphila (NT5021)")
        self.assertEqual(objective["app_threshold"], 0.04374140128493309)
        self.assertEqual(objective["min_nkill"], 10)
        self.assertEqual(objective["tau"], 0.02)

    def test_load_objective_spec_resolves_single_strain_index(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "objective.json"
            path.write_text('{"mode":"single_strain","strain_index":1}', encoding="utf-8")
            objective = load_objective_spec(path)
        self.assertEqual(objective["strain"], "Akkermansia muciniphila (NT5021)")
        self.assertNotIn("strain_index", objective)

    def test_load_objective_spec_resolves_multiple_strain_indices(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "objective.json"
            path.write_text(
                '{"mode":"single_strain","strain_indices":[1,6,16],"weights":[0.5,0.2,0.3]}',
                encoding="utf-8",
            )
            objective = load_objective_spec(path)
        self.assertEqual(
            objective["strains"],
            [
                "Akkermansia muciniphila (NT5021)",
                "Bacteroides thetaiotaomicron (NT5004)",
                "Clostridium perfringens (NT5032)",
            ],
        )
        self.assertEqual(objective["weights"], [0.5, 0.2, 0.3])

    def test_load_objective_spec_rejects_unknown_strain_index(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "objective.json"
            path.write_text('{"mode":"single_strain","strain_index":99}', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Unknown strain index"):
                load_objective_spec(path)

    def test_load_long_run_spec_resolves_monitor_panel(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "long_run.json"
            path.write_text(
                (
                    '{"experiment_name":"Pathogen Group A",'
                    '"objective":{"mode":"single_strain","strain_indices":[15,16,3,27],"weights":[0.3,0.3,0.2,0.2]},'
                    '"chunk_steps":25,"min_total_steps":100,"max_total_steps":400,'
                    '"success_high":0.8,"success_top5_mean":0.6,"success_consecutive_chunks":3,'
                    '"monitor_top_n":20,"commensal_monitor_indices":[1,10,11]}'
                ),
                encoding="utf-8",
            )
            spec = load_long_run_spec(path)
        self.assertEqual(spec["min_chunks"], 4)
        self.assertEqual(spec["max_chunks"], 16)
        self.assertEqual(
            spec["objective"]["strains"],
            [
                "Clostridium difficile (NT5083)",
                "Clostridium perfringens (NT5032)",
                "Bacteroides fragilis (ET) (NT5033)",
                "Fusobacterium nucleatum (NT5025)",
            ],
        )
        self.assertEqual(
            spec["commensal_monitor_strains"],
            [
                "Akkermansia muciniphila (NT5021)",
                "Bifidobacterium adolescentis (NT5022)",
                "Bifidobacterium longum (NT5028)",
            ],
        )

    def test_load_long_run_spec_rejects_non_divisible_steps(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "long_run.json"
            path.write_text(
                (
                    '{"experiment_name":"Bad Steps",'
                    '"objective":{"mode":"single_strain","strain_index":15},'
                    '"chunk_steps":25,"min_total_steps":90,"max_total_steps":300}'
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "divisible by 'chunk_steps'"):
                load_long_run_spec(path)

    def test_build_score_request_uses_expected_chem_ids(self) -> None:
        objective = {
            "mode": "single_strain",
            "strain": "Akkermansia muciniphila (NT5021)",
            "app_threshold": 0.04374140128493309,
            "min_nkill": 10,
            "tau": 0.02,
            "label": "single_akkermansia",
        }
        payload = build_score_request(["CCO", "CCN"], objective)
        self.assertEqual(payload["chem_id"], ["mol1", "mol2"])
        self.assertEqual(payload["objective"]["strain"], objective["strain"])

    def test_external_process_payload_preserves_input_order(self) -> None:
        response = {
            "objective": {"mode": "single_strain"},
            "items": [
                {"chem_id": "mol2", "score": 0.2},
                {"chem_id": "mol1", "score": 0.8},
            ],
        }
        payload = external_process_payload(response, ["mol1", "mol2"])
        self.assertEqual(payload["payload"]["predictions"], [0.8, 0.2])
        self.assertEqual(payload["payload"]["chem_ids"], ["mol1", "mol2"])

    def test_extract_top_unique_smiles_prefers_high_scores(self) -> None:
        rows = [
            {"SMILES": "CCO", "total_score": "0.2", "step": "1"},
            {"SMILES": "CCN", "total_score": "0.9", "step": "1"},
            {"SMILES": "CCO", "total_score": "0.8", "step": "2"},
        ]
        selected = extract_top_unique_smiles(rows, top_n=2)
        self.assertEqual([item["smiles"] for item in selected], ["CCN", "CCO"])
        self.assertEqual(selected[0]["reinvent_total_score"], 0.9)

    def test_chunk_reward_summary_and_unique_murcko_count(self) -> None:
        rows = [
            {"SMILES": "CCO", "Score": "0.0", "Antimicrobial Reward": "0.0"},
            {"SMILES": "CCN", "Score": "0.8", "Antimicrobial Reward": "0.8"},
            {"SMILES": "CCN", "Score": "0.6", "Antimicrobial Reward": "0.6"},
            {"SMILES": "c1ccccc1O", "Score": "0.4", "Antimicrobial Reward": "0.4"},
        ]
        summary = chunk_reward_summary(rows)
        self.assertAlmostEqual(summary["chunk_max_reward"], 0.8)
        self.assertAlmostEqual(summary["chunk_top5_mean_reward"], 0.45)
        self.assertAlmostEqual(summary["zero_score_fraction"], 0.25)
        self.assertAlmostEqual(summary["unique_smiles_ratio"], 0.75)
        self.assertGreaterEqual(unique_murcko_count(["CCN", "c1ccccc1O"]), 1)

    def test_render_context_contains_bridge_script_and_objective_file(self) -> None:
        settings = {
            "prior_file": Path("/tmp/libinvent.prior"),
            "agent_file": Path("/tmp/libinvent.prior"),
            "score_url": "http://127.0.0.1:8000/score",
            "device": "cpu",
            "batch_size": 64,
            "sigma": 128,
            "learning_rate": 0.0001,
            "bucket_size": 25,
            "diversity_minscore": 0.4,
            "max_score": 0.8,
            "min_steps": 25,
            "max_steps": 500,
            "request_timeout_seconds": 120,
        }
        objective = {
            "label": "single_demo",
            "mode": "single_strain",
            "strain": "Akkermansia muciniphila (NT5021)",
            "app_threshold": 0.04374140128493309,
            "min_nkill": 10,
            "tau": 0.02,
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir) / "run"
            run_dir.mkdir()
            objective_file = Path(tmpdir) / "objective.json"
            objective_file.write_text("{}", encoding="utf-8")
            scaffold_file = Path(tmpdir) / "scaffold.smi"
            scaffold_file.write_text("[*:0]CC[*:1]\n", encoding="utf-8")
            context = build_reinvent_context(settings, objective, scaffold_file, objective_file, run_dir)
            template = Path(tmpdir) / "template.toml.tpl"
            template.write_text("params.args = ${score_bridge_args}\nname = ${score_component_name}\n", encoding="utf-8")
            rendered = render_template_file(template, context)
        self.assertIn("score_bridge.py", rendered)
        self.assertIn("objective.json", rendered)
        self.assertIn("Antimicrobial Reward", rendered)

    def test_load_runtime_settings_and_build_command(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            reinvent_root = root / "REINVENT4"
            prior = root / "libinvent.prior"
            reinvent_root.mkdir()
            prior.write_text("prior", encoding="utf-8")
            env_file = root / "local.env"
            env_file.write_text(
                "\n".join(
                    [
                        f"REINVENT4_ROOT={reinvent_root}",
                        "REINVENT4_RUNNER=python -m reinvent.Reinvent",
                        f"REINVENT4_PRIOR_FILE={prior}",
                        "MOLE_SCORE_API_URL=http://127.0.0.1:8000/score",
                    ]
                ),
                encoding="utf-8",
            )
            settings = load_runtime_settings(env_file)
            command = build_reinvent_command(settings, root / "reinvent.toml")
        self.assertEqual(command[:3], ["python", "-m", "reinvent.Reinvent"])
        self.assertTrue(command[-1].endswith("reinvent.toml"))

    def test_validate_scaffold_file_skips_comments(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            scaffold_file = Path(tmpdir) / "scaffolds.smi"
            scaffold_file.write_text("# comment\n\n[*:0]CC[*:1]\n", encoding="utf-8")
            records = validate_scaffold_file(scaffold_file)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["attachment_count"], 2)

    def test_create_run_directory_uses_sanitized_experiment_name(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = create_run_directory(tmpdir, "Single Strain Demo", run_id="fixed")
        self.assertTrue(str(run_dir).endswith("single_strain_demo/fixed"))

    def test_run_long_rl_dry_run_writes_manifest(self) -> None:
        module_path = (
            Path(__file__).resolve().parents[1]
            / "workflows"
            / "reinvent4"
            / "scripts"
            / "run_long_rl.py"
        )
        spec = util.spec_from_file_location("run_long_rl", module_path)
        assert spec is not None and spec.loader is not None
        module = util.module_from_spec(spec)
        spec.loader.exec_module(module)

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            reinvent_root = root / "REINVENT4"
            prior = root / "libinvent.prior"
            reinvent_root.mkdir()
            prior.write_text("prior", encoding="utf-8")
            env_file = root / "local.env"
            env_file.write_text(
                "\n".join(
                    [
                        f"REINVENT4_ROOT={reinvent_root}",
                        "REINVENT4_RUNNER=python -m reinvent.Reinvent",
                        f"REINVENT4_PRIOR_FILE={prior}",
                        "MOLE_SCORE_API_URL=http://127.0.0.1:8000/score",
                    ]
                ),
                encoding="utf-8",
            )
            long_run_file = root / "long_run.json"
            long_run_file.write_text(
                (
                    '{"experiment_name":"Pathogen Group A",'
                    '"objective":{"mode":"single_strain","strain_indices":[15,16,3,27],"weights":[0.3,0.3,0.2,0.2]},'
                    '"chunk_steps":25,"min_total_steps":100,"max_total_steps":400,'
                    '"success_high":0.8,"success_top5_mean":0.6,"success_consecutive_chunks":3,'
                    '"monitor_top_n":20,"commensal_monitor_indices":[1,10,11]}'
                ),
                encoding="utf-8",
            )
            output_root = root / "runs"

            argv = [
                "run_long_rl.py",
                "--template",
                str(
                    Path(__file__).resolve().parents[1]
                    / "workflows"
                    / "reinvent4"
                    / "configs"
                    / "templates"
                    / "multi_strain_rl_macrocycle.toml.tpl"
                ),
                "--env-file",
                str(env_file),
                "--experiment-spec",
                str(long_run_file),
                "--scaffold-file",
                str(
                    Path(__file__).resolve().parents[1]
                    / "workflows"
                    / "reinvent4"
                    / "inputs"
                    / "scaffolds"
                    / "mother_scaffold.template.smi"
                ),
                "--output-root",
                str(output_root),
                "--run-id",
                "dryrun",
                "--dry-run",
            ]
            with mock.patch.object(module, "http_health_check", return_value={"status": "ready"}):
                with mock.patch("sys.argv", argv):
                    exit_code = module.main()

            self.assertEqual(exit_code, 0)
            run_dir = output_root / "pathogen_group_a" / "dryrun"
            manifest = (run_dir / "run_manifest.json").read_text(encoding="utf-8")
            self.assertIn("dry_run_rendered_first_chunk", manifest)
            self.assertTrue((run_dir / "chunk_001" / "reinvent.toml").is_file())


if __name__ == "__main__":
    unittest.main()
