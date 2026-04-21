from __future__ import annotations

import tempfile
import tarfile
import sqlite3
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pandas as pd

import mole_cli
from src.batch_screening import load_screening_input


class BatchScreeningTestCase(unittest.TestCase):
    def test_load_screening_input_autogenerates_chem_id_when_missing(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "input.tsv"
            path.write_text("smiles\nCCO\nCCN\n", encoding="utf-8")

            frame = load_screening_input(path)

        self.assertEqual(frame["chem_id"].tolist(), ["mol1", "mol2"])
        self.assertEqual(frame["smiles"].tolist(), ["CCO", "CCN"])

    def test_load_screening_input_autogenerates_chem_id_inside_archive(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            csv_path = root / "tylosin_scheme_b_unique_products.csv"
            csv_path.write_text("product_smiles_canonical\nCCO\nCCN\n", encoding="utf-8")
            archive_path = root / "input.tar.gz"
            with tarfile.open(archive_path, "w:gz") as archive:
                archive.add(csv_path, arcname="2026-04-21_screening/tylosin/scheme_b_fix_pos13/tylosin_scheme_b_unique_products.csv")

            frame = load_screening_input(archive_path)

        self.assertEqual(frame["chem_id"].tolist(), ["tylosin_scheme_b_unique_products::1", "tylosin_scheme_b_unique_products::2"])
        self.assertEqual(frame["source_group"].tolist(), ["tylosin", "tylosin"])

    def test_load_screening_input_autogenerates_chem_id_inside_sqlite(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            db_path = root / "input.sqlite"
            with sqlite3.connect(db_path) as connection:
                connection.execute("CREATE TABLE molecules (smiles TEXT)")
                connection.executemany("INSERT INTO molecules (smiles) VALUES (?)", [("CCO",), ("CCN",)])
                connection.commit()

            frame = load_screening_input(db_path)

        self.assertEqual(frame["chem_id"].tolist(), ["molecules::1", "molecules::2"])
        self.assertEqual(frame["source_group"].tolist(), ["molecules", "molecules"])

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
                chunk_size=256,
                output_format="tsv",
            )
            predicted = pd.DataFrame(
                [{"chem_id": "mol1", "score": 0.7, "source_group": "input"}]
            )

            with mock.patch.object(mole_cli, "screen_path", return_value=(pd.DataFrame([{"smiles": "CCO", "chem_id": "mol1"}]), predicted)):
                exit_code = mole_cli._command_screen(args)

            self.assertEqual(exit_code, 0)
            self.assertTrue(output_dir.exists())
            self.assertTrue((output_dir / "predictions_all.tsv").is_file())
            self.assertTrue((output_dir / "normalized_input.tsv").is_file())
            self.assertTrue((output_dir / "manifest.json").is_file())


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
