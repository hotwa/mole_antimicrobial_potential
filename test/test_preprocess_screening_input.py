from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

import mole_cli
from src.screening_planner import WorkUnit
from src.screening_sources import process_work_unit


class PreprocessScreeningInputTestCase(unittest.TestCase):
    def test_preprocess_csv_to_parquet_shards_keeps_required_columns(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_path = tmp / "input.csv"
            output_dir = tmp / "prepared"
            input_path.write_text(
                "product_smiles_canonical,example_combo_id,unused\n"
                "CCO,m1,x\n"
                "CCN,m2,y\n",
                encoding="utf-8",
            )

            from src.preprocess_screening_input import preprocess_to_parquet

            manifest = preprocess_to_parquet(
                input_path=input_path,
                output_dir=output_dir,
                smiles_colname="product_smiles_canonical",
                chem_id_colname="example_combo_id",
                source_group="demo",
                rows_per_shard=1,
                row_group_size=1,
            )

            self.assertEqual(manifest["shard_count"], 2)
            self.assertEqual(manifest["rows_written"], 2)

            shard_path = output_dir / "demo" / "shard_0001.parquet"
            self.assertTrue(shard_path.is_file())

            shard = pd.read_parquet(shard_path)
            self.assertEqual(
                list(shard.columns),
                ["smiles", "chem_id", "source_group", "source_file", "source_row"],
            )
            self.assertEqual(
                shard.to_dict(orient="records"),
                [
                    {
                        "smiles": "CCO",
                        "chem_id": "m1",
                        "source_group": "demo",
                        "source_file": str(input_path.resolve()),
                        "source_row": 1,
                    }
                ],
            )

    def test_preprocess_multiple_csvs_writes_one_source_group_per_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            a_path = tmp / "tylosin.csv"
            b_path = tmp / "tilmicosin.csv"
            output_dir = tmp / "prepared"
            a_path.write_text("product_smiles_canonical,example_combo_id\nCCO,t1\n", encoding="utf-8")
            b_path.write_text("product_smiles_canonical,example_combo_id\nCCN,m1\n", encoding="utf-8")

            from src.preprocess_screening_input import preprocess_many_to_parquet

            manifest = preprocess_many_to_parquet(
                input_paths=[a_path, b_path],
                output_dir=output_dir,
                smiles_colname="product_smiles_canonical",
                chem_id_colname="example_combo_id",
                rows_per_shard=1,
                row_group_size=1,
            )

            self.assertEqual(manifest["source_groups"], ["tilmicosin", "tylosin"])
            self.assertTrue((output_dir / "tylosin" / "shard_0001.parquet").is_file())
            self.assertTrue((output_dir / "tilmicosin" / "shard_0001.parquet").is_file())

    def test_process_work_unit_reads_parquet_shard_metadata_columns(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_path = tmp / "demo.csv"
            output_dir = tmp / "prepared"
            input_path.write_text(
                "product_smiles_canonical,example_combo_id,unused\n"
                "CCO,m1,x\n"
                "CCN,m2,y\n",
                encoding="utf-8",
            )

            from src.preprocess_screening_input import preprocess_to_parquet

            preprocess_to_parquet(
                input_path=input_path,
                output_dir=output_dir,
                smiles_colname="product_smiles_canonical",
                chem_id_colname="example_combo_id",
                source_group="demo",
                rows_per_shard=10,
                row_group_size=1,
            )

            shard_path = output_dir / "demo" / "shard_0001.parquet"
            unit = WorkUnit(
                source_type="parquet",
                source_path=str(shard_path),
                group_id="demo",
                source_group="fallback-group",
            )

            frames = list(
                process_work_unit(
                    unit,
                    "smiles",
                    "chem_id",
                    "product_smiles_canonical",
                    "example_combo_id",
                    chunk_size=1,
                )
            )

        combined = pd.concat(frames, ignore_index=True)
        self.assertEqual(combined["chem_id"].tolist(), ["m1", "m2"])
        self.assertEqual(combined["source_group"].tolist(), ["demo", "demo"])
        self.assertEqual(combined["source_file"].tolist(), [str(input_path.resolve()), str(input_path.resolve())])
        self.assertEqual(combined["source_row"].tolist(), [1, 2])

    def test_preprocess_generates_chem_id_when_missing_for_tsv(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_path = tmp / "input.tsv"
            output_dir = tmp / "prepared"
            input_path.write_text("smiles\tunused\nCCO\tx\nCCN\ty\n", encoding="utf-8")

            from src.preprocess_screening_input import preprocess_to_parquet

            manifest = preprocess_to_parquet(
                input_path=input_path,
                output_dir=output_dir,
                smiles_colname="smiles",
                chem_id_colname="chem_id",
                source_group="demo",
                rows_per_shard=10,
                row_group_size=2,
            )

            self.assertEqual(manifest["rows_written"], 2)
            self.assertEqual(manifest["shard_count"], 1)

            shard = pd.read_parquet(output_dir / "demo" / "shard_0001.parquet")
            self.assertEqual(list(shard["chem_id"]), ["demo__1", "demo__2"])

    def test_cli_parser_exposes_preprocess_screening_input(self) -> None:
        parser = mole_cli._build_parser()
        args = parser.parse_args(
            [
                "preprocess-screening-input",
                "--input-path",
                "input.csv",
                "--output-dir",
                "prepared",
            ]
        )
        self.assertEqual(args.command, "preprocess-screening-input")

    def test_command_preprocess_writes_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_path = tmp / "input.csv"
            output_dir = tmp / "prepared"
            input_path.write_text("smiles\nCCO\n", encoding="utf-8")

            args = mock.Mock(
                input_path=str(input_path),
                output_dir=str(output_dir),
                smiles_colname="smiles",
                chem_id_colname="chem_id",
                source_group="demo",
                rows_per_shard=100,
                row_group_size=10,
            )

            with mock.patch.object(mole_cli, "_dump_json") as dump_json:
                exit_code = mole_cli._command_preprocess_screening_input(args)

            self.assertEqual(exit_code, 0)
            dump_json.assert_called_once()


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
