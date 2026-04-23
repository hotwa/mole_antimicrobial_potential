import tempfile
import tarfile
import sqlite3
import unittest
from pathlib import Path
import pandas as pd

from src.screening_planner import WorkUnit
from src.screening_sources import process_work_unit, stream_screening_input


class ScreeningSourcesTestCase(unittest.TestCase):
    def test_load_screening_input_reads_parquet_chunks(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            shard = Path(tmpdir) / "demo.parquet"
            pd.DataFrame(
                {
                    "smiles": ["CCO", "CCN", "CCC"],
                    "chem_id": ["m1", "m2", "m3"],
                    "source_group": ["demo", "demo", "demo"],
                }
            ).to_parquet(shard, index=False, engine="pyarrow", row_group_size=2)

            frames = list(stream_screening_input(shard, chunk_size=1))

        self.assertEqual(len(frames), 3)
        self.assertEqual([len(frame) for frame in frames], [1, 1, 1])
        self.assertEqual(
            pd.concat(frames, ignore_index=True)["source_row"].tolist(),
            [1, 2, 3],
        )
        self.assertTrue(all(frame.iloc[0]["source_file"] == str(shard) for frame in frames))
        self.assertTrue(all(frame.iloc[0]["source_group"] == "demo" for frame in frames))

    def test_load_screening_input_autogenerates_chem_id_when_missing(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "input.tsv"
            path.write_text("smiles\nCCO\nCCN\n", encoding="utf-8")

            frames = list(stream_screening_input(path))
            frame = pd.concat(frames)

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

            frames = list(stream_screening_input(archive_path))
            frame = pd.concat(frames)

        self.assertEqual(frame["chem_id"].tolist(), ["tylosin_scheme_b_unique_products__1", "tylosin_scheme_b_unique_products__2"])
        self.assertEqual(frame["source_group"].tolist(), ["tylosin", "tylosin"])

    def test_load_screening_input_autogenerates_chem_id_inside_sqlite(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            db_path = root / "input.sqlite"
            with sqlite3.connect(db_path) as connection:
                connection.execute("CREATE TABLE molecules (smiles TEXT)")
                connection.executemany("INSERT INTO molecules (smiles) VALUES (?)", [("CCO",), ("CCN",)])
                connection.commit()

            frames = list(stream_screening_input(db_path))
            frame = pd.concat(frames)

        self.assertEqual(frame["chem_id"].tolist(), ["molecules__1", "molecules__2"])
        self.assertEqual(frame["source_group"].tolist(), ["molecules", "molecules"])

    def test_tabular_input_closes_file_on_early_exit(self) -> None:
        import warnings
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "large_input.tsv"
            # Write enough data so it reads in multiple chunks
            lines = ["smiles\n"] + [f"C{i}C\n" for i in range(100)]
            path.write_text("".join(lines), encoding="utf-8")

            # Request chunk_size=10 so it yields 10 chunks
            # Catch warnings to ensure no ResourceWarning is emitted
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always", ResourceWarning)

                gen = stream_screening_input(path, chunk_size=10)
                # Read just the first chunk and then we exit
                next(gen)
                # Close the generator explicitly, simulating what generator garbage collection
                # or a break in the loop would initiate (which triggers GeneratorExit)
                gen.close()

                for warning in w:
                    if issubclass(warning.category, ResourceWarning):
                        self.fail(f"ResourceWarning was emitted: {warning.message}")

    def test_process_work_unit_reads_parquet_chunks(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            shard = Path(tmpdir) / "shard_0001.parquet"
            pd.DataFrame(
                {
                    "smiles": ["CCO", "CCN"],
                    "chem_id": ["m1", "m2"],
                    "source_group": ["demo", "demo"],
                }
            ).to_parquet(shard, index=False, engine="pyarrow", row_group_size=1)

            unit = WorkUnit(
                source_type="parquet",
                source_path=str(shard),
                group_id="demo",
                source_group="demo",
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

        self.assertEqual(len(frames), 2)
        self.assertEqual([frame.iloc[0]["chem_id"] for frame in frames], ["m1", "m2"])
        self.assertTrue(all(frame.iloc[0]["source_group"] == "demo" for frame in frames))

    def test_process_work_unit_autogenerates_chem_id_for_parquet(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            shard = Path(tmpdir) / "shard_0007.parquet"
            pd.DataFrame(
                {
                    "smiles": ["CCO", "CCN"],
                    "source_group": ["demo", "demo"],
                }
            ).to_parquet(shard, index=False, engine="pyarrow", row_group_size=1)

            unit = WorkUnit(
                source_type="parquet",
                source_path=str(shard),
                group_id="demo_part1",
                source_group="demo",
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
        self.assertEqual(combined["chem_id"].tolist(), ["shard_0007__1", "shard_0007__2"])
        self.assertEqual(combined["source_group"].tolist(), ["demo", "demo"])
        self.assertEqual(combined["source_file"].tolist(), [str(shard), str(shard)])
        self.assertEqual(combined["source_row"].tolist(), [1, 2])


if __name__ == "__main__":
    unittest.main()
