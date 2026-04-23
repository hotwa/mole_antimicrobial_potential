import os
import sqlite3
import tarfile
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from src.screening_planner import ScreeningPlanner, PlannerConfig
from src.screening_sources import process_work_unit


class TestScreeningPlanner(unittest.TestCase):

    def setUp(self):
        self.tmpdir_obj = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tmpdir_obj.name)

    def tearDown(self):
        self.tmpdir_obj.cleanup()

    def test_auto_cpu_workers_selects_reasonable_value(self):
        # Auto mode should pick a value between 2 and 8
        planner = ScreeningPlanner(PlannerConfig(cpu_workers="auto"))
        self.assertGreaterEqual(planner.cpu_workers, 2)
        self.assertLessEqual(planner.cpu_workers, 8)

    def test_plan_tabular_large_file_splits_into_chunks(self):
        # Create a "large" dummy tabular file
        file_path = self.tmpdir / "large.tsv"
        with open(file_path, "w", encoding="utf-8") as f:
            f.write("smiles\tchem_id\n")
            for i in range(100):
                f.write(f"C{i}\tC{i}\n")

        # Force low limits to trigger chunking
        planner = ScreeningPlanner(PlannerConfig(
            grouping_mode="auto",
            target_rows_per_group=1, # 1 row per chunk to force maximum splits
            target_bytes_per_group=10
        ))

        units = planner.plan(file_path, archive_pattern="*.csv")
        self.assertGreater(len(units), 1)
        for unit in units:
            self.assertEqual(unit.source_type, "tabular")
            self.assertTrue(unit.group_id.startswith("large"))
            self.assertIsNotNone(unit.start_row)
            self.assertIsNotNone(unit.max_rows)

    def test_plan_tabular_chunking_does_not_truncate_short_line_file(self):
        file_path = self.tmpdir / "short_lines.tsv"
        with open(file_path, "w", encoding="utf-8") as f:
            f.write("smiles\tchem_id\n")
            for i in range(1000):
                f.write(f"C\tm{i}\n")

        planner = ScreeningPlanner(
            PlannerConfig(
                grouping_mode="chunk",
                target_rows_per_group=50,
                target_bytes_per_group=10**9,
            )
        )

        units = planner.plan(file_path, archive_pattern="*.csv")
        frames = []
        for unit in units:
            frames.extend(
                process_work_unit(
                    unit,
                    "smiles",
                    "chem_id",
                    "product_smiles_canonical",
                    "example_combo_id",
                    chunk_size=25,
                )
            )

        combined = pd.concat(frames, ignore_index=True)
        self.assertEqual(len(combined), 1000)
        self.assertEqual(combined["chem_id"].tolist()[0], "m0")
        self.assertEqual(combined["chem_id"].tolist()[-1], "m999")

    def test_plan_archive_identifies_members(self):
        archive_path = self.tmpdir / "bundle.tar.gz"
        with tarfile.open(archive_path, "w:gz") as tar:
            f1 = self.tmpdir / "m1_scheme_b_unique_products.csv"
            f1.write_text("smiles,chem_id\nA,A")
            f2 = self.tmpdir / "m2_scheme_b_unique_products.csv"
            f2.write_text("smiles,chem_id\nB,B")

            tar.add(f1, arcname=f1.name)
            tar.add(f2, arcname=f2.name)

        planner = ScreeningPlanner(PlannerConfig(
            grouping_mode="source",
            target_rows_per_group=100000,
            target_bytes_per_group=1000000
        ))

        units = planner.plan(archive_path, archive_pattern="*_scheme_b_unique_products.csv")
        self.assertEqual(len(units), 2)
        members = {u.source_member for u in units}
        self.assertEqual(members, {"m1_scheme_b_unique_products.csv", "m2_scheme_b_unique_products.csv"})
        for u in units:
            self.assertEqual(u.source_type, "archive")
            self.assertFalse("part" in u.group_id)

    def test_plan_archive_does_not_split_large_member(self):
        archive_path = self.tmpdir / "bundle_large.tar.gz"
        with tarfile.open(archive_path, "w:gz") as tar:
            f1 = self.tmpdir / "huge_scheme_b_unique_products.csv"
            # 100 lines
            f1.write_text("smiles,chem_id\n" + "\n".join(f"C{i},C{i}" for i in range(100)))
            tar.add(f1, arcname=f1.name)

        planner = ScreeningPlanner(PlannerConfig(
            grouping_mode="auto", # auto might try to split
            target_rows_per_group=1,
            target_bytes_per_group=10
        ))

        units = planner.plan(archive_path, archive_pattern="*_scheme_b_unique_products.csv")
        # should NOT be split into multiple chunks
        self.assertEqual(len(units), 1)
        u = units[0]
        self.assertEqual(u.source_type, "archive")
        self.assertEqual(u.source_member, "huge_scheme_b_unique_products.csv")
        self.assertIsNone(u.start_row)
        self.assertIsNone(u.max_rows)

    def test_plan_sqlite_query_no_nameerror(self):
        db_path = self.tmpdir / "data.sqlite"
        with sqlite3.connect(db_path) as conn:
            conn.execute("CREATE TABLE t1 (smiles TEXT, chem_id TEXT)")
            conn.executemany("INSERT INTO t1 VALUES (?, ?)", [(f"C{i}", str(i)) for i in range(25)])
            conn.commit()

        planner = ScreeningPlanner(PlannerConfig(
            grouping_mode="chunk",
            target_rows_per_group=10
        ))

        units = planner.plan(str(db_path), archive_pattern="", sqlite_query="SELECT * FROM t1", smiles_colname="smiles")
        # 25 rows, 10 per chunk => 3 chunks. The test proves pandas NameError is fixed.
        self.assertEqual(len(units), 3)
        for i, u in enumerate(units):
            self.assertEqual(u.source_type, "sqlite")
            self.assertEqual(u.source_member, "query")
            self.assertEqual(u.start_row, i * 9)

    def test_plan_sqlite_splits_large_table(self):
        db_path = self.tmpdir / "data.sqlite"
        with sqlite3.connect(db_path) as conn:
            conn.execute("CREATE TABLE t1 (smiles TEXT, chem_id TEXT)")
            conn.executemany("INSERT INTO t1 VALUES (?, ?)", [(f"C{i}", str(i)) for i in range(50)])
            conn.commit()

        planner = ScreeningPlanner(PlannerConfig(
            grouping_mode="chunk",
            target_rows_per_group=10
        ))

        units = planner.plan(db_path, archive_pattern="", sqlite_table="t1")
        # 50 rows, 10 per chunk => 5 chunks
        self.assertEqual(len(units), 5)
        for i, u in enumerate(units):
            self.assertEqual(u.source_type, "sqlite")
            self.assertEqual(u.source_member, "t1")
            self.assertEqual(u.start_row, i * 10)
            self.assertEqual(u.max_rows, 10)

    def test_plan_parquet_file_emits_single_parquet_work_unit(self):
        shard = self.tmpdir / "shard_0001.parquet"
        pd.DataFrame(
            {
                "smiles": ["CCO", "CCN"],
                "chem_id": ["m1", "m2"],
                "source_group": ["demo", "demo"],
            }
        ).to_parquet(shard, index=False, engine="pyarrow", row_group_size=1)

        planner = ScreeningPlanner(PlannerConfig(grouping_mode="auto"))
        units = planner.plan(shard, archive_pattern="*.csv")

        self.assertEqual(len(units), 1)
        self.assertEqual(units[0].source_type, "parquet")
        self.assertEqual(units[0].source_path, str(shard.resolve()))
        self.assertEqual(units[0].source_group, "demo")

    def test_plan_parquet_file_without_source_group_uses_file_stem(self):
        shard = self.tmpdir / "single_input.parquet"
        pd.DataFrame(
            {
                "smiles": ["CCO", "CCN"],
                "chem_id": ["m1", "m2"],
            }
        ).to_parquet(shard, index=False, engine="pyarrow", row_group_size=1)

        planner = ScreeningPlanner(PlannerConfig(grouping_mode="auto"))
        units = planner.plan(shard, archive_pattern="*.csv")

        self.assertEqual(len(units), 1)
        self.assertEqual(units[0].source_type, "parquet")
        self.assertEqual(units[0].source_group, "single_input")

    def test_plan_parquet_directory_emits_one_work_unit_per_shard(self):
        shard_dir = self.tmpdir / "prepared" / "demo"
        shard_dir.mkdir(parents=True)
        for idx, smiles in enumerate(["CCO", "CCN"], start=1):
            pd.DataFrame(
                {
                    "smiles": [smiles],
                    "chem_id": [f"m{idx}"],
                    "source_group": ["demo"],
                }
            ).to_parquet(
                shard_dir / f"shard_{idx:04d}.parquet",
                index=False,
                engine="pyarrow",
                row_group_size=1,
            )

        planner = ScreeningPlanner(PlannerConfig(grouping_mode="auto"))
        units = planner.plan(shard_dir, archive_pattern="*.csv")

        self.assertEqual(len(units), 2)
        self.assertTrue(all(unit.source_type == "parquet" for unit in units))
        self.assertEqual({unit.source_group for unit in units}, {"demo"})

    def test_plan_parquet_directory_without_source_group_uses_parent_directory(self):
        shard_dir = self.tmpdir / "prepared" / "demo_no_column"
        shard_dir.mkdir(parents=True)
        pd.DataFrame(
            {
                "smiles": ["CCO"],
                "chem_id": ["m1"],
            }
        ).to_parquet(
            shard_dir / "shard_0001.parquet",
            index=False,
            engine="pyarrow",
            row_group_size=1,
        )

        planner = ScreeningPlanner(PlannerConfig(grouping_mode="auto"))
        units = planner.plan(shard_dir, archive_pattern="*.csv")

        self.assertEqual(len(units), 1)
        self.assertEqual(units[0].source_type, "parquet")
        self.assertEqual(units[0].source_group, "demo_no_column")

if __name__ == "__main__":
    unittest.main()
