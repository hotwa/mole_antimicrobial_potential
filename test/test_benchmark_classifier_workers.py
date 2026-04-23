from __future__ import annotations

import unittest

from scripts import benchmark_classifier_workers as bench


class BenchmarkClassifierWorkersHelpersTest(unittest.TestCase):
    def test_parse_worker_value_accepts_auto_and_int(self) -> None:
        self.assertEqual(bench._parse_worker_value("auto"), "auto")
        self.assertEqual(bench._parse_worker_value("6"), 6)

    def test_parse_worker_value_rejects_non_positive_int(self) -> None:
        with self.assertRaises(ValueError):
            bench._parse_worker_value("0")

    def test_summarize_candidate_runs_uses_median_wall_time(self) -> None:
        summary = bench._summarize_candidate_runs(
            mode="predictor",
            worker=6,
            sample_size=100,
            runs=[
                {"wall_seconds": 3.0},
                {"wall_seconds": 1.0},
                {"wall_seconds": 2.0},
            ],
        )

        self.assertEqual(summary["mode"], "predictor")
        self.assertEqual(summary["worker"], 6)
        self.assertEqual(summary["sample_size"], 100)
        self.assertEqual(summary["median_wall_seconds"], 2.0)
        self.assertEqual(summary["median_molecules_per_second"], 50.0)

    def test_select_best_candidate_prefers_higher_throughput(self) -> None:
        best = bench._select_best_candidate(
            [
                {"worker": 2, "median_molecules_per_second": 500.0, "median_wall_seconds": 2.0},
                {"worker": 6, "median_molecules_per_second": 700.0, "median_wall_seconds": 1.5},
                {"worker": 8, "median_molecules_per_second": 650.0, "median_wall_seconds": 1.7},
            ]
        )

        self.assertEqual(best["worker"], 6)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
