from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
import os
from unittest import mock

import numpy as np

import src.classifier_backend as classifier_backend
from src.classifier_backend import (
    DEFAULT_PICKLE_MODEL_PATH,
    DEFAULT_TIMBER_MODEL_DIR,
    PickleXGBoostBackend,
    TimberCompiledArtifactBackend,
    inspect_classifier_backends,
    resolve_classifier_backend,
)


class TestClassifierBackend(unittest.TestCase):
    def test_auto_prefers_pickle_when_pickle_is_available(self):
        probe = inspect_classifier_backends()
        if not probe.pickle_available:
            self.skipTest("Pickle model is not available on this machine")

        backend = resolve_classifier_backend()
        self.assertEqual(backend.name, "pickle")

    def test_auto_falls_back_to_timber_when_pickle_missing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            timber_dir = Path(tmpdir)
            (timber_dir / "compiled").mkdir(parents=True)
            timber_lib = timber_dir / "compiled" / "libtimber_model.so"
            timber_lib.write_bytes(b"stub")

            probe = inspect_classifier_backends(
                pickle_path=timber_dir / "missing.pkl",
                timber_model_dir=timber_dir,
            )

        self.assertFalse(probe.pickle_available)
        self.assertTrue(probe.timber_available)
        self.assertEqual(probe.selected_backend, "timber")

    def test_auto_falls_back_to_pickle_when_timber_missing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            probe = inspect_classifier_backends(
                pickle_path=DEFAULT_PICKLE_MODEL_PATH,
                timber_model_dir=tmpdir,
            )
            self.assertTrue(probe.pickle_available)
            self.assertFalse(probe.timber_available)
            self.assertEqual(probe.selected_backend, "pickle")

            backend = resolve_classifier_backend(
                pickle_path=DEFAULT_PICKLE_MODEL_PATH,
                timber_model_dir=tmpdir,
            )
            self.assertEqual(backend.name, "pickle")

    def test_auto_falls_back_to_pickle_when_timber_initialization_fails(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            timber_dir = Path(tmpdir)
            (timber_dir / "compiled").mkdir(parents=True)
            timber_lib = timber_dir / "compiled" / "libtimber_model.so"
            timber_lib.write_bytes(b"stub")
            with mock.patch.object(
                classifier_backend,
                "TimberCompiledArtifactBackend",
                side_effect=OSError("boom"),
            ):
                backend = resolve_classifier_backend(
                    pickle_path=DEFAULT_PICKLE_MODEL_PATH,
                    timber_model_dir=timber_dir,
                )
            self.assertEqual(backend.name, "pickle")

    def test_default_pickle_path_is_repo_root_relative(self):
        repo_root = Path(__file__).resolve().parents[1]
        expected = repo_root / "data" / "03.model_evaluation" / "MolE-XGBoost-08.03.2024_14.20.pkl"
        with tempfile.TemporaryDirectory() as tmpdir:
            with mock.patch.object(os, "getcwd", return_value=tmpdir):
                probe = inspect_classifier_backends(timber_model_dir=tmpdir)
        self.assertEqual(probe.pickle_path, expected)

    def test_timber_and_pickle_backends_match_on_random_features(self):
        if not DEFAULT_TIMBER_MODEL_DIR.is_dir():
            self.skipTest("Timber model cache is not available on this machine")
        if not (DEFAULT_TIMBER_MODEL_DIR / "lib" / "timber_model.so").is_file():
            self.skipTest("Timber compiled artifact is not available on this machine")

        rng = np.random.default_rng(0)
        X = rng.random((4, 1040), dtype=np.float32)

        pickle_backend = PickleXGBoostBackend(DEFAULT_PICKLE_MODEL_PATH)
        timber_backend = TimberCompiledArtifactBackend(DEFAULT_TIMBER_MODEL_DIR)

        pickle_probs = pickle_backend.predict_proba(X)
        timber_probs = timber_backend.predict_proba(X)

        self.assertEqual(pickle_probs.shape, (4, 2))
        self.assertEqual(timber_probs.shape, (4, 2))
        np.testing.assert_allclose(pickle_probs, timber_probs, rtol=2e-2, atol=2e-2)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
