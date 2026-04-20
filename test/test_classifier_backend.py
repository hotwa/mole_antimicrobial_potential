from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from src.classifier_backend import (
    DEFAULT_PICKLE_MODEL_PATH,
    DEFAULT_TIMBER_MODEL_DIR,
    PickleXGBoostBackend,
    TimberCompiledArtifactBackend,
    inspect_classifier_backends,
    resolve_classifier_backend,
)


class TestClassifierBackend(unittest.TestCase):
    def test_auto_prefers_timber_when_artifact_exists(self):
        probe = inspect_classifier_backends()
        if not probe.timber_available:
            self.skipTest("Timber artifact is not available on this machine")

        backend = resolve_classifier_backend()
        self.assertEqual(backend.name, "timber")

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
