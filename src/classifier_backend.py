"""Classifier backend adapters for antimicrobial prediction.

The prediction pipeline uses MolE to build embeddings and a classifier to turn
those embeddings into strain-level probabilities.  This module isolates the
classifier side so the application can switch between the original pickle /
XGBoost model and the compiled Timber artifact without changing the rest of
the predictor stack.
"""

from __future__ import annotations

import ctypes
import json
import logging
import os
import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Protocol

import numpy as np


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PICKLE_MODEL_PATH = REPO_ROOT / "data" / "03.model_evaluation" / "MolE-XGBoost-08.03.2024_14.20.pkl"
DEFAULT_TIMBER_MODEL_DIR = Path.home() / ".timber" / "models" / "mole-antimicrobial"
LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class BackendStatus:
    """Human-readable backend discovery / load status."""

    backend: str
    available: bool
    path: str
    detail: str = ""


@dataclass(frozen=True)
class BackendProbe:
    """Non-loading probe of the available classifier backends."""

    preference: str
    pickle_path: Path
    pickle_available: bool
    timber_model_dir: Path
    timber_library_path: Optional[Path]
    timber_available: bool
    selected_backend: str


class ClassifierBackend(Protocol):
    """Protocol shared by classifier backends."""

    name: str

    def predict_proba(self, X: np.ndarray) -> np.ndarray: ...

    def status(self) -> BackendStatus: ...


def _resolve_timber_library_path(model_dir: Path) -> Optional[Path]:
    candidates = (
        model_dir / "compiled" / "libtimber_model.so",
        model_dir / "compiled" / "timber_model.so",
        model_dir / "lib" / "libtimber_model.so",
        model_dir / "lib" / "timber_model.so",
        model_dir / "lib" / "model.so",
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def inspect_classifier_backends(
    preference: Optional[str] = None,
    pickle_path: str | Path | None = None,
    timber_model_dir: str | Path | None = None,
) -> BackendProbe:
    """Probe backend availability without loading heavy model objects."""

    backend_preference = (preference or os.environ.get("MOLE_CLASSIFIER_BACKEND", "auto")).strip().lower()
    if backend_preference not in {"auto", "timber", "pickle"}:
        raise ValueError(
            "MOLE_CLASSIFIER_BACKEND must be one of: auto, timber, pickle"
        )

    resolved_pickle_path = Path(
        pickle_path
        or os.environ.get("MOLE_PICKLE_MODEL_PATH", DEFAULT_PICKLE_MODEL_PATH)
    ).expanduser().resolve()
    resolved_timber_dir = Path(
        timber_model_dir
        or os.environ.get("MOLE_TIMBER_MODEL_DIR", DEFAULT_TIMBER_MODEL_DIR)
    ).expanduser().resolve()

    timber_lib_path = _resolve_timber_library_path(resolved_timber_dir)
    timber_available = timber_lib_path is not None and resolved_timber_dir.is_dir()
    pickle_available = resolved_pickle_path.is_file()

    if backend_preference == "timber":
        selected_backend = "timber"
    elif backend_preference == "pickle":
        selected_backend = "pickle"
    else:
        # The original MolE prediction path uses the pickle/XGBoost model.
        # Keep auto mode canonical and fast; Timber remains an explicit
        # portability option when the Python XGBoost runtime is unavailable.
        selected_backend = "pickle" if pickle_available else "timber"

    return BackendProbe(
        preference=backend_preference,
        pickle_path=resolved_pickle_path,
        pickle_available=pickle_available,
        timber_model_dir=resolved_timber_dir,
        timber_library_path=timber_lib_path,
        timber_available=timber_available,
        selected_backend=selected_backend,
    )


class PickleXGBoostBackend:
    """Original sklearn / XGBoost pickle backend."""

    name = "pickle"

    def __init__(self, model_path: str | Path | None = None):
        self.model_path = Path(
            model_path
            or os.environ.get("MOLE_PICKLE_MODEL_PATH", DEFAULT_PICKLE_MODEL_PATH)
        ).expanduser().resolve()
        if not self.model_path.is_file():
            raise FileNotFoundError(f"Pickle model file not found: {self.model_path}")
        with self.model_path.open("rb") as handle:
            self._model = pickle.load(handle)

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        X = np.asarray(X, dtype=np.float32)
        return self._model.predict_proba(X)

    def status(self) -> BackendStatus:
        return BackendStatus(
            backend=self.name,
            available=True,
            path=str(self.model_path),
            detail="pickle XGBoost classifier loaded",
        )


class TimberCompiledArtifactBackend:
    """Timber compiled shared-library backend.

    This backend uses the compiled `timber_model.so` artifact emitted by the
    Timber toolchain.  It does not import Timber runtime Python modules, so the
    main application can remain self-contained.
    """

    name = "timber"

    def __init__(self, model_dir: str | Path | None = None):
        self.model_dir = Path(
            model_dir
            or os.environ.get("MOLE_TIMBER_MODEL_DIR", DEFAULT_TIMBER_MODEL_DIR)
        ).expanduser().resolve()
        if not self.model_dir.is_dir():
            raise FileNotFoundError(f"Timber model directory not found: {self.model_dir}")

        self.model_info_path = self.model_dir / "model_info.json"
        self.library_path = _resolve_timber_library_path(self.model_dir)
        if self.library_path is None:
            raise FileNotFoundError(
                f"No Timber shared library found under {self.model_dir / 'lib'}"
            )

        self.model_info = self._load_model_info()
        self.n_features = int(self.model_info.get("n_features", 1040))
        self.n_outputs = int(self.model_info.get("n_outputs", 1))
        self._lib: Optional[ctypes.CDLL] = None
        self._ctx: Optional[ctypes.c_void_p] = None
        self._load_library()

    def _load_model_info(self) -> dict:
        if not self.model_info_path.is_file():
            return {}
        try:
            return json.loads(self.model_info_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            return {}

    def _load_library(self) -> None:
        self._lib = ctypes.CDLL(str(self.library_path))

        self._lib.timber_init.argtypes = [ctypes.POINTER(ctypes.c_void_p)]
        self._lib.timber_init.restype = ctypes.c_int

        self._lib.timber_free.argtypes = [ctypes.c_void_p]
        self._lib.timber_free.restype = None

        self._lib.timber_infer.argtypes = [
            ctypes.POINTER(ctypes.c_float),
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_float),
            ctypes.c_void_p,
        ]
        self._lib.timber_infer.restype = ctypes.c_int

        self._lib.timber_infer_single.argtypes = [
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_float),
            ctypes.c_void_p,
        ]
        self._lib.timber_infer_single.restype = ctypes.c_int

        ctx = ctypes.c_void_p()
        rc = self._lib.timber_init(ctypes.byref(ctx))
        if rc != 0:
            raise RuntimeError(f"timber_init failed with code {rc}")
        self._ctx = ctx

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        if self._lib is None or self._ctx is None:
            raise RuntimeError("Timber backend is not loaded")

        X = np.ascontiguousarray(X, dtype=np.float32)
        single_sample = False
        if X.ndim == 1:
            X = X.reshape(1, -1)
            single_sample = True

        if X.shape[1] != self.n_features:
            raise ValueError(
                f"Timber backend expected {self.n_features} features, got {X.shape[1]}"
            )

        outputs = np.zeros((X.shape[0], self.n_outputs), dtype=np.float32)
        inputs_ptr = X.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        outputs_ptr = outputs.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        rc = self._lib.timber_infer(inputs_ptr, X.shape[0], outputs_ptr, self._ctx)
        if rc != 0:
            raise RuntimeError(f"timber_infer failed with code {rc}")

        if self.n_outputs == 1:
            positive = outputs.reshape(-1)
            outputs = np.column_stack([1.0 - positive, positive]).astype(np.float32, copy=False)

        if single_sample:
            return outputs.reshape(1, -1)
        return outputs

    def status(self) -> BackendStatus:
        return BackendStatus(
            backend=self.name,
            available=True,
            path=str(self.library_path),
            detail=f"Timber artifact loaded from {self.model_dir}",
        )

    def close(self) -> None:
        if self._lib is not None and self._ctx is not None:
            self._lib.timber_free(self._ctx)
            self._ctx = None

    def __del__(self):
        self.close()


def resolve_classifier_backend(
    preference: Optional[str] = None,
    pickle_path: str | Path | None = None,
    timber_model_dir: str | Path | None = None,
) -> ClassifierBackend:
    """Resolve and instantiate the preferred classifier backend."""

    probe = inspect_classifier_backends(
        preference=preference,
        pickle_path=pickle_path,
        timber_model_dir=timber_model_dir,
    )

    if probe.selected_backend == "timber":
        if not probe.timber_available or probe.timber_library_path is None:
            if probe.preference == "timber":
                raise FileNotFoundError(
                    "Timber backend requested but no compiled artifact was found. "
                    f"Expected a shared library under {probe.timber_model_dir / 'lib'}"
                )
            if not probe.pickle_available:
                raise FileNotFoundError(
                    "No usable classifier backend found: Timber artifact missing and pickle model missing."
                )
            return PickleXGBoostBackend(probe.pickle_path)
        try:
            return TimberCompiledArtifactBackend(probe.timber_model_dir)
        except Exception as exc:
            if probe.preference == "timber":
                raise
            if not probe.pickle_available:
                raise
            LOGGER.warning(
                "Timber backend failed to load (%s); falling back to pickle backend.",
                exc,
            )
            return PickleXGBoostBackend(probe.pickle_path)

    if not probe.pickle_available:
        if probe.preference == "pickle":
            raise FileNotFoundError(
                f"Pickle backend requested but file not found: {probe.pickle_path}"
            )
        if probe.timber_available:
            try:
                return TimberCompiledArtifactBackend(probe.timber_model_dir)
            except Exception as exc:
                if probe.preference == "timber":
                    raise
                raise RuntimeError(
                    "Timber backend appears available but failed to load and no pickle backend exists"
                ) from exc
        raise FileNotFoundError(
            "No usable classifier backend found: pickle model missing and Timber artifact missing."
        )
    try:
        return PickleXGBoostBackend(probe.pickle_path)
    except Exception as exc:
        if probe.timber_available:
            if probe.preference == "timber":
                raise
            LOGGER.warning(
                "Pickle backend failed to load (%s); falling back to Timber backend.",
                exc,
            )
            return TimberCompiledArtifactBackend(probe.timber_model_dir)
        raise
