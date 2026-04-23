"""Core antimicrobial prediction logic."""

from __future__ import annotations

import asyncio
import copy
import logging
import os
import re
import time
import threading
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import torch
from scipy.stats.mstats import gmean
from sklearn.preprocessing import OneHotEncoder

from src.mole_representation import load_pretrained_model, process_representation
from src.classifier_backend import inspect_classifier_backends, resolve_classifier_backend
from src.models import MoleculeInfo, MoleculeInput, StatusModel
from workflow.dataset.dataset_representation import iter_batch_representation

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[1]
NT_NUMBER_PATTERN = re.compile(r".*?\((NT\d+)\)")
GRAPH_BUILD_PROFILE_SUM_FIELDS = (
    "graph_items",
    "rdkit_parse_seconds",
    "add_hs_seconds",
    "atom_feature_seconds",
    "bond_feature_seconds",
    "graph_total_seconds",
    "dataloader_setup_seconds",
    "dataloader_iter_seconds",
    "model_forward_seconds",
)
DEFAULT_CLASSIFIER_WORKERS: int | str = "auto"
DEFAULT_CLASSIFIER_INFLIGHT_BATCHES: int | str = "auto"


def _format_strain_feature_name(strain_name: str) -> str:
    """Match the legacy one-hot column names stored in the pickled XGBoost model."""
    return str((strain_name,))


class AntimicrobialPredictor:
    """Predict antimicrobial activity with lazy async loading."""

    def __init__(self) -> None:
        # Respect CUDA_VISIBLE_DEVICES remapping. In a one-process-per-GPU
        # launch pattern the selected device inside the masked process is
        # always cuda:0, regardless of the original physical GPU index.
        if torch.cuda.is_available():
            cuda_visible = os.environ.get("CUDA_VISIBLE_DEVICES", "")
            if cuda_visible:
                self.device = "cuda:0"
            else:
                self.device = "cuda:0"
        else:
            self.device = "cpu"
        self.classifier_backend_probe = inspect_classifier_backends()
        self.classifier_backend_preference = self.classifier_backend_probe.preference
        self.xgboost_model_path = str(self.classifier_backend_probe.pickle_path)
        self.timber_model_dir = str(self.classifier_backend_probe.timber_model_dir)
        default_mole_model = REPO_ROOT / "pretrained_model" / "model_ginconcat_btwin_100k_d8000_l0.0001"
        self.mole_model_path = os.environ.get("MOLE_MOLE_MODEL_PATH", str(default_mole_model))
        self.strain_categories = str(REPO_ROOT / "data" / "01.prepare_training_data" / "maier_screening_results.tsv.gz")
        self.strain_index_path = str(REPO_ROOT / "workflows" / "reinvent4" / "inputs" / "strain_index.tsv")
        self.gram_information = str(REPO_ROOT / "raw_data" / "maier_microbiome" / "strain_info_SF2.xlsx")

        self._model_loaded = False
        self._load_task: Optional[asyncio.Task[None]] = None
        self._loading_lock = asyncio.Lock()

        self.classifier_backend = None
        self.model = None
        self.maier_screen = None
        self.strain_ohe = None
        self._gram_dict: Optional[Dict[str, str]] = None
        self.mole_model = None
        self._mole_model_lock = threading.Lock()
        self._last_profile: Optional[Dict[str, Any]] = None

    async def ensure_loaded(self) -> None:
        """Load models and metadata once under concurrency."""
        if self._model_loaded:
            return
        async with self._loading_lock:
            if self._model_loaded:
                return
            if self._load_task is None:
                self._load_task = asyncio.create_task(self._async_load())
        try:
            await self._load_task
        except Exception:
            async with self._loading_lock:
                if self._load_task and self._load_task.done():
                    self._load_task = None
            raise

    def get_status(self) -> StatusModel:
        """Return service status without triggering load."""
        return StatusModel(
            loaded=self._model_loaded,
            device=self.device,
            model_path=self.model.status().path if self._model_loaded and self.model else self.xgboost_model_path,
            classifier_backend=self.model.status().backend if self._model_loaded and self.model else None,
            classifier_backend_path=self.model.status().path if self._model_loaded and self.model else None,
            classifier_backend_preference=self.classifier_backend_preference,
        )

    async def _async_load(self) -> None:
        loop = asyncio.get_running_loop()

        start = time.monotonic()
        self.classifier_backend = await loop.run_in_executor(
            None,
            resolve_classifier_backend,
            self.classifier_backend_preference,
            self.classifier_backend_probe.pickle_path,
            self.classifier_backend_probe.timber_model_dir,
        )
        self.model = self.classifier_backend
        logger.info(
            "Loaded classifier backend (%s) in %.2fs",
            self.classifier_backend.status().backend,
            time.monotonic() - start,
        )

        start = time.monotonic()
        self.maier_screen, selected_strain_source = await loop.run_in_executor(
            None, self._load_strain_panel
        )
        self.strain_categories = selected_strain_source
        logger.info("Loaded maier_screen data in %.2fs", time.monotonic() - start)

        start = time.monotonic()
        self.strain_ohe = await loop.run_in_executor(
            None, self._prep_ohe, self.maier_screen.columns
        )
        logger.info("Prepared strain OHE in %.2fs", time.monotonic() - start)

        start = time.monotonic()
        self._gram_dict = await loop.run_in_executor(None, self._load_gram_dict)
        logger.info("Loaded gram excel data in %.2fs", time.monotonic() - start)

        start = time.monotonic()
        await loop.run_in_executor(None, self._load_mole_model)
        logger.info("Loaded MolE representation model in %.2fs", time.monotonic() - start)

        self._model_loaded = True

    def _prep_ohe(self, categories):
        try:
            ohe = OneHotEncoder(sparse_output=False)
        except TypeError:
            ohe = OneHotEncoder(sparse=False)
        ohe.fit(pd.DataFrame(categories))
        feature_columns = [_format_strain_feature_name(category) for category in categories]
        return pd.DataFrame(
            ohe.transform(pd.DataFrame(categories)),
            columns=feature_columns,
            index=categories,
        )

    def _load_strain_panel(self):
        primary_path = os.path.expanduser(self.strain_categories)
        if os.path.isfile(primary_path):
            return pd.read_csv(primary_path, sep="\t", index_col=0), primary_path

        fallback_path = os.path.expanduser(self.strain_index_path)
        if os.path.isfile(fallback_path):
            strain_index = pd.read_csv(fallback_path, sep="\t")
            if "strain_name" not in strain_index.columns:
                raise ValueError(
                    f"Fallback strain index file missing 'strain_name' column: {fallback_path}"
                )
            columns = [str(strain) for strain in strain_index["strain_name"].tolist()]
            return pd.DataFrame(columns=columns), fallback_path

        raise FileNotFoundError(
            "No strain panel file found. Expected either "
            f"{primary_path} or fallback {fallback_path}"
        )

    def _load_gram_dict(self) -> Dict[str, str]:
        maier_strains = pd.read_excel(
            self.gram_information,
            skiprows=[0, 1, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54],
            index_col="NT data base",
        )
        return maier_strains[["Gram stain"]].to_dict()["Gram stain"]

    def _load_mole_model(self):
        if self.mole_model is not None:
            return self.mole_model

        with self._mole_model_lock:
            if self.mole_model is not None:
                return self.mole_model

            self.mole_model = load_pretrained_model(
                pretrained_model_dir=self.mole_model_path,
                device=self.device,
            )
            return self.mole_model

    def _get_mole_representation(
        self,
        molecules: List[MoleculeInfo],
        num_graph_workers: int | str = "auto",
        graph_batch_size: int = 1024,
        prefetch_batches: int = 2,
        enable_profiling: bool = False,
        deterministic_representation: bool = False,
    ):
        smiles = [mol.smiles for mol in molecules]
        chem_ids = [mol.chem_id or f"mol{index + 1}" for index, mol in enumerate(molecules)]
        df = pd.DataFrame({"smiles": smiles, "chem_id": chem_ids})
        return process_representation(
            dataset_path=df,
            smile_column_str="smiles",
            id_column_str="chem_id",
            pretrained_dir=self.mole_model_path,
            device=self.device,
            model=self._load_mole_model(),
            num_graph_workers=num_graph_workers,
            graph_batch_size=graph_batch_size,
            prefetch_batches=prefetch_batches,
            validate_smiles=False,
            enable_profiling=enable_profiling,
            deterministic_representation=deterministic_representation,
        )

    def _iter_mole_representation_batches(
        self,
        molecules: List[MoleculeInfo],
        num_graph_workers: int | str = "auto",
        graph_batch_size: int = 1024,
        prefetch_batches: int = 2,
        enable_profiling: bool = False,
        deterministic_representation: bool = False,
    ):
        smiles = [mol.smiles for mol in molecules]
        chem_ids = [mol.chem_id or f"mol{index + 1}" for index, mol in enumerate(molecules)]
        df = pd.DataFrame({"smiles": smiles, "chem_id": chem_ids})
        return iter_batch_representation(
            smile_df=df,
            dl_model=self._load_mole_model(),
            column_str="smiles",
            id_str="chem_id",
            batch_size=graph_batch_size,
            device=self.device,
            num_graph_workers=num_graph_workers,
            graph_batch_size=graph_batch_size,
            prefetch_batches=prefetch_batches,
            enable_profiling=enable_profiling,
            deterministic_representation=deterministic_representation,
        )

    def last_profile_snapshot(self) -> Optional[Dict[str, Any]]:
        if self._last_profile is None:
            return None
        return copy.deepcopy(self._last_profile)

    def _add_strains(self, chemfeats_df: pd.DataFrame) -> pd.DataFrame:
        if self.strain_ohe is None:
            raise RuntimeError("Strain OHE data not loaded")

        chem_ids = chemfeats_df.index.astype(str).to_numpy(copy=False)
        strain_names = self.strain_ohe.index.astype(str).to_numpy(copy=False)
        xpred_values = self._strain_feature_array(chemfeats_df)

        pred_ids = np.char.add(
            np.char.add(np.repeat(chem_ids, len(strain_names)), ":"),
            np.tile(strain_names, len(chem_ids)),
        )

        columns = [str(column) for column in chemfeats_df.columns] + [
            str(column) for column in self.strain_ohe.columns
        ]
        return pd.DataFrame(
            xpred_values,
            index=pd.Index(pred_ids, name="pred_id"),
            columns=columns,
        )

    def _strain_feature_array(self, chemfeats_df: pd.DataFrame) -> np.ndarray:
        if self.strain_ohe is None:
            raise RuntimeError("Strain OHE data not loaded")

        chem_values = chemfeats_df.to_numpy(dtype=np.float32, copy=False)
        strain_values = self.strain_ohe.to_numpy(dtype=np.float32, copy=False)
        strain_count = len(self.strain_ohe.index)
        molecule_count = len(chemfeats_df.index)

        expanded_chem = np.repeat(chem_values, strain_count, axis=0)
        expanded_strains = np.tile(strain_values, (molecule_count, 1))
        return np.concatenate([expanded_chem, expanded_strains], axis=1)

    def _gram_stain(self, label_df: pd.DataFrame) -> pd.DataFrame:
        if self._gram_dict is None:
            raise RuntimeError("Gram stain data not loaded")

        df_label = label_df.copy()
        df_label["nt_number"] = df_label["strain_name"].apply(
            lambda x: NT_NUMBER_PATTERN.search(x).group(1)
        )
        df_label["gram_stain"] = df_label["nt_number"].map(self._gram_dict)
        return df_label

    def _resolve_strain_grams(self, strain_names: np.ndarray) -> np.ndarray:
        if self._gram_dict is None:
            raise RuntimeError("Gram stain data not loaded")

        gram_by_strain: dict[str, str | None] = {}
        resolved = np.empty(len(strain_names), dtype=object)
        for index, strain_name in enumerate(strain_names):
            strain_name_key = str(strain_name)
            gram = gram_by_strain.get(strain_name_key)
            if gram is None and strain_name_key not in gram_by_strain:
                nt_number = NT_NUMBER_PATTERN.search(strain_name_key).group(1)
                gram = self._gram_dict.get(nt_number)
                gram_by_strain[strain_name_key] = gram
            resolved[index] = gram_by_strain[strain_name_key]
        return resolved

    def _aggregate_scores_from_matrix(
        self,
        probability_matrix: np.ndarray,
        growth_inhibition_matrix: np.ndarray,
        chem_ids: np.ndarray,
        strain_names: np.ndarray,
    ) -> pd.DataFrame:
        probabilities = np.asarray(probability_matrix)
        growth_inhibition = np.asarray(growth_inhibition_matrix)
        chem_ids_array = np.asarray(chem_ids, dtype=object)
        strain_names_array = np.asarray(strain_names, dtype=object)

        if probabilities.ndim != 2:
            raise ValueError("probability_matrix must be a 2D array")
        if growth_inhibition.shape != probabilities.shape:
            raise ValueError("growth_inhibition_matrix must match probability_matrix shape")
        if probabilities.shape[0] != len(chem_ids_array):
            raise ValueError("chem_ids length must match probability_matrix row count")
        if probabilities.shape[1] != len(strain_names_array):
            raise ValueError("strain_names length must match probability_matrix column count")

        gram_by_strain = self._resolve_strain_grams(strain_names_array)
        negative_mask = gram_by_strain == "negative"
        positive_mask = gram_by_strain == "positive"
        if not negative_mask.any():
            raise KeyError("apscore_gnegative")
        if not positive_mask.any():
            raise KeyError("apscore_gpositive")

        order = np.argsort(chem_ids_array.astype(str), kind="stable")
        sorted_chem_ids = chem_ids_array[order].astype(str)
        sorted_probabilities = probabilities[order]
        sorted_growth_inhibition = growth_inhibition[order]
        inhibition_dtype = growth_inhibition.dtype

        def _log2_gmean_rows(values: np.ndarray) -> np.ndarray:
            output_dtype = values.dtype if np.issubdtype(values.dtype, np.floating) else np.float64
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.fromiter(
                    (float(np.log2(gmean(row))) for row in values),
                    dtype=output_dtype,
                    count=values.shape[0],
                )

        with np.errstate(divide="ignore", invalid="ignore"):
            apscore_total = _log2_gmean_rows(sorted_probabilities)
            apscore_gnegative = _log2_gmean_rows(sorted_probabilities[:, negative_mask])
            apscore_gpositive = _log2_gmean_rows(sorted_probabilities[:, positive_mask])

        result = pd.DataFrame(
            {
                "apscore_total": apscore_total,
                "apscore_gnegative": apscore_gnegative,
                "apscore_gpositive": apscore_gpositive,
                "ginhib_total": sorted_growth_inhibition.sum(axis=1, dtype=np.int64).astype(inhibition_dtype, copy=False),
                "ginhib_gnegative": sorted_growth_inhibition[:, negative_mask].sum(axis=1, dtype=np.int64).astype(inhibition_dtype, copy=False),
                "ginhib_gpositive": sorted_growth_inhibition[:, positive_mask].sum(axis=1, dtype=np.int64).astype(inhibition_dtype, copy=False),
            },
            index=pd.Index(sorted_chem_ids, name="chem_id"),
        )
        return result

    def _antimicrobial_potential(self, score_df: pd.DataFrame) -> pd.DataFrame:
        split = score_df["pred_id"].astype(str).str.rsplit(":", n=1, expand=True)
        score_df["chem_id"] = split[0]
        score_df["strain_name"] = split[1]

        pred_df = self._gram_stain(score_df)

        apscore_total = (
            pred_df.groupby("chem_id")["1"].apply(gmean).to_frame().rename(columns={"1": "apscore_total"})
        )
        apscore_total["apscore_total"] = np.log2(apscore_total["apscore_total"])

        apscore_gram = (
            pred_df.groupby(["chem_id", "gram_stain"])["1"]
            .apply(gmean)
            .unstack()
            .rename(columns={"negative": "apscore_gnegative", "positive": "apscore_gpositive"})
        )
        apscore_gram["apscore_gnegative"] = np.log2(apscore_gram["apscore_gnegative"])
        apscore_gram["apscore_gpositive"] = np.log2(apscore_gram["apscore_gpositive"])

        inhibted_total = (
            pred_df.groupby("chem_id")["growth_inhibition"]
            .sum()
            .to_frame()
            .rename(columns={"growth_inhibition": "ginhib_total"})
        )

        inhibted_gram = (
            pred_df.groupby(["chem_id", "gram_stain"])["growth_inhibition"]
            .sum()
            .unstack()
            .rename(columns={"negative": "ginhib_gnegative", "positive": "ginhib_gpositive"})
        )

        return apscore_total.join(apscore_gram).join(inhibted_total).join(inhibted_gram)

    def _new_graph_build_profile(self) -> Dict[str, Any]:
        return {
            "graph_items": 0,
            "rdkit_parse_seconds": 0.0,
            "add_hs_seconds": 0.0,
            "atom_feature_seconds": 0.0,
            "bond_feature_seconds": 0.0,
            "graph_total_seconds": 0.0,
            "dataloader_setup_seconds": 0.0,
            "dataloader_iter_seconds": 0.0,
            "model_forward_seconds": 0.0,
            "graph_batch_size": None,
            "graph_workers": None,
        }

    def _merge_graph_build_profile(
        self,
        total_profile: Optional[Dict[str, Any]],
        batch_profile: Optional[Dict[str, Any]],
    ) -> Optional[Dict[str, Any]]:
        if batch_profile is None:
            return total_profile

        if total_profile is None:
            total_profile = self._new_graph_build_profile()

        for field in GRAPH_BUILD_PROFILE_SUM_FIELDS:
            total_profile[field] += batch_profile.get(field, 0.0)

        if "graph_batch_size" in batch_profile:
            total_profile["graph_batch_size"] = batch_profile["graph_batch_size"]
        if "graph_workers" in batch_profile:
            total_profile["graph_workers"] = batch_profile["graph_workers"]

        return total_profile

    def _predict_aggregate_streaming_sync(
        self,
        normalized: MoleculeInput,
        num_graph_workers: int | str,
        graph_batch_size: int,
        prefetch_batches: int,
        enable_profiling: bool,
        deterministic_representation: bool,
        classifier_workers: int | str = DEFAULT_CLASSIFIER_WORKERS,
        classifier_inflight_batches: int | str = DEFAULT_CLASSIFIER_INFLIGHT_BATCHES,
    ) -> tuple[List[Dict[str, Any]], Optional[Dict[str, Any]]]:
        molecules = normalized.molecules or []
        if not molecules:
            raise ValueError("No molecules provided for prediction")

        strain_names = self.strain_ohe.index.astype(str).to_numpy(copy=False)
        resolved_classifier_workers = self._resolve_classifier_workers(classifier_workers)
        resolved_classifier_inflight_batches = self._resolve_classifier_inflight_batches(
            classifier_inflight_batches,
            resolved_classifier_workers,
        )
        iterator = self._iter_mole_representation_batches(
            molecules,
            num_graph_workers=num_graph_workers,
            graph_batch_size=graph_batch_size,
            prefetch_batches=prefetch_batches,
            enable_profiling=enable_profiling,
            deterministic_representation=deterministic_representation,
        )

        representation_seconds = 0.0
        expand_seconds = 0.0
        xgboost_seconds = 0.0
        prediction_frame_seconds = 0.0
        growth_inhibition_seconds = 0.0
        aggregate_scores_seconds = 0.0
        result_records_seconds = 0.0
        classifier_stage_seconds = 0.0
        first_result_latency_seconds: Optional[float] = None
        graph_build_profile: Optional[Dict[str, Any]] = None
        streaming_batches_by_index: dict[int, Dict[str, Any]] = {}
        aggregate_batches_by_index: dict[int, pd.DataFrame] = {}
        pending_batches: dict[Future[Dict[str, Any]], Dict[str, Any]] = {}

        pipeline_start = time.perf_counter()
        with ThreadPoolExecutor(
            max_workers=resolved_classifier_workers,
            thread_name_prefix="predictor-classifier",
        ) as classifier_executor:
            while True:
                if len(pending_batches) >= resolved_classifier_inflight_batches:
                    (
                        first_result_latency_seconds,
                        expand_seconds,
                        xgboost_seconds,
                        growth_inhibition_seconds,
                        aggregate_scores_seconds,
                        classifier_stage_seconds,
                    ) = self._collect_classifier_batch_results(
                        pending_batches=pending_batches,
                        aggregate_batches_by_index=aggregate_batches_by_index,
                        streaming_batches_by_index=streaming_batches_by_index,
                        wait_for_result=True,
                        pipeline_start=pipeline_start,
                        enable_profiling=enable_profiling,
                        first_result_latency_seconds=first_result_latency_seconds,
                        strain_expand_seconds=expand_seconds,
                        xgboost_seconds=xgboost_seconds,
                        growth_inhibition_seconds=growth_inhibition_seconds,
                        aggregate_scores_seconds=aggregate_scores_seconds,
                        classifier_stage_seconds=classifier_stage_seconds,
                    )

                representation_batch_start = time.perf_counter()
                try:
                    batch = next(iterator)
                except StopIteration:
                    break
                representation_batch_seconds = time.perf_counter() - representation_batch_start
                representation_seconds += representation_batch_seconds

                if enable_profiling:
                    graph_build_profile = self._merge_graph_build_profile(
                        graph_build_profile,
                        batch.get("profiling"),
                    )

                future = classifier_executor.submit(
                    self._classify_aggregate_batch_sync,
                    batch["batch_index"],
                    batch["embedding_batch"],
                    strain_names,
                    normalized.app_threshold,
                    normalized.min_nkill,
                )
                pending_batches[future] = {
                    "batch_index": batch["batch_index"],
                    "chem_ids": list(batch["chem_ids"]),
                    "representation_batch_production_seconds": representation_batch_seconds,
                }
                (
                    first_result_latency_seconds,
                    expand_seconds,
                    xgboost_seconds,
                    growth_inhibition_seconds,
                    aggregate_scores_seconds,
                    classifier_stage_seconds,
                ) = self._collect_classifier_batch_results(
                    pending_batches=pending_batches,
                    aggregate_batches_by_index=aggregate_batches_by_index,
                    streaming_batches_by_index=streaming_batches_by_index,
                    wait_for_result=False,
                    pipeline_start=pipeline_start,
                    enable_profiling=enable_profiling,
                    first_result_latency_seconds=first_result_latency_seconds,
                    strain_expand_seconds=expand_seconds,
                    xgboost_seconds=xgboost_seconds,
                    growth_inhibition_seconds=growth_inhibition_seconds,
                    aggregate_scores_seconds=aggregate_scores_seconds,
                    classifier_stage_seconds=classifier_stage_seconds,
                )

            while pending_batches:
                (
                    first_result_latency_seconds,
                    expand_seconds,
                    xgboost_seconds,
                    growth_inhibition_seconds,
                    aggregate_scores_seconds,
                    classifier_stage_seconds,
                ) = self._collect_classifier_batch_results(
                    pending_batches=pending_batches,
                    aggregate_batches_by_index=aggregate_batches_by_index,
                    streaming_batches_by_index=streaming_batches_by_index,
                    wait_for_result=True,
                    pipeline_start=pipeline_start,
                    enable_profiling=enable_profiling,
                    first_result_latency_seconds=first_result_latency_seconds,
                    strain_expand_seconds=expand_seconds,
                    xgboost_seconds=xgboost_seconds,
                    growth_inhibition_seconds=growth_inhibition_seconds,
                    aggregate_scores_seconds=aggregate_scores_seconds,
                    classifier_stage_seconds=classifier_stage_seconds,
                )

        if not aggregate_batches_by_index:
            raise ValueError("No molecules provided for prediction")

        ordered_batch_indices = sorted(aggregate_batches_by_index)
        agg_df = pd.concat([aggregate_batches_by_index[index] for index in ordered_batch_indices])
        order = np.argsort(agg_df.index.astype(str).to_numpy(copy=False), kind="stable")
        agg_df = agg_df.take(order)

        result_records_start = time.perf_counter()
        records = agg_df.reset_index().to_dict(orient="records")
        result_records_seconds += time.perf_counter() - result_records_start

        profile = None
        if enable_profiling:
            profile = {
                "representation_seconds": representation_seconds,
                "strain_expand_seconds": expand_seconds,
                "xgboost_seconds": xgboost_seconds,
                "prediction_frame_seconds": prediction_frame_seconds,
                "growth_inhibition_seconds": growth_inhibition_seconds,
                "aggregate_scores_seconds": aggregate_scores_seconds,
                "classifier_stage_seconds": classifier_stage_seconds,
                "classifier_workers": resolved_classifier_workers,
                "classifier_inflight_batches": resolved_classifier_inflight_batches,
                "result_records_seconds": result_records_seconds,
                "graph_build": graph_build_profile,
                "first_result_latency_seconds": first_result_latency_seconds or 0.0,
                "streaming_batches": [
                    streaming_batches_by_index[index] for index in sorted(streaming_batches_by_index)
                ],
            }

        return records, profile

    def _resolve_classifier_workers(self, classifier_workers: int | str) -> int:
        if classifier_workers == "auto":
            return 1

        resolved = int(classifier_workers)
        if resolved < 1:
            raise ValueError("classifier_workers must be >= 1 or 'auto'")
        return resolved

    def _resolve_classifier_inflight_batches(
        self,
        classifier_inflight_batches: int | str,
        classifier_workers: int,
    ) -> int:
        if classifier_inflight_batches == "auto":
            return max(2, classifier_workers + 1)

        resolved = int(classifier_inflight_batches)
        if resolved < 1:
            raise ValueError("classifier_inflight_batches must be >= 1 or 'auto'")
        return resolved

    def _classify_aggregate_batch_sync(
        self,
        batch_index: int,
        representation_df: pd.DataFrame,
        strain_names: np.ndarray,
        app_threshold: float,
        min_nkill: int,
    ) -> Dict[str, Any]:
        classifier_stage_start = time.perf_counter()

        expand_start = time.perf_counter()
        x_input_values = self._strain_feature_array(representation_df)
        strain_expand_seconds = time.perf_counter() - expand_start

        classifier_batch_start = time.perf_counter()
        y_pred = self.model.predict_proba(x_input_values)
        classifier_batch_seconds = time.perf_counter() - classifier_batch_start

        y_pred_array = np.asarray(y_pred)
        probability_values = y_pred_array[:, 1]

        growth_inhibition_start = time.perf_counter()
        chem_ids = representation_df.index.astype(str).to_numpy(copy=False)
        probability_matrix = probability_values.reshape(len(chem_ids), len(strain_names))
        growth_inhibition_matrix = (probability_matrix >= app_threshold).astype(np.int64, copy=False)
        growth_inhibition_seconds = time.perf_counter() - growth_inhibition_start

        aggregate_scores_start = time.perf_counter()
        batch_agg_df = self._aggregate_scores_from_matrix(
            probability_matrix=probability_matrix,
            growth_inhibition_matrix=growth_inhibition_matrix,
            chem_ids=chem_ids,
            strain_names=strain_names,
        )
        batch_agg_df["broad_spectrum"] = (
            batch_agg_df["ginhib_total"].to_numpy(copy=False) >= min_nkill
        ).astype(np.int64, copy=False)
        aggregate_batch_seconds = time.perf_counter() - aggregate_scores_start
        classifier_stage_total_seconds = time.perf_counter() - classifier_stage_start

        return {
            "batch_index": batch_index,
            "agg_df": batch_agg_df,
            "profiling": {
                "strain_expand_seconds": strain_expand_seconds,
                "xgboost_seconds": classifier_batch_seconds,
                "prediction_frame_seconds": 0.0,
                "growth_inhibition_seconds": growth_inhibition_seconds,
                "aggregate_scores_seconds": aggregate_batch_seconds,
                "classifier_stage_seconds": classifier_stage_total_seconds,
                "classifier_batch_inference_seconds": classifier_batch_seconds,
                "classifier_batch_stage_seconds": classifier_stage_total_seconds,
                "aggregate_accumulate_seconds": aggregate_batch_seconds,
            },
        }

    def _collect_classifier_batch_results(
        self,
        pending_batches: Dict[Future[Dict[str, Any]], Dict[str, Any]],
        aggregate_batches_by_index: Dict[int, pd.DataFrame],
        streaming_batches_by_index: Dict[int, Dict[str, Any]],
        wait_for_result: bool,
        pipeline_start: float,
        enable_profiling: bool,
        first_result_latency_seconds: Optional[float],
        strain_expand_seconds: float,
        xgboost_seconds: float,
        growth_inhibition_seconds: float,
        aggregate_scores_seconds: float,
        classifier_stage_seconds: float,
    ) -> tuple[Optional[float], float, float, float, float, float]:
        if not pending_batches:
            return (
                first_result_latency_seconds,
                strain_expand_seconds,
                xgboost_seconds,
                growth_inhibition_seconds,
                aggregate_scores_seconds,
                classifier_stage_seconds,
            )

        if wait_for_result:
            done, _ = wait(tuple(pending_batches), return_when=FIRST_COMPLETED)
        else:
            done = {future for future in pending_batches if future.done()}
            if not done:
                return (
                    first_result_latency_seconds,
                    strain_expand_seconds,
                    xgboost_seconds,
                    growth_inhibition_seconds,
                    aggregate_scores_seconds,
                    classifier_stage_seconds,
                )

        for future in done:
            metadata = pending_batches.pop(future)
            batch_result = future.result()
            batch_index = int(batch_result["batch_index"])
            aggregate_batches_by_index[batch_index] = batch_result["agg_df"]

            batch_profile = batch_result["profiling"]
            strain_expand_seconds += float(batch_profile["strain_expand_seconds"])
            xgboost_seconds += float(batch_profile["xgboost_seconds"])
            growth_inhibition_seconds += float(batch_profile["growth_inhibition_seconds"])
            aggregate_scores_seconds += float(batch_profile["aggregate_scores_seconds"])
            classifier_stage_seconds += float(batch_profile["classifier_stage_seconds"])

            if first_result_latency_seconds is None:
                first_result_latency_seconds = time.perf_counter() - pipeline_start

            if enable_profiling:
                streaming_batches_by_index[batch_index] = {
                    "batch_index": batch_index,
                    "chem_ids": list(metadata["chem_ids"]),
                    "representation_batch_production_seconds": metadata[
                        "representation_batch_production_seconds"
                    ],
                    "classifier_batch_inference_seconds": batch_profile[
                        "classifier_batch_inference_seconds"
                    ],
                    "classifier_batch_stage_seconds": batch_profile[
                        "classifier_batch_stage_seconds"
                    ],
                    "aggregate_accumulate_seconds": batch_profile[
                        "aggregate_accumulate_seconds"
                    ],
                }

        return (
            first_result_latency_seconds,
            strain_expand_seconds,
            xgboost_seconds,
            growth_inhibition_seconds,
            aggregate_scores_seconds,
            classifier_stage_seconds,
        )

    async def predict(
        self,
        input_data: MoleculeInput,
        num_graph_workers: int | str = "auto",
        graph_batch_size: int = 1024,
        prefetch_batches: int = 2,
        already_normalized: bool = False,
        enable_profiling: bool = False,
        deterministic_representation: bool = False,
        classifier_workers: int | str = DEFAULT_CLASSIFIER_WORKERS,
        classifier_inflight_batches: int | str = DEFAULT_CLASSIFIER_INFLIGHT_BATCHES,
    ) -> List[Dict[str, Any]]:
        """Predict antimicrobial activity and return items only."""
        if not self._model_loaded:
            raise RuntimeError("Model not loaded. Call ensure_loaded() first.")

        normalized = input_data if already_normalized else input_data.normalize()
        molecules = normalized.molecules or []
        if not molecules:
            raise ValueError("No molecules provided for prediction")

        self._last_profile = None
        loop = asyncio.get_running_loop()
        if normalized.aggregate_scores:
            records, profile = await loop.run_in_executor(
                None,
                self._predict_aggregate_streaming_sync,
                normalized,
                num_graph_workers,
                graph_batch_size,
                prefetch_batches,
                enable_profiling,
                deterministic_representation,
                classifier_workers,
                classifier_inflight_batches,
            )

            if enable_profiling:
                self._last_profile = profile
            return records

        try:
            representation_start = time.perf_counter()
            mole_representation = await asyncio.wait_for(
                loop.run_in_executor(
                    None,
                    self._get_mole_representation,
                    molecules,
                    num_graph_workers,
                    graph_batch_size,
                    prefetch_batches,
                    enable_profiling,
                    deterministic_representation,
                ),
                timeout=60,
            )
            representation_seconds = time.perf_counter() - representation_start
        except asyncio.TimeoutError as exc:
            raise TimeoutError("MolE representation generation timeout") from exc

        expand_start = time.perf_counter()
        x_input = self._add_strains(mole_representation)
        expand_seconds = time.perf_counter() - expand_start

        try:
            xgboost_start = time.perf_counter()
            x_input_values = x_input.to_numpy(dtype=np.float32, copy=False)
            y_pred = await asyncio.wait_for(
                loop.run_in_executor(None, self.model.predict_proba, x_input_values),
                timeout=30,
            )
            xgboost_seconds = time.perf_counter() - xgboost_start
        except asyncio.TimeoutError as exc:
            raise TimeoutError("Prediction timeout") from exc

        y_pred_array = np.asarray(y_pred)
        probability_values = y_pred_array[:, 1]
        chem_ids = mole_representation.index.astype(str).to_numpy(copy=False)
        aggregate_scores_seconds = 0.0
        result_records_seconds = 0.0

        prediction_frame_start = time.perf_counter()
        pred_df = pd.DataFrame(y_pred_array, columns=["0", "1"], index=x_input.index)
        prediction_frame_seconds = time.perf_counter() - prediction_frame_start

        growth_inhibition_start = time.perf_counter()
        growth_inhibition = (
            probability_values >= normalized.app_threshold
        ).astype(np.int64, copy=False)
        pred_df["growth_inhibition"] = growth_inhibition
        growth_inhibition_seconds = time.perf_counter() - growth_inhibition_start

        pred_df = pred_df.drop(columns=["0"]).rename(
            columns={"1": "antimicrobial_predictive_probability"}
        )
        result_records_start = time.perf_counter()
        records = pred_df.reset_index().to_dict(orient="records")
        result_records_seconds = time.perf_counter() - result_records_start

        if enable_profiling:
            self._last_profile = {
                "representation_seconds": representation_seconds,
                "strain_expand_seconds": expand_seconds,
                "xgboost_seconds": xgboost_seconds,
                "prediction_frame_seconds": prediction_frame_seconds,
                "growth_inhibition_seconds": growth_inhibition_seconds,
                "aggregate_scores_seconds": aggregate_scores_seconds,
                "result_records_seconds": result_records_seconds,
                "graph_build": copy.deepcopy(mole_representation.attrs.get("profiling")),
            }

        return records
