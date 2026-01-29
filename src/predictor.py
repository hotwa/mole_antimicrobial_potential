"""Core antimicrobial prediction logic."""

from __future__ import annotations

import asyncio
import logging
import os
import pickle
import re
import time
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import torch
from scipy.stats.mstats import gmean
from sklearn.preprocessing import OneHotEncoder

from mole_representation import process_representation
from src.models import MoleculeInfo, MoleculeInput, StatusModel

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class AntimicrobialPredictor:
    """Predict antimicrobial activity with lazy async loading."""

    def __init__(self) -> None:
        # 动态选择 GPU，根据 CUDA_VISIBLE_DEVICES
        if torch.cuda.is_available():
            cuda_visible = os.environ.get("CUDA_VISIBLE_DEVICES", "")
            if cuda_visible:
                # 取第一个可用的 GPU
                gpu_id = cuda_visible.split(",")[0].strip()
                self.device = f"cuda:{gpu_id}"
            else:
                self.device = "cuda:0"
        else:
            self.device = "cpu"
        self.xgboost_model_path = (
            "data/03.model_evaluation/MolE-XGBoost-08.03.2024_14.20.pkl"
        )
        self.mole_model_path = "pretrained_model/model_ginconcat_btwin_100k_d8000_l0.0001"
        self.strain_categories = "data/01.prepare_training_data/maier_screening_results.tsv.gz"
        self.gram_information = "raw_data/maier_microbiome/strain_info_SF2.xlsx"

        self._model_loaded = False
        self._load_task: Optional[asyncio.Task[None]] = None
        self._loading_lock = asyncio.Lock()

        self.model = None
        self.maier_screen = None
        self.strain_ohe = None
        self._gram_dict: Optional[Dict[str, str]] = None

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
            model_path=self.xgboost_model_path,
        )

    async def _async_load(self) -> None:
        loop = asyncio.get_running_loop()

        start = time.monotonic()
        self.model = await loop.run_in_executor(None, self._load_xgb_model)
        logger.info("Loaded xgboost model in %.2fs", time.monotonic() - start)

        start = time.monotonic()
        self.maier_screen = await loop.run_in_executor(
            None,
            lambda: pd.read_csv(self.strain_categories, sep="\t", index_col=0),
        )
        logger.info("Loaded maier_screen data in %.2fs", time.monotonic() - start)

        start = time.monotonic()
        self.strain_ohe = await loop.run_in_executor(
            None, self._prep_ohe, self.maier_screen.columns
        )
        logger.info("Prepared strain OHE in %.2fs", time.monotonic() - start)

        start = time.monotonic()
        self._gram_dict = await loop.run_in_executor(None, self._load_gram_dict)
        logger.info("Loaded gram excel data in %.2fs", time.monotonic() - start)

        self._model_loaded = True

    def _load_xgb_model(self):
        logger.info("Loading XGBoost model from %s", self.xgboost_model_path)
        if not os.path.exists(self.xgboost_model_path):
            raise FileNotFoundError(f"Model file not found: {self.xgboost_model_path}")
        with open(self.xgboost_model_path, "rb") as file:
            return pickle.load(file)

    def _prep_ohe(self, categories):
        ohe = OneHotEncoder(sparse=False)
        ohe.fit(pd.DataFrame(categories))
        return pd.DataFrame(
            ohe.transform(pd.DataFrame(categories)),
            columns=categories,
            index=categories,
        )

    def _load_gram_dict(self) -> Dict[str, str]:
        maier_strains = pd.read_excel(
            self.gram_information,
            skiprows=[0, 1, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54],
            index_col="NT data base",
        )
        return maier_strains[["Gram stain"]].to_dict()["Gram stain"]

    def _get_mole_representation(self, molecules: List[MoleculeInfo]):
        smiles = [mol.smiles for mol in molecules]
        chem_ids = [mol.chem_id or f"mol{index + 1}" for index, mol in enumerate(molecules)]
        df = pd.DataFrame({"smiles": smiles, "chem_id": chem_ids})
        return process_representation(
            dataset_path=df,
            smile_column_str="smiles",
            id_column_str="chem_id",
            pretrained_dir=self.mole_model_path,
            device=self.device,
        )

    def _add_strains(self, chemfeats_df: pd.DataFrame) -> pd.DataFrame:
        chemfe = chemfeats_df.reset_index().rename(columns={"index": "chem_id"})
        chemfe["chem_id"] = chemfe["chem_id"].astype(str)

        sohe = self.strain_ohe.reset_index().rename(columns={"index": "strain_name"})
        xpred = chemfe.merge(sohe, how="cross")
        xpred["pred_id"] = xpred["chem_id"].str.cat(xpred["strain_name"], sep=":")

        xpred = xpred.set_index("pred_id")
        xpred = xpred.drop(columns=["chem_id", "strain_name"])
        return xpred

    def _gram_stain(self, label_df: pd.DataFrame) -> pd.DataFrame:
        if self._gram_dict is None:
            raise RuntimeError("Gram stain data not loaded")

        df_label = label_df.copy()
        df_label["nt_number"] = df_label["strain_name"].apply(
            lambda x: re.search(r".*?\((NT\d+)\)", x).group(1)
        )
        df_label["gram_stain"] = df_label["nt_number"].map(self._gram_dict)
        return df_label

    def _antimicrobial_potential(self, score_df: pd.DataFrame) -> pd.DataFrame:
        score_df["chem_id"] = score_df["pred_id"].str.split(":", expand=True)[0]
        score_df["strain_name"] = score_df["pred_id"].str.split(":", expand=True)[1]

        pred_df = self._gram_stain(score_df)

        apscore_total = (
            pred_df.groupby("chem_id")["1"].apply(gmean).to_frame().rename(columns={"1": "apscore_total"})
        )
        apscore_total["apscore_total"] = np.log(apscore_total["apscore_total"])

        apscore_gram = (
            pred_df.groupby(["chem_id", "gram_stain"])["1"]
            .apply(gmean)
            .unstack()
            .rename(columns={"negative": "apscore_gnegative", "positive": "apscore_gpositive"})
        )
        apscore_gram["apscore_gnegative"] = np.log(apscore_gram["apscore_gnegative"])
        apscore_gram["apscore_gpositive"] = np.log(apscore_gram["apscore_gpositive"])

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

    async def predict(self, input_data: MoleculeInput) -> List[Dict[str, Any]]:
        """Predict antimicrobial activity and return items only."""
        if not self._model_loaded:
            raise RuntimeError("Model not loaded. Call ensure_loaded() first.")

        normalized = input_data.normalize()
        molecules = normalized.molecules or []
        if not molecules:
            raise ValueError("No molecules provided for prediction")

        loop = asyncio.get_running_loop()
        try:
            mole_representation = await asyncio.wait_for(
                loop.run_in_executor(None, self._get_mole_representation, molecules),
                timeout=60,
            )
        except asyncio.TimeoutError as exc:
            raise TimeoutError("MolE representation generation timeout") from exc

        x_input = self._add_strains(mole_representation)

        try:
            y_pred = await asyncio.wait_for(
                loop.run_in_executor(None, self.model.predict_proba, x_input),
                timeout=30,
            )
        except asyncio.TimeoutError as exc:
            raise TimeoutError("Prediction timeout") from exc

        pred_df = pd.DataFrame(y_pred, columns=["0", "1"], index=x_input.index)
        pred_df["growth_inhibition"] = pred_df["1"].apply(
            lambda x: 1 if x >= normalized.app_threshold else 0
        )

        if normalized.aggregate_scores:
            pred_df = pred_df.reset_index()
            agg_df = self._antimicrobial_potential(pred_df)
            agg_df["broad_spectrum"] = agg_df["ginhib_total"].apply(
                lambda x: 1 if x >= normalized.min_nkill else 0
            )
            return agg_df.reset_index().to_dict(orient="records")

        pred_df = pred_df.drop(columns=["0"]).rename(
            columns={"1": "antimicrobial_predictive_probability"}
        )
        return pred_df.reset_index().to_dict(orient="records")
