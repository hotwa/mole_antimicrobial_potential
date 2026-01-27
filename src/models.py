"""
Shared data models for MCP and FastAPI.
Avoids import-time side effects.
"""

from typing import Any, Dict, List, Optional, Sequence, Union

from pydantic import BaseModel, Field, field_validator

try:
    from rdkit import Chem

    RDKIT_AVAILABLE = True
except ImportError:  # pragma: no cover - optional dependency
    RDKIT_AVAILABLE = False


class MoleculeInfo(BaseModel):
    """Single molecule input."""

    smiles: str
    chem_id: Optional[str] = None

    @field_validator("smiles")
    @classmethod
    def validate_smiles(cls, value: str) -> str:
        if not value or not isinstance(value, str):
            raise ValueError("smiles must be a non-empty string")
        if RDKIT_AVAILABLE and Chem.MolFromSmiles(value) is None:
            raise ValueError(f"invalid SMILES: {value}")
        return value


class MoleculeInput(BaseModel):
    """
    Unified molecule input.

    Priority order:
    1) molecules
    2) smiles + chem_id
    3) smiles
    """

    molecules: Optional[List[MoleculeInfo]] = Field(
        default=None,
        description="List of {smiles, chem_id}. Preferred structured input.",
    )
    smiles: Optional[Union[str, List[str]]] = Field(
        default=None,
        description="SMILES string or list of SMILES strings.",
    )
    chem_id: Optional[Union[str, List[str]]] = Field(
        default=None,
        description="Compound id or list of ids aligned with smiles.",
    )
    aggregate_scores: bool = Field(
        default=False,
        description="True: aggregated antimicrobial potential. False: per-strain probabilities.",
    )
    app_threshold: float = Field(
        default=0.04374140128493309,
        ge=0,
        le=1,
        description="Threshold (0-1) for binarizing growth inhibition.",
    )
    min_nkill: int = Field(
        default=10,
        ge=0,
        description="Minimum inhibited strain count to label broad_spectrum=1.",
    )

    @field_validator("smiles")
    @classmethod
    def validate_smiles_input(cls, value: Optional[Union[str, List[str]]]) -> Optional[Union[str, List[str]]]:
        if value is None:
            return value
        if isinstance(value, str):
            if RDKIT_AVAILABLE and Chem.MolFromSmiles(value) is None:
                raise ValueError(f"invalid SMILES: {value}")
            return value
        if not isinstance(value, list) or len(value) == 0:
            raise ValueError("smiles list cannot be empty")
        if RDKIT_AVAILABLE:
            for smiles in value:
                if Chem.MolFromSmiles(smiles) is None:
                    raise ValueError(f"invalid SMILES in list: {smiles}")
        return value

    @field_validator("chem_id")
    @classmethod
    def validate_chem_id(
        cls, value: Optional[Union[str, List[str]]], info
    ) -> Optional[Union[str, List[str]]]:
        if value is None:
            return value
        smiles = info.data.get("smiles")
        if smiles is None:
            return value
        smiles_count = 1 if isinstance(smiles, str) else len(smiles)
        chem_id_count = 1 if isinstance(value, str) else len(value)
        if chem_id_count != smiles_count:
            raise ValueError(
                f"chem_id length ({chem_id_count}) must match smiles length ({smiles_count})"
            )
        return value

    @field_validator("molecules")
    @classmethod
    def validate_molecules(
        cls, value: Optional[Sequence[MoleculeInfo]]
    ) -> Optional[Sequence[MoleculeInfo]]:
        if value is None:
            return value
        if len(value) == 0:
            raise ValueError("molecules list cannot be empty")
        return value

    def normalize(self) -> "MoleculeInput":
        """Normalize input into a molecules list with default chem_id values."""
        if self.molecules is not None:
            normalized = []
            for index, molecule in enumerate(self.molecules, start=1):
                chem_id = molecule.chem_id or f"mol{index}"
                normalized.append(MoleculeInfo(smiles=molecule.smiles, chem_id=chem_id))
            return MoleculeInput(
                molecules=normalized,
                aggregate_scores=self.aggregate_scores,
                app_threshold=self.app_threshold,
                min_nkill=self.min_nkill,
            )

        if self.smiles is None:
            raise ValueError("Either molecules or smiles must be provided")

        smiles_list = [self.smiles] if isinstance(self.smiles, str) else list(self.smiles)
        if self.chem_id is None:
            chem_id_list = [f"mol{i + 1}" for i in range(len(smiles_list))]
        else:
            chem_id_list = (
                [self.chem_id] if isinstance(self.chem_id, str) else list(self.chem_id)
            )
        if len(chem_id_list) != len(smiles_list):
            raise ValueError(
                f"chem_id length ({len(chem_id_list)}) must match smiles length ({len(smiles_list)})"
            )

        molecules = [
            MoleculeInfo(smiles=smiles, chem_id=chem_id)
            for smiles, chem_id in zip(smiles_list, chem_id_list)
        ]
        return MoleculeInput(
            molecules=molecules,
            aggregate_scores=self.aggregate_scores,
            app_threshold=self.app_threshold,
            min_nkill=self.min_nkill,
        )


class StatusModel(BaseModel):
    """Service status response."""

    service: str = "mole-antimicrobial-prediction"
    loaded: bool = False
    device: Optional[str] = None
    version: str = "1.0.0"
    model_path: str = "data/03.model_evaluation/MolE-XGBoost-08.03.2024_14.20.pkl"
