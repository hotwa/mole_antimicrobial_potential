"""
Shared data models for MCP and FastAPI.
Avoids import-time side effects.
"""

from typing import Any, Dict, List, Literal, Optional, Sequence, Union

from pydantic import BaseModel, ConfigDict, Field, field_validator

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

class SiteRewardConfig(BaseModel):
    """Optional structural auxiliary reward for multi-site scaffold decoration."""

    model_config = ConfigDict(populate_by_name=True)

    enabled: bool = Field(
        default=False,
        description="Enable experimental multi-site structural reward shaping.",
    )
    scaffold_smiles: Optional[str] = Field(
        default=None,
        description="Fixed scaffold used for R-group decomposition.",
    )
    range_min: int = Field(
        default=4,
        ge=0,
        description="Preferred lower bound for per-site heavy atom counts.",
    )
    range_max: int = Field(
        default=27,
        ge=0,
        description="Preferred upper bound for per-site heavy atom counts.",
    )
    alpha: float = Field(
        default=1.5,
        gt=0,
        description="Soft lower-bound transition width.",
    )
    beta: float = Field(
        default=2.5,
        gt=0,
        description="Soft upper-bound transition width.",
    )
    coverage_weight: float = Field(
        default=0.7,
        ge=0,
        description="Weight for the site coverage term inside site_reward.",
    )
    balance_weight: float = Field(
        default=0.3,
        ge=0,
        description="Weight for the site balance term inside site_reward.",
    )
    lambda_weight: float = Field(
        default=0.85,
        ge=0,
        le=1,
        alias="lambda",
        description="MolE reward weight in final_score = lambda * Mole_reward + (1-lambda) * site_reward.",
    )

    @field_validator("scaffold_smiles")
    @classmethod
    def validate_scaffold_smiles(
        cls, value: Optional[str]
    ) -> Optional[str]:
        if value is None:
            return value
        cleaned = value.strip()
        if not cleaned:
            raise ValueError("scaffold_smiles must be a non-empty string")
        return cleaned

    def validate_configuration(self) -> None:
        if not self.enabled:
            return
        if self.scaffold_smiles is None:
            raise ValueError("site_reward.enabled=true requires scaffold_smiles")
        if self.range_max < self.range_min:
            raise ValueError("site_reward.range_max must be >= site_reward.range_min")
        if (self.coverage_weight + self.balance_weight) <= 0:
            raise ValueError("site_reward coverage_weight + balance_weight must be positive")

    def normalized_component_weights(self) -> tuple[float, float]:
        total = self.coverage_weight + self.balance_weight
        if total <= 0:
            raise ValueError("site_reward component weights must sum to a positive value")
        return (self.coverage_weight / total, self.balance_weight / total)

    def as_payload(self) -> Dict[str, Any]:
        return {
            "enabled": self.enabled,
            "scaffold_smiles": self.scaffold_smiles,
            "range_min": self.range_min,
            "range_max": self.range_max,
            "alpha": self.alpha,
            "beta": self.beta,
            "coverage_weight": self.coverage_weight,
            "balance_weight": self.balance_weight,
            "lambda": self.lambda_weight,
        }


class ScoreObjective(BaseModel):
    """Reward configuration for REINVENT4 scoring."""

    mode: Literal["single_strain", "broad_spectrum_soft"] = Field(
        description="Reward mode for REINVENT4."
    )
    strain: Optional[str] = Field(
        default=None,
        description="Single target strain shortcut for single_strain mode.",
    )
    strains: Optional[List[str]] = Field(
        default=None,
        description="Explicit target strain list. For broad_spectrum_soft, omit to use all strains.",
    )
    weights: Optional[List[float]] = Field(
        default=None,
        description="Optional weights aligned with strains for single_strain mode.",
    )
    tau: float = Field(
        default=0.02,
        gt=0,
        description="Soft threshold temperature for broad_spectrum_soft reward.",
    )
    site_reward: Optional[SiteRewardConfig] = Field(
        default=None,
        description="Optional experimental structural auxiliary reward.",
    )

    @field_validator("strains")
    @classmethod
    def validate_strains(
        cls, value: Optional[List[str]]
    ) -> Optional[List[str]]:
        if value is None:
            return value
        if len(value) == 0:
            raise ValueError("strains list cannot be empty")
        cleaned = [item.strip() for item in value]
        if any(not item for item in cleaned):
            raise ValueError("strains entries must be non-empty strings")
        if len(set(cleaned)) != len(cleaned):
            raise ValueError("strains must be unique")
        return cleaned

    @field_validator("strain")
    @classmethod
    def validate_strain(
        cls, value: Optional[str]
    ) -> Optional[str]:
        if value is None:
            return value
        cleaned = value.strip()
        if not cleaned:
            raise ValueError("strain must be a non-empty string")
        return cleaned

    @field_validator("weights")
    @classmethod
    def validate_weights(
        cls, value: Optional[List[float]]
    ) -> Optional[List[float]]:
        if value is None:
            return value
        if len(value) == 0:
            raise ValueError("weights list cannot be empty")
        if any(weight < 0 for weight in value):
            raise ValueError("weights must be non-negative")
        if sum(value) <= 0:
            raise ValueError("weights must sum to a positive value")
        return value

    def resolved_strains(self) -> Optional[List[str]]:
        """Return the explicit strain panel after applying the strain shortcut."""
        if self.strain and self.strains:
            raise ValueError("Provide either strain or strains, not both")
        if self.strain:
            return [self.strain]
        return self.strains

    def validate_configuration(self) -> None:
        """Validate mode-specific options."""
        selected = self.resolved_strains()
        if self.mode == "single_strain" and not selected:
            raise ValueError("single_strain objective requires strain or strains")
        if self.mode == "broad_spectrum_soft" and self.weights is not None:
            raise ValueError("weights are only supported for single_strain mode")
        if self.site_reward is not None:
            self.site_reward.validate_configuration()

    def normalized_weights(self, panel_size: int) -> List[float]:
        """Return weights aligned to the selected panel."""
        if panel_size <= 0:
            raise ValueError("panel_size must be positive")
        if self.weights is None:
            return [1.0 / panel_size] * panel_size
        if len(self.weights) != panel_size:
            raise ValueError(
                f"weights length ({len(self.weights)}) must match selected strains ({panel_size})"
            )
        weight_sum = sum(self.weights)
        return [weight / weight_sum for weight in self.weights]


class ReinventScoreRequest(BaseModel):
    """FastAPI request model for the dedicated REINVENT4 scoring endpoint."""

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
    objective: ScoreObjective = Field(
        description="Reward definition consumed by REINVENT4."
    )
    app_threshold: float = Field(
        default=0.04374140128493309,
        ge=0,
        le=1,
        description="Shared inhibition threshold used to center broad_spectrum_soft rewards.",
    )
    min_nkill: int = Field(
        default=10,
        ge=0,
        description="Reported for reference in the score response metadata.",
    )

    def to_molecule_input(self) -> MoleculeInput:
        """Convert score requests into the predictor input format."""
        return MoleculeInput(
            molecules=self.molecules,
            smiles=self.smiles,
            chem_id=self.chem_id,
            aggregate_scores=False,
            app_threshold=self.app_threshold,
            min_nkill=self.min_nkill,
        ).normalize()


class StatusModel(BaseModel):
    """Service status response."""

    service: str = "mole-antimicrobial-prediction"
    loaded: bool = False
    device: Optional[str] = None
    version: str = "1.0.0"
    model_path: str = "data/03.model_evaluation/MolE-XGBoost-08.03.2024_14.20.pkl"
    classifier_backend: Optional[str] = None
    classifier_backend_path: Optional[str] = None
    classifier_backend_preference: Optional[str] = None
