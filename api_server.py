"""FastAPI server that shares the predictor singleton with MCP."""

from fastapi import FastAPI, HTTPException

from src.models import MoleculeInput, ReinventScoreRequest
from src.reinvent_scoring import score_reinvent_predictions
from src.service import get_predictor, get_scheduler

app = FastAPI(title="Antimicrobial Prediction API")


@app.get("/health")
async def health() -> dict:
    predictor = get_predictor()
    status = predictor.get_status()
    return {
        "status": "ready" if status.loaded else "loading",
        "service": status.service,
        "version": status.version,
        "loaded": status.loaded,
        "device": status.device,
        "model_path": status.model_path,
        "classifier_backend": status.classifier_backend,
        "classifier_backend_path": status.classifier_backend_path,
        "classifier_backend_preference": status.classifier_backend_preference,
    }


@app.post("/predict")
async def predict(input_data: MoleculeInput) -> dict:
    try:
        scheduler = get_scheduler()
        items = await scheduler.predict_molecules(
            input_data=input_data,
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc

    return {
        "mode": "aggregate" if input_data.aggregate_scores else "per_strain",
        "items": items,
    }


@app.post("/score")
async def score(input_data: ReinventScoreRequest) -> dict:
    try:
        input_data.objective.validate_configuration()
        normalized = input_data.to_molecule_input()
        smiles_by_chem_id = {
            molecule.chem_id: molecule.smiles
            for molecule in normalized.molecules or []
        }
        scheduler = get_scheduler()
        raw_items = await scheduler.predict_molecules(
            input_data=normalized,
            already_normalized=True,
        )
        items = score_reinvent_predictions(
            raw_items,
            input_data,
            smiles_by_chem_id=smiles_by_chem_id,
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc

    return {
        "mode": "reinvent_score",
        "objective": {
            "mode": input_data.objective.mode,
            "strains": input_data.objective.resolved_strains(),
            "app_threshold": input_data.app_threshold,
            "min_nkill": input_data.min_nkill,
            "tau": input_data.objective.tau,
            "site_reward": (
                input_data.objective.site_reward.as_payload()
                if input_data.objective.site_reward is not None
                else None
            ),
        },
        "items": items,
    }
