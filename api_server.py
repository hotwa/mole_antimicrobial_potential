"""FastAPI server that shares the predictor singleton with MCP."""

from fastapi import FastAPI, HTTPException

from src.models import MoleculeInput, ReinventScoreRequest
from src.reinvent_scoring import score_reinvent_predictions
from src.service import get_predictor

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
    }


@app.post("/predict")
async def predict(input_data: MoleculeInput) -> dict:
    try:
        normalized = input_data.normalize()
        predictor = get_predictor()
        await predictor.ensure_loaded()
        items = await predictor.predict(normalized)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc

    return {
        "mode": "aggregate" if normalized.aggregate_scores else "per_strain",
        "items": items,
    }


@app.post("/score")
async def score(input_data: ReinventScoreRequest) -> dict:
    try:
        input_data.objective.validate_configuration()
        normalized = input_data.to_molecule_input()
        predictor = get_predictor()
        await predictor.ensure_loaded()
        raw_items = await predictor.predict(normalized)
        items = score_reinvent_predictions(raw_items, input_data)
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
        },
        "items": items,
    }
