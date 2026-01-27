"""FastAPI server that shares the predictor singleton with MCP."""

from fastapi import FastAPI, HTTPException

from src.models import MoleculeInput
from src.service import get_predictor

app = FastAPI(title="Antimicrobial Prediction API")

predictor = get_predictor()


@app.get("/health")
async def health() -> dict:
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
