#!/usr/bin/env python3
"""Enhanced FastMCP server (streamable HTTP)."""

from __future__ import annotations

import os
import time
from typing import Any, Dict, List, Optional, Union

import pandas as pd
from fastmcp import FastMCP
from pydantic import Field

from src.models import MoleculeInput
from src.service import get_predictor, get_scheduler

mcp = FastMCP("mole-antimicrobial-prediction")

VERSION = "1.0.0"
SERVICE_NAME = "mole-antimicrobial-prediction"
STRAIN_FILE = "data/01.prepare_training_data/maier_screening_results.tsv.gz"

_strain_cache: Optional[List[str]] = None


def _load_strain_names() -> List[str]:
    global _strain_cache
    if _strain_cache is not None:
        return _strain_cache
    if not os.path.exists(STRAIN_FILE):
        _strain_cache = []
        return _strain_cache
    data = pd.read_csv(STRAIN_FILE, sep="\t", index_col=0)
    _strain_cache = list(data.columns)
    return _strain_cache


@mcp.tool
async def status() -> Dict[str, Any]:
    """Return service status without forcing model load."""
    predictor = get_predictor()
    status_info = predictor.get_status()
    return {
        "service": SERVICE_NAME,
        "version": VERSION,
        "loaded": status_info.loaded,
        "device": status_info.device,
        "model_path": status_info.model_path,
    }


@mcp.tool
async def predict_antimicrobial_potential(
    molecules: Optional[List[Dict[str, Any]]] = Field(
        default=None,
        description="List of {smiles, chem_id}. Preferred structured input.",
    ),
    smiles: Optional[Union[str, List[str]]] = Field(
        default=None,
        description="SMILES string or list of SMILES strings.",
    ),
    chem_id: Optional[Union[str, List[str]]] = Field(
        default=None,
        description="Compound id or list of ids aligned with smiles.",
    ),
    aggregate_scores: bool = Field(
        default=False,
        description="If true, return aggregated antimicrobial potential scores.",
    ),
    app_threshold: float = Field(
        default=0.04374140128493309,
        ge=0,
        le=1,
        description="Threshold (0-1) for binarizing growth inhibition.",
    ),
    min_nkill: int = Field(
        default=10,
        ge=0,
        description="Minimum inhibited strain count to label broad_spectrum=1.",
    ),
) -> Dict[str, Any]:
    """Predict antimicrobial potential for molecules."""
    start_time = time.monotonic()

    input_data = MoleculeInput(
        molecules=molecules,
        smiles=smiles,
        chem_id=chem_id,
        aggregate_scores=aggregate_scores,
        app_threshold=app_threshold,
        min_nkill=min_nkill,
    )

    scheduler = get_scheduler()
    items = await scheduler.predict_molecules(
        input_data=input_data,
    )

    latency_ms = int((time.monotonic() - start_time) * 1000)
    # Status comes from predictor
    status_info = get_predictor().get_status()

    return {
        "mode": "aggregate" if input_data.aggregate_scores else "per_strain",
        "items": items,
        "meta": {
            "service": SERVICE_NAME,
            "version": VERSION,
            "loaded": status_info.loaded,
            "device": status_info.device,
            "latency_ms": latency_ms,
        },
    }


@mcp.resource("resource://schema")
async def get_schema_info() -> Dict[str, Any]:
    """Schema and usage notes."""
    return {
        "service": SERVICE_NAME,
        "version": VERSION,
        "tools": {
            "status": {
                "description": "Check service status and load state",
                "returns": ["service", "version", "loaded", "device", "model_path"],
            },
            "predict_antimicrobial_potential": {
                "description": "Predict antimicrobial activity",
                "input": ["molecules", "smiles", "chem_id"],
                "options": ["aggregate_scores", "app_threshold", "min_nkill"],
                "returns": ["mode", "items", "meta"],
            },
        },
        "examples": {
            "tools/list": {
                "jsonrpc": "2.0",
                "id": "1",
                "method": "tools/list",
            },
            "predict": {
                "jsonrpc": "2.0",
                "id": "2",
                "method": "tools/call",
                "params": {
                    "name": "predict_antimicrobial_potential",
                    "arguments": {"smiles": "CCO", "aggregate_scores": True},
                },
            },
        },
    }


@mcp.resource("resource://about")
async def get_about() -> Dict[str, Any]:
    """Service metadata."""
    return {
        "service": SERVICE_NAME,
        "version": VERSION,
        "description": "AI-powered antimicrobial activity prediction using MolE + XGBoost",
        "paper": {
            "title": "Pre-trained molecular representations enable antimicrobial discovery",
            "journal": "Nature Communications",
            "year": 2025,
            "url": "https://www.nature.com/articles/s41467-025-58804-4",
        },
        "training_data": {
            "source": "Maier et al. (2018)",
            "journal": "Nature",
            "strains": "100+ bacterial strains",
        },
        "framework": {
            "mole": "https://github.com/rolayoalarcon/MolE",
            "xgboost": "XGBoost classifier",
        },
    }


@mcp.resource("resource://strains")
async def get_strains_summary() -> Dict[str, Any]:
    """Strain summary (total + sample)."""
    strains = _load_strain_names()
    sample = strains[:5]
    return {
        "total": len(strains),
        "sample": sample,
        "source": "Maier et al. (2018)",
    }


@mcp.resource("resource://strains/all")
async def get_strains_all() -> Dict[str, Any]:
    """Full strain list."""
    strains = _load_strain_names()
    return {"total": len(strains), "strains": strains}


def main() -> None:
    """Start the FastMCP server."""
    host = os.getenv("HOST", "0.0.0.0")
    port = int(os.getenv("PORT_MCP", "8001"))
    transport = os.getenv("MCP_TRANSPORT", "http")
    json_response = os.getenv("MCP_JSON_RESPONSE", "true").lower() == "true"

    print("[INFO] Starting Enhanced FastMCP server")
    print(f"[INFO] Service: {SERVICE_NAME}")
    print(f"[INFO] Version: {VERSION}")
    print(f"[INFO] Listening on {host}:{port} via {transport}")

    mcp.run(transport=transport, host=host, port=port, json_response=json_response)


if __name__ == "__main__":
    main()
