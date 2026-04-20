#!/usr/bin/env python3
"""Simple API-mode request test script."""

import json
import os
import sys

import httpx

API_URL = os.getenv("API_URL", "http://localhost:8000")


def _print_response(label: str, response: httpx.Response) -> None:
    print(f"{label} status: {response.status_code}")
    try:
        payload = response.json()
    except json.JSONDecodeError:
        print(response.text)
        return
    print(json.dumps(payload, indent=2, ensure_ascii=False))


def main() -> int:
    try:
        health_resp = httpx.get(f"{API_URL}/health", timeout=10.0)
        health_resp.raise_for_status()
        _print_response("/health", health_resp)
    except Exception as exc:
        print(f"Health check failed: {exc}")
        return 1

    payload = {"smiles": "CCCCS(=N)(=O)CC[C@H](N)C(=O)Oc1ccc(S(=O)(=O)NC(=O)c2ccccc2)cc1", "aggregate_scores": True}
    try:
        predict_resp = httpx.post(
            f"{API_URL}/predict",
            json=payload,
            timeout=120.0,
        )
        predict_resp.raise_for_status()
        _print_response("/predict", predict_resp)
    except Exception as exc:
        print(f"Predict request failed: {exc}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
