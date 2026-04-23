from __future__ import annotations

import unittest
import importlib
import sys
from unittest import mock

from fastapi.testclient import TestClient

import api_server
from src.models import MoleculeInput, ReinventScoreRequest


class _RecordingScheduler:
    def __init__(self, return_value=None):
        self.calls = []
        self.return_value = return_value or []

    async def predict_molecules(self, *args, **kwargs):
        self.calls.append((args, kwargs))
        return self.return_value


class ApiNormalizationPathTestCase(unittest.TestCase):
    def test_predict_endpoint_does_not_normalize_before_scheduler(self) -> None:
        scheduler = _RecordingScheduler(return_value=[{"chem_id": "mol1", "apscore_total": -1.2}])
        client = TestClient(api_server.app)

        with mock.patch("api_server.get_scheduler", return_value=scheduler), mock.patch.object(
            MoleculeInput,
            "normalize",
            autospec=True,
            wraps=MoleculeInput.normalize,
        ) as normalize_mock:
            response = client.post(
                "/predict",
                json={"smiles": "CCO", "aggregate_scores": True},
            )

        self.assertEqual(response.status_code, 200, response.text)
        self.assertEqual(normalize_mock.call_count, 0)
        _, kwargs = scheduler.calls[0]
        self.assertIn("input_data", kwargs)
        self.assertIsInstance(kwargs["input_data"], MoleculeInput)

    def test_mcp_tool_does_not_normalize_before_scheduler(self) -> None:
        scheduler = _RecordingScheduler(return_value=[{"chem_id": "mol1", "apscore_total": -1.2}])

        async def _run() -> None:
            def fake_tool(*args, **kwargs):
                if args and callable(args[0]) and len(args) == 1 and not kwargs:
                    return args[0]
                if len(args) >= 2 and callable(args[1]):
                    return args[1]
                return lambda func: func

            with mock.patch("fastmcp.FastMCP.tool", new=fake_tool), mock.patch(
                "fastmcp.FastMCP.resource", new=fake_tool
            ):
                sys.modules.pop("mcp_server_enhanced", None)
                mcp_server_enhanced = importlib.import_module("mcp_server_enhanced")

            with mock.patch("mcp_server_enhanced.get_scheduler", return_value=scheduler), mock.patch.object(
                MoleculeInput,
                "normalize",
                autospec=True,
                wraps=MoleculeInput.normalize,
            ) as normalize_mock:
                result = await mcp_server_enhanced.predict_antimicrobial_potential(
                    molecules=None,
                    smiles="CCO",
                    chem_id=None,
                    aggregate_scores=True,
                    app_threshold=0.04374140128493309,
                    min_nkill=10,
                )

            self.assertEqual(result["mode"], "aggregate")
            self.assertEqual(normalize_mock.call_count, 0)
            _, kwargs = scheduler.calls[0]
            self.assertIn("input_data", kwargs)
            self.assertIsInstance(kwargs["input_data"], MoleculeInput)

        import asyncio

        asyncio.run(_run())

    def test_score_endpoint_converts_request_once(self) -> None:
        scheduler = _RecordingScheduler(
            return_value=[
                {
                    "pred_id": "mol1:Strain A (NT1)",
                    "antimicrobial_predictive_probability": 0.8,
                    "growth_inhibition": 1,
                }
            ]
        )
        client = TestClient(api_server.app)

        with mock.patch("api_server.get_scheduler", return_value=scheduler), mock.patch.object(
            ReinventScoreRequest,
            "to_molecule_input",
            autospec=True,
            side_effect=ReinventScoreRequest.to_molecule_input,
        ) as to_molecule_input_mock:
            response = client.post(
                "/score",
                json={
                    "smiles": "CCO",
                    "chem_id": "mol1",
                    "objective": {
                        "mode": "single_strain",
                        "strain": "Strain A (NT1)",
                    },
                },
            )

        self.assertEqual(response.status_code, 200, response.text)
        self.assertEqual(to_molecule_input_mock.call_count, 1)
        _, kwargs = scheduler.calls[0]
        self.assertTrue(kwargs["already_normalized"])
        self.assertIsInstance(kwargs["input_data"], MoleculeInput)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
