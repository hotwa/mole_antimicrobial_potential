#!/usr/bin/env python3
"""Regression tests for MolE representation helpers."""

from __future__ import annotations

import tempfile
from pathlib import Path
from unittest import TestCase
from unittest.mock import MagicMock, patch

from src.mole_representation import load_pretrained_model


class LoadPretrainedModelTestCase(TestCase):
    def test_load_pretrained_model_closes_config_file_via_context_manager(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            model_dir = Path(tmpdir)
            (model_dir / "config.yaml").write_text("model: {}\n", encoding="utf-8")
            (model_dir / "model.pth").write_bytes(b"unused")

            config_stream = MagicMock(name="config_stream")
            entered_stream = MagicMock(name="entered_stream")
            config_stream.__enter__.return_value = entered_stream

            fake_model = MagicMock(name="fake_model")
            fake_model.to.return_value = fake_model

            with patch("builtins.open", return_value=config_stream) as open_mock, patch(
                "src.mole_representation.yaml.load",
                return_value={"model": {}},
            ) as yaml_load_mock, patch(
                "src.mole_representation.GINet",
                return_value=fake_model,
            ), patch(
                "src.mole_representation.torch.load",
                return_value={},
            ):
                load_pretrained_model(str(model_dir), device="cpu")

            open_mock.assert_called_once_with(str(model_dir / "config.yaml"), "r", encoding="utf-8")
            self.assertIs(
                yaml_load_mock.call_args.args[0],
                entered_stream,
                "yaml.load() should read the context-managed config stream",
            )
            config_stream.__exit__.assert_called_once()
