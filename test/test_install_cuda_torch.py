#!/usr/bin/env python3
"""Tests for the CUDA torch installer helper."""

from __future__ import annotations

from unittest import TestCase

from scripts.install_cuda_torch import build_install_command


class InstallCudaTorchTestCase(TestCase):
    def test_build_install_command_uses_cuda_tag_and_version(self) -> None:
        command = build_install_command(
            torch_version="2.5.1+cu121",
            cuda_tag="cu121",
            extra_packages=("torch-geometric==2.5.0",),
        )

        self.assertEqual(command[:4], ["pip", "install", "--upgrade", "--force-reinstall"])
        self.assertIn("--index-url", command)
        self.assertIn("https://download.pytorch.org/whl/cu121", command)
        self.assertIn("torch==2.5.1+cu121", command)
        self.assertIn("torch-geometric==2.5.0", command)

    def test_build_install_command_supports_other_cuda_tags(self) -> None:
        command = build_install_command(
            torch_version="2.5.1+cu124",
            cuda_tag="cu124",
            extra_packages=(),
        )

        self.assertIn("https://download.pytorch.org/whl/cu124", command)
        self.assertIn("torch==2.5.1+cu124", command)
