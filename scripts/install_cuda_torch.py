#!/usr/bin/env python3
"""Install a CUDA-enabled PyTorch wheel into the active Pixi environment."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from typing import Sequence

DEFAULT_TORCH_VERSION = "2.5.1+cu121"
DEFAULT_CUDA_TAG = "cu121"
DEFAULT_INDEX_URL_TEMPLATE = "https://download.pytorch.org/whl/{cuda_tag}"


def _default_cuda_tag(torch_version: str) -> str:
    if "+" not in torch_version:
        return DEFAULT_CUDA_TAG
    return torch_version.split("+", 1)[1]


def build_install_command(
    torch_version: str,
    cuda_tag: str,
    extra_packages: Sequence[str] = (),
    index_url: str | None = None,
) -> list[str]:
    resolved_index_url = index_url or DEFAULT_INDEX_URL_TEMPLATE.format(cuda_tag=cuda_tag)
    command = [
        "pip",
        "install",
        "--upgrade",
        "--force-reinstall",
        "--no-cache-dir",
        "--index-url",
        resolved_index_url,
        f"torch=={torch_version}",
    ]
    command.extend(extra_packages)
    return command


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--torch-version",
        default=os.environ.get("MOLE_TORCH_VERSION", DEFAULT_TORCH_VERSION),
        help="Torch wheel version, including the CUDA suffix (default: %(default)s).",
    )
    parser.add_argument(
        "--cuda-tag",
        default=os.environ.get("MOLE_TORCH_CUDA_TAG"),
        help="PyTorch wheel CUDA tag such as cu121 or cu124 (default: inferred from torch version).",
    )
    parser.add_argument(
        "--index-url",
        default=os.environ.get("MOLE_TORCH_INDEX_URL"),
        help="Custom PyTorch wheel index URL (default: download.pytorch.org for the selected CUDA tag).",
    )
    parser.add_argument(
        "--extra-package",
        action="append",
        default=None,
        help="Optional extra pip package to install together with torch.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the command without installing anything.",
    )
    return parser.parse_args()


def _verify_installation() -> None:
    import torch

    cuda_version = torch.version.cuda
    print(f"torch.__version__={torch.__version__}")
    print(f"torch.version.cuda={cuda_version}")
    print(f"torch.cuda.is_available()={torch.cuda.is_available()}")
    print(f"torch.cuda.device_count()={torch.cuda.device_count()}")
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA torch installation completed but CUDA is still unavailable")


def main() -> int:
    args = _parse_args()
    torch_version = args.torch_version
    cuda_tag = args.cuda_tag or _default_cuda_tag(torch_version)
    extra_packages = tuple(args.extra_package or ())
    command = build_install_command(
        torch_version=torch_version,
        cuda_tag=cuda_tag,
        extra_packages=extra_packages,
        index_url=args.index_url,
    )

    print("Running:", " ".join(command))
    if args.dry_run:
        return 0

    subprocess.run([sys.executable, "-m", *command], check=True)
    _verify_installation()
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
