#!/usr/bin/env python3
"""Validate a LibInvent scaffold file before launching REINVENT4."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.reinvent4_workflow import validate_scaffold_file  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scaffold_file", help="Path to the LibInvent scaffold file")
    return parser.parse_args()


def main() -> int:
    records = validate_scaffold_file(Path(parse_args().scaffold_file))
    sys.stdout.write(json.dumps(records, indent=2, ensure_ascii=False) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
