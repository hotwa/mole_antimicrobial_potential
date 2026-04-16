#!/usr/bin/env python3
"""Bridge REINVENT4 ExternalProcess scoring to the local /score API."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.reinvent4_workflow import (  # noqa: E402
    build_score_request,
    external_process_payload,
    http_json_request,
    load_objective_spec,
    write_json_file,
)
from src.site_reward import scaffold_from_file  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score-url", required=True, help="Absolute URL of POST /score")
    parser.add_argument(
        "--objective-file",
        required=True,
        help="Path to the objective JSON file consumed by /score",
    )
    parser.add_argument(
        "--scaffold-file",
        default="",
        help="Optional scaffold .smi file used to fill site_reward.scaffold_smiles.",
    )
    parser.add_argument(
        "--audit-file",
        default="",
        help="Optional JSON file or JSONL path where bridge requests are recorded",
    )
    parser.add_argument(
        "--request-timeout",
        type=int,
        default=120,
        help="HTTP timeout in seconds for /score requests",
    )
    return parser.parse_args()


def read_smiles_from_stdin() -> list[str]:
    return [line.strip() for line in sys.stdin if line.strip()]


def append_audit_record(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".jsonl":
        with path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(payload, ensure_ascii=False) + "\n")
        return
    write_json_file(path, payload)


def main() -> int:
    args = parse_args()
    smiles = read_smiles_from_stdin()
    if not smiles:
        raise SystemExit("No SMILES received on stdin")

    objective = load_objective_spec(args.objective_file)
    site_reward = objective.get("site_reward")
    if site_reward and site_reward.get("enabled") and not site_reward.get("scaffold_smiles"):
        if not args.scaffold_file:
            raise SystemExit("site_reward.enabled=true requires --scaffold-file or site_reward.scaffold_smiles")
        objective["site_reward"] = {
            **site_reward,
            "scaffold_smiles": scaffold_from_file(args.scaffold_file),
        }
    request_body = build_score_request(smiles, objective)
    response_body = http_json_request(
        args.score_url,
        request_body,
        timeout_seconds=args.request_timeout,
    )
    payload = external_process_payload(
        response_body,
        request_body["chem_id"],
    )

    if args.audit_file:
        append_audit_record(
            Path(args.audit_file).expanduser(),
            {
                "request": request_body,
                "response": response_body,
                "external_payload": payload,
            },
        )

    sys.stdout.write(json.dumps(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
