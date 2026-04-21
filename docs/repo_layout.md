# Repository Layout

This repository now keeps only a few canonical entrypoints at the root. Most
implementation logic lives under `src/`, and older compatibility scripts live
under `scripts/legacy/`.

## Canonical Entrypoints

Use these first for new work:

- [`mole_cli.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/mole_cli.py)
  - Unified CLI for MolE, Timber, and REINVENT4
  - Recommended usage: `pixi run mole ...`
- [`api_server.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/api_server.py)
  - FastAPI app exposing `/health`, `/predict`, and `/score`
  - Recommended usage: `python api_server.py`
- [`mcp_server_enhanced.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/mcp_server_enhanced.py)
  - Streamable HTTP FastMCP server
  - Recommended usage: `python mcp_server_enhanced.py`

## Core Library Modules

Import these from other code instead of duplicating behavior:

- [`src/mole_representation.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/mole_representation.py)
  - MolE SMILES-to-embedding pipeline
- [`src/predictor.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/predictor.py)
  - Shared antimicrobial predictor service
- [`src/classifier_backend.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/classifier_backend.py)
  - Timber/pickle backend selection
- [`src/models.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/models.py)
  - Request schemas, objective normalization, and validation
- [`src/reinvent_scoring.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/reinvent_scoring.py)
  - REINVENT4 reward computation
- [`src/site_reward.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/site_reward.py)
  - Optional site-aware auxiliary reward prototype
- [`src/reinvent4_workflow.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/reinvent4_workflow.py)
  - REINVENT4 config rendering, chunk orchestration, and validation
- [`src/service.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/src/service.py)
  - Predictor singleton shared by API and MCP servers

## Legacy Compatibility Scripts

These scripts are kept for reproduction or reference. Prefer the canonical
entrypoints above for new work.

- [`scripts/legacy/mole_representation.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/scripts/legacy/mole_representation.py)
  - Legacy MolE embedding CLI
- [`scripts/legacy/mole_antimicrobial_prediction.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/scripts/legacy/mole_antimicrobial_prediction.py)
  - Legacy TSV batch predictor
- [`scripts/legacy/predict_api.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/scripts/legacy/predict_api.py)
  - Legacy FastAPI app
- [`scripts/legacy/mcp_server.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/scripts/legacy/mcp_server.py)
  - Legacy FastMCP app
- [`scripts/legacy/health_check.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/scripts/legacy/health_check.py)
  - Legacy health probe helper

## REINVENT4 Workflow Scripts

The optimization workflow lives under:

- [`workflows/reinvent4/scripts/`](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts)

Main utilities:

- [`run_reinvent4.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/run_reinvent4.py)
  - Render config and launch a single REINVENT4 run
- [`run_long_rl.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/run_long_rl.py)
  - Chunked long-run RL orchestration
- [`score_bridge.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/score_bridge.py)
  - Bridge REINVENT4 batches to local `/score`
- [`evaluate_results.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/evaluate_results.py)
  - Re-score and summarize generated molecules
- [`analyze_rgroup_sites.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/analyze_rgroup_sites.py)
  - Inspect per-site decoration sizes
- [`validate_scaffold.py`](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/validate_scaffold.py)
  - Validate scaffold attachment points

## Quick Rules

- Use `pixi run mole doctor` before running anything on a new machine.
- Use `pixi run mole predict` and `pixi run mole score` for day-to-day work.
- Use `python api_server.py` only when you specifically need the HTTP server.
- Use `python mcp_server_enhanced.py` only when you specifically need MCP.
- Use the scripts in `scripts/legacy/` only when reproducing older flows.
- Use `workflows/reinvent4/scripts/` for all REINVENT4-specific automation.

