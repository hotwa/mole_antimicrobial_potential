# Example Client Script Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Create a unified, robust, and flexible example client script for querying the Antimicrobial Prediction API.

**Architecture:** A single Python script using `httpx` and `asyncio` for high-performance concurrent requests, with `tqdm` for progress tracking.

**Tech Stack:** Python, httpx, pandas, tqdm.

---

### Task 1: Create Consolidated Client Script

**Files:**
- Create: `examples/async_predict_client.py`

**Step 1: Implement the script**

Combine the logic from the user's two snippets into a single file.
- **Imports**: `argparse`, `asyncio`, `httpx`, `pandas`, `tqdm`.
- **Arguments**: Add all CLI args (`--input`, `--output`, `--api-url`, `--concurrency`, `--retries`, etc.) with sensible defaults.
- **Logic**:
  - `_build_payload`: Construct the JSON body.
  - `_fetch_one`: Handle single request with retry loop and backoff.
  - `_run_async`: Manage the `asyncio` task pool and `tqdm` progress bar.
  - `main`: Handle file I/O, column detection, and orchestration.
- **Improvements**: Ensure it handles missing columns gracefully (e.g., auto-generate IDs) and supports flexible input formats.

**Step 2: Commit**

```bash
git add examples/async_predict_client.py
git commit -m "feat: add consolidated async client example script"
```

### Task 2: Create AGENTS.md

**Files:**
- Create: `AGENTS.md`

**Step 1: Write AGENTS.md**

Create the documentation for agents.
- **API Reference**: `POST /predict`.
- **Schema**: Request/Response JSON structure.
- **Client Usage**: How to run `examples/async_predict_client.py`.
  - Example: `python examples/async_predict_client.py --input data.csv --api-url http://localhost:8000`

**Step 2: Commit**

```bash
git add AGENTS.md
git commit -m "docs: add AGENTS.md"
```

### Task 3: Finalize READMEs

**Files:**
- Modify: `README.md`
- Modify: `README_CN.md`

**Step 1: Update README.md (English)**

- Add "Example Client" section referencing `examples/async_predict_client.py`.
- Ensure "Deployment Options" section is present (Single vs Multi GPU).

**Step 2: Update README_CN.md (Chinese)**

- Translate the new sections.

**Step 3: Commit**

```bash
git add README.md README_CN.md
git commit -m "docs: finalize documentation with client script usage"
```
