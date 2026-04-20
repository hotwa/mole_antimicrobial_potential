# Final Cleanup & Documentation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Finalize documentation, provide single-GPU fallback options, and create agent/client examples.

**Architecture:** Documentation and example scripts only.

**Tech Stack:** Docker, Markdown, Python.

---

### Task 1: Create Single-GPU Docker Compose File

**Files:**
- Create: `docker-compose.single-gpu.yml`

**Step 1: Create the file**

Create a simplified version of the current `docker-compose.yml` that:
- Removes the Nginx load balancer.
- Removes `mole_api_1`.
- Keeps only `mole_api_0` (renamed to `mole_api`) and `mole_mcp`.
- Maps port 8000 directly to the container.
- Sets `NVIDIA_VISIBLE_DEVICES=0`.
- Keeps `MODEL_LOAD_MODE=resident` (it works fine on single GPU too).

**Step 2: Commit**

```bash
git add docker-compose.single-gpu.yml
git commit -m "chore: add single-gpu docker-compose file"
```

### Task 2: Create Example Client Script

**Files:**
- Create: `examples/async_predict_client.py`

**Step 1: Implement the script**

Adapt the user's provided Python script into a clean, standalone example.
- Use `httpx` and `asyncio`.
- Include `_calibrate_concurrency`-like logic (client side throttling) or just simple semaphore.
- Add comments explaining how to tune `concurrency` to match the server's capabilities.

**Step 2: Commit**

```bash
git add examples/async_predict_client.py
git commit -m "docs: add async client example script"
```

### Task 3: Create AGENTS.md

**Files:**
- Create: `AGENTS.md`

**Step 1: Write AGENTS.md**

Create a guide specifically for AI agents (like me!) on how to use this tool.
- **Tool Description**: "Antimicrobial Potential Prediction Service".
- **API Specification**:
  - `POST /predict`
  - JSON Schema for Request/Response.
- **MCP Usage**: Explain it exposes an MCP endpoint at `/mcp`.
- **Example Usage**: Reference `examples/async_predict_client.py`.

**Step 2: Commit**

```bash
git add AGENTS.md
git commit -m "docs: add AGENTS.md for AI integration"
```

### Task 4: Update READMEs with Multi-GPU/Single-GPU Guide

**Files:**
- Modify: `README.md`
- Modify: `README_CN.md`

**Step 1: Update README.md (English)**

Add a section "Deployment Options":
- **Multi-GPU (Recommended)**: `docker-compose.yml` (Default).
- **Single GPU**: `docker-compose.single-gpu.yml`.
  - Usage: `docker-compose -f docker-compose.single-gpu.yml up -d`

**Step 2: Update README_CN.md (Chinese)**

Add the same section translated into Chinese.

**Step 3: Commit**

```bash
git add README.md README_CN.md
git commit -m "docs: document single vs multi-gpu deployment options"
```
