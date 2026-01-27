# README Update Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Create a comprehensive README.md that explains the project to beginners, documents the new high-performance architecture, and provides clear configuration guides.

**Architecture:** Documentation update only.

**Tech Stack:** Markdown, Mermaid Diagrams.

---

### Task 1: Create Comprehensive README.md

**Files:**
- Create: `README.md` (Overwriting existing or creating new if logic dictates, but we will overwrite/update the main README)

**Step 1: Write the Project Introduction (Beginner Friendly)**

Explain that this tool uses AI to predict if a chemical molecule can kill bacteria.
- Input: A chemical structure (SMILES).
- Output: "Yes/No" (Inhibition) and probability scores.
- Analogy: Like a spellchecker but for potential antibiotics.

**Step 2: Add System Architecture Diagram**

Include the Mermaid diagram visualizing:
- Client -> Nginx (Load Balancer)
- Nginx -> GPU Container 0 & GPU Container 1
- GPU Containers -> Resident AI Models

**Step 3: Add Configuration & Performance Tuning Guide**

Document the new environment variables:
- `MODEL_LOAD_MODE`: "resident" (fast, default) vs "on_demand" (slow, memory saving).
- `MAX_CONCURRENT_REQUESTS`: How many parallel calculations per GPU (default 2).
- `WORKERS_API`: Gunicorn workers (keep low if using GPU).

**Step 4: Add Quick Start Guide**

- `docker-compose up -d --build`
- `curl localhost:8000/health`

**Step 5: Commit**

```bash
git add README.md
git commit -m "docs: update README with architecture diagram and configuration guide"
```
