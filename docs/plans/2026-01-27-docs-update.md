# Documentation Update Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Update `README.md` and create `README_CN.md` to document the new `MODEL_LOAD_MODE` functionality and auto-tuning behavior.

**Architecture:** Documentation only.

**Tech Stack:** Markdown.

---

### Task 1: Update README.md (English)

**Files:**
- Modify: `README.md`

**Step 1: Update Configuration Section**

Update the `MODEL_LOAD_MODE` explanation to explicitly mention the auto-tuning behavior:
- `resident` (Default):
  - Models load once at startup and stay in VRAM.
  - **Auto-Tuning:** System automatically calculates available VRAM and sets the maximum concurrent request limit to fully utilize the GPU without OOM.
  - Recommended for high-performance production use.
- `on_demand`:
  - Models load/unload for every individual request.
  - No auto-tuning (uses fixed `MAX_CONCURRENT_REQUESTS`).
  - Useful for low-resource environments or when GPU is shared with other heavy tasks.

**Step 2: Commit**

```bash
git add README.md
git commit -m "docs: update README.md with detailed MODEL_LOAD_MODE auto-tuning explanation"
```

### Task 2: Create README_CN.md (Chinese)

**Files:**
- Create: `README_CN.md`

**Step 1: Translate Content**

Create a Chinese version of the `README.md`, ensuring technical terms are accurately translated but `MODEL_LOAD_MODE`, `resident`, `on_demand` etc. remain in English for clarity.
- Project Introduction (项目介绍)
- System Architecture (系统架构 - Mermaid图保持不变)
- Configuration & Performance Tuning (配置与性能调优)
  - 重点解释 `MODEL_LOAD_MODE=resident` 时的“自动显存探测与并发调优”机制。
- Quick Start (快速开始)

**Step 2: Commit**

```bash
git add README_CN.md
git commit -m "docs: add README_CN.md (Chinese documentation)"
```
