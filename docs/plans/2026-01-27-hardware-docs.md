# Documentation Update Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Update documentation to explicitly warn users about PyTorch CUDA compatibility limits (sm_50...sm_90) and the consequences of using unsupported hardware (CPU fallback).

**Architecture:** Documentation only.

**Tech Stack:** Markdown.

---

### Task 1: Update README.md (English)

**Files:**
- Modify: `README.md`

**Step 1: Add Hardware Compatibility Section**

Add a warning box or dedicated section "Hardware Compatibility":
- State clearly: The current PyTorch installation supports CUDA compute capabilities **sm_50 to sm_90**.
- **Important Warning**: Newer GPUs like **RTX 5090 (sm_120)** are NOT currently supported by the built-in PyTorch version.
- **Consequence**: If an unsupported GPU is detected, the system will fallback to **CPU mode**, significantly reducing performance.
- **Action**: Users with unsupported GPUs should either upgrade the base image PyTorch version or accept CPU performance.

**Step 2: Commit**

```bash
git add README.md
git commit -m "docs: add GPU compatibility warning to README"
```

### Task 2: Update README_CN.md (Chinese)

**Files:**
- Modify: `README_CN.md`

**Step 1: Add Hardware Compatibility Section (Translated)**

Add the same section in Chinese "硬件兼容性说明":
- 说明当前 PyTorch 版本支持的 CUDA 算力范围 (sm_50 - sm_90)。
- **警告**: RTX 5090 (sm_120) 等最新显卡暂不支持。
- **后果**: 将自动回退到 CPU 模式，导致推理速度变慢。

**Step 2: Commit**

```bash
git add README_CN.md
git commit -m "docs: add GPU compatibility warning to README_CN"
```
