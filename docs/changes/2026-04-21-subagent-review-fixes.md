# Subagent Review Fix Summary

This note summarizes the fixes applied after the independent subagent review of the unified MolE / Timber / CLI / REINVENT4 work.

## What the review found

### 1. Timber auto-backend could fail hard instead of falling back

**Problem**
- `auto` backend selection only checked whether Timber artifacts existed on disk.
- If `TimberCompiledArtifactBackend` failed during initialization, the process could abort even when the pickle/XGBoost backend was available.

**Fix**
- `src/classifier_backend.py`
  - wrapped Timber initialization in fallback-aware selection logic
  - `auto` now falls back to the pickle backend when Timber initialization fails
  - explicit `timber` preference still fails fast, which is correct

### 2. Default pickle model path was not repo-root relative

**Problem**
- The fallback pickle path was effectively resolved from the current working directory.
- Running the CLI or API from a different directory could point at the wrong `data/...` path.

**Fix**
- `src/classifier_backend.py`
  - resolved the default pickle model path relative to the repository root
  - the app no longer depends on the caller's current directory for the default fallback path

### 3. Duplicate `chem_id` values were silently accepted

**Problem**
- `/score` accepted duplicate `chem_id` values.
- Downstream parsing keyed results only by `chem_id`, so later rows could overwrite earlier rows without warning.

**Fix**
- `src/models.py`
  - added duplicate `chem_id` validation during request normalization
  - duplicate IDs now return a validation error instead of silently corrupting the output

## Resulting behavior

- Timber remains the preferred accelerated backend when available and valid.
- The app is safer to run from any working directory.
- `/score` now rejects ambiguous batch inputs instead of silently losing data.

## Verification

The following checks passed after the fixes:

- `pixi run test-cli`
- `pixi run test-classifier-backend`
- `pixi run test-score`
- `pixi run test-reinvent4-workflow`
- `MOLE_MOLE_MODEL_PATH=/home/lingyuzeng/project/mole_antimicrobial_potential/pretrained_model/model_ginconcat_btwin_100k_d8000_l0.0001 pixi run mole doctor`
- `MOLE_MOLE_MODEL_PATH=/home/lingyuzeng/project/mole_antimicrobial_potential/pretrained_model/model_ginconcat_btwin_100k_d8000_l0.0001 pixi run mole predict --smiles CCO`
- `MOLE_MOLE_MODEL_PATH=/home/lingyuzeng/project/mole_antimicrobial_potential/pretrained_model/model_ginconcat_btwin_100k_d8000_l0.0001 pixi run mole score --objective-file workflows/reinvent4/inputs/objectives/pathogen_group_a.site_reward.prototype.json --smiles CCO`

## Commit reference

The fixes were committed as:

- `bdf98e7` `fix: harden backend fallback and request validation`

