# QSAR-104 MolE/pMIC Correlation Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` (recommended) or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate a reproducible MolE/XGBoost prediction and correlation analysis for `data/09.qsar_smiles/QSAR-104-SMILES.xlsx`, comparing 104 molecule pMIC values against 40 strain-level MolE probabilities and 3 aggregate antimicrobial scores.

**Architecture:** Add one standalone analysis script that reuses the existing predictor singleton and does not modify core prediction code. The script writes a 44-column prediction matrix, metadata-rich companion matrix, correlation statistics, plots, and a Markdown experiment report under `data/09.qsar_smiles/analysis/`.

**Tech Stack:** pixi, Python 3.10, pandas, NumPy, SciPy, matplotlib, seaborn, existing MolE predictor in `src.service`.

---

## Current Context

Input:

```text
data/09.qsar_smiles/QSAR-104-SMILES.xlsx
```

Known Excel structure:

```text
sheet: 104-smiles_2d
rows: 104
columns: SMILES, name, pMIC
```

Known MolE/XGBoost panel:

```text
strain count: 40
E. coli strains:
- Escherichia coli ED1a (NT5078)
- Escherichia coli IAI1 (NT5077)
```

Important score direction:

```text
antimicrobial_predictive_probability: higher = stronger predicted inhibition
pMIC: higher = stronger real activity
apscore_total/apscore_gnegative/apscore_gpositive: lower = stronger predicted activity
```

Therefore:

```text
strain probabilities should be expected to correlate positively with pMIC
raw apscore values should be expected to correlate negatively with pMIC
plots and directional ranking should use -apscore_* for direction alignment
CSV output must keep the raw apscore_* values
```

---

## Output Contract

Create this output directory:

```text
data/09.qsar_smiles/analysis/
```

Required output files:

```text
data/09.qsar_smiles/analysis/qsar104_mole_prediction_matrix.csv
data/09.qsar_smiles/analysis/qsar104_mole_prediction_matrix_with_metadata.csv
data/09.qsar_smiles/analysis/qsar104_pmic_correlations.csv
data/09.qsar_smiles/analysis/qsar104_ecoli_correlations.csv
data/09.qsar_smiles/analysis/qsar104_feature_summary.csv
data/09.qsar_smiles/analysis/qsar104_analysis_manifest.json
data/09.qsar_smiles/analysis/EXPERIMENT_REPORT.md
data/09.qsar_smiles/analysis/qsar104_trend_original_order_aligned_zscore.png
data/09.qsar_smiles/analysis/qsar104_trend_pmic_sorted_aligned_zscore.png
data/09.qsar_smiles/analysis/qsar104_spearman_correlation_barplot.png
data/09.qsar_smiles/analysis/qsar104_key_scatterplots.png
data/09.qsar_smiles/analysis/qsar104_mole_feature_spearman_heatmap.png
```

Main matrix contract:

```text
file: qsar104_mole_prediction_matrix.csv
shape: 104 rows x 44 columns
column 1: pMIC
columns 2-41: 40 strain antimicrobial_predictive_probability values
columns 42-44: apscore_total, apscore_gnegative, apscore_gpositive
```

Metadata matrix contract:

```text
file: qsar104_mole_prediction_matrix_with_metadata.csv
rows: 104
must include: input_order, chem_id, name, SMILES, pMIC, 40 strain columns, 3 apscore columns
```

Correlation table contract:

```text
file: qsar104_pmic_correlations.csv
rows: 43
must include:
- feature
- n
- pearson_r
- pearson_p
- pearson_q
- pearson_ci95_low
- pearson_ci95_high
- spearman_rho
- spearman_p
- spearman_q
- spearman_ci95_low
- spearman_ci95_high
- kendall_tau
- kendall_p
- kendall_q
- feature_type
- gram_group
- nt_number
- is_ecoli
- expected_correlation_with_pMIC
- directional_pearson_r
- directional_spearman_rho
- direction_matches_expected
- rank_by_directional_spearman
- rank_by_abs_spearman
```

Report contract:

```text
file: EXPERIMENT_REPORT.md
language: Chinese
must include:
- 实验目的
- 实验数据
- 实验方案
- 实验过程
- 结果解读
- E. coli 专项解读
- 分数方向
- 预测值分布
- 结论边界
```

---

## Task 1: Implement Standalone Analysis Script

**Files:**

```text
Create: scripts/analyze_qsar104_mole_correlation.py
Do not modify: src/predictor.py, src/models.py, src/service.py
```

- [ ] **Step 1: Add imports and constants**

Use `matplotlib.use("Agg")` before importing pyplot. Add `REPO_ROOT` to `sys.path` so the script can import project modules when run from the repo root.

Constants:

```python
DEFAULT_INPUT = REPO_ROOT / "data" / "09.qsar_smiles" / "QSAR-104-SMILES.xlsx"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "data" / "09.qsar_smiles" / "analysis"
NT_NUMBER_PATTERN = re.compile(r".*?\((NT\d+)\)")
ECOLI_PATTERN = re.compile(r"escherichia\s+coli|e\.\s*coli|ecoli", re.IGNORECASE)
AGGREGATE_COLUMNS = ["apscore_total", "apscore_gnegative", "apscore_gpositive"]
```

- [ ] **Step 2: Add CLI**

Arguments:

```text
--input default data/09.qsar_smiles/QSAR-104-SMILES.xlsx
--sheet default 0
--output-dir default data/09.qsar_smiles/analysis
--app-threshold default 0.04374140128493309
--min-nkill default 10
--prediction-batch-size default 16
--graph-batch-size default 256
--num-graph-workers default auto
--prefetch-batches default 2
--classifier-workers default auto
--classifier-inflight-batches default auto
--bootstrap-iterations default 1000
--bootstrap-seed default 20260424
--device choices auto,cpu,cuda,cuda:0 default auto
--deterministic-representation action store_true
```

- [ ] **Step 3: Read and validate Excel**

Function: `read_qsar_excel(path: Path, sheet: str | int) -> pd.DataFrame`

Requirements:

```text
read only SMILES, name, pMIC
strip SMILES and name
coerce pMIC to numeric with errors="raise"
fail on missing SMILES/name/pMIC
fail on duplicate name
insert chem_id as QSAR001...QSAR104
insert input_order as 1...104
preserve original row order
```

- [ ] **Step 4: Run per-strain MolE predictions**

Function: `async predict_per_strain(qsar_df, args) -> tuple[pd.DataFrame, list[str], dict[str, str]]`

Use:

```python
from src.models import MoleculeInfo, MoleculeInput
from src.service import get_predictor
```

Implementation requirements:

```text
get predictor with get_predictor()
if --device is not auto, set predictor.device before ensure_loaded()
await predictor.ensure_loaded()
extract strain names from predictor.strain_ohe.index
extract gram dictionary from predictor._gram_dict
iterate QSAR molecules in prediction batches
call predictor.predict(... aggregate_scores=False ...)
pass num_graph_workers, graph_batch_size, prefetch_batches, classifier_workers, classifier_inflight_batches
split pred_id into chem_id and strain_name
assert output row count = 104 * 40
```

- [ ] **Step 5: Build 44-column matrix**

Function: `build_prediction_matrix(qsar_df, pred_df, strain_names, gram_dict)`

Requirements:

```text
pivot per-strain predictions to one row per chem_id
columns must follow predictor strain_names order
calculate apscore_total = mean(log(probability_all_strains))
calculate apscore_gnegative = mean(log(probability_gram_negative_strains))
calculate apscore_gpositive = mean(log(probability_gram_positive_strains))
use NT number in strain name to look up gram stain
return:
  matrix: pMIC + 40 strain columns + 3 apscore columns
  matrix_with_metadata: input_order, chem_id, name, SMILES + matrix columns
  feature_meta: metadata for 43 MolE features
```

Do not include `ginhib_*` or `broad_spectrum` in the 44-column matrix, because the requested 44 groups are `pMIC + 40 probabilities + 3 aggregate scores`.

- [ ] **Step 6: Calculate correlation statistics**

Function: `calculate_correlations(matrix, feature_meta, bootstrap_iterations, bootstrap_seed)`

For each of the 43 MolE features, calculate:

```text
Pearson r and p-value
Spearman rho and p-value
Kendall tau and p-value
bootstrap 95% CI for Pearson
bootstrap 95% CI for Spearman
Benjamini-Hochberg FDR q-values for Pearson/Spearman/Kendall p-values
directional_pearson_r
directional_spearman_rho
direction_matches_expected
rank_by_directional_spearman
rank_by_abs_spearman
```

Expected sign:

```text
strain_probability: positive
aggregate_apscore: negative
```

Directional score:

```text
directional correlation = raw correlation * +1 for strain probabilities
directional correlation = raw correlation * -1 for apscore values
```

- [ ] **Step 7: Build feature summary**

Function: `build_feature_summary(matrix, feature_meta)`

For `pMIC` and all 43 MolE features, report:

```text
n, mean, std, min, q25, median, q75, max
```

Merge feature metadata where applicable.

- [ ] **Step 8: Generate plots**

Plot 1:

```text
qsar104_trend_original_order_aligned_zscore.png
```

Requirements:

```text
z-score pMIC and direction-aligned MolE features
use -apscore_* for plotting only
highlight pMIC as thick black line
highlight two E. coli strains as colored lines
highlight -apscore_total, -apscore_gnegative, -apscore_gpositive as dashed colored lines
draw all other strains as thin gray transparent lines
```

Plot 2:

```text
qsar104_trend_pmic_sorted_aligned_zscore.png
```

Same as plot 1, but sort molecules by ascending pMIC before z-scoring or plotting.

Plot 3:

```text
qsar104_spearman_correlation_barplot.png
```

Requirements:

```text
horizontal bar plot ordered by directional_spearman_rho
E. coli bars orange
aggregate apscore bars green
other strain bars gray
mark spearman_q < 0.05 with an asterisk
label aggregate features as -apscore_* in the plot
```

Plot 4:

```text
qsar104_key_scatterplots.png
```

Requirements:

```text
scatter/regression panels for:
- Escherichia coli ED1a (NT5078)
- Escherichia coli IAI1 (NT5077)
- top directional Spearman feature if not already included
- -apscore_gnegative
- -apscore_total
panel titles must include Spearman rho, Spearman q, and directional Spearman
```

Plot 5:

```text
qsar104_mole_feature_spearman_heatmap.png
```

Requirements:

```text
Spearman correlation heatmap across the 43 direction-aligned MolE features
use -apscore_* columns for the aggregate features
```

- [ ] **Step 9: Generate manifest**

File:

```text
data/09.qsar_smiles/analysis/qsar104_analysis_manifest.json
```

Include:

```text
input path
output directory
n_molecules
n_strains
n_mole_features
app_threshold
min_nkill
bootstrap_iterations
bootstrap_seed
python version
platform
top 10 directional Spearman rows
all output paths
```

- [ ] **Step 10: Generate Chinese Markdown report**

File:

```text
data/09.qsar_smiles/analysis/EXPERIMENT_REPORT.md
```

Must include:

```markdown
# QSAR-104 MolE/pMIC 相关性分析报告

## 实验目的
解释：评估 QSAR-104 的真实 E. coli pMIC 与 MolE/XGBoost 的 40 菌株概率和 3 个聚合分数之间是否存在探索性相关趋势。

## 实验数据
说明输入文件、104 个分子、pMIC 范围、40 个菌株、两个 E. coli 菌株。

## 实验方案
说明 MolE 预测、矩阵构建、apscore 计算、相关性分析、FDR、bootstrap、方向统一。

## 实验过程
记录运行命令、输出目录、输出文件。

## 结果解读
给出 Top directional Spearman 表格和 FDR 显著数量。

### E. coli 专项解读
列出两个 E. coli 菌株的 Pearson/Spearman/q-value/CI/rank。

### 分数方向
明确概率越高越强，apscore 越低越强，图中用 -apscore，CSV 保留原始 apscore。

### 预测值分布
说明 40 个菌株概率范围，提示概率分布过窄会限制相关系数。

## 结论边界
强调这是探索性相关性分析，不能证明 MolE 准确回归 pMIC，也不能把相关性解释为因果。
```

- [ ] **Step 11: Wire main flow**

Main flow:

```text
parse args
create output directory
read Excel
run per-strain prediction
build matrices
calculate correlations
build feature summary
write CSVs
write plots
write manifest
write report
print output directory and key files
```

---

## Task 2: Run Analysis

**Files:**

```text
Read: data/09.qsar_smiles/QSAR-104-SMILES.xlsx
Create: data/09.qsar_smiles/analysis/*
```

- [ ] **Step 1: Check project environment**

Run:

```bash
pixi run mole doctor
```

Expected:

```text
MolE model directory: ... [ok]
Strain panel: primary ok
Pickle backend: ok
Selected classifier backend: pickle
```

- [ ] **Step 2: Run analysis**

Run:

```bash
pixi run python scripts/analyze_qsar104_mole_correlation.py
```

If GPU memory is constrained, rerun:

```bash
pixi run python scripts/analyze_qsar104_mole_correlation.py --device cpu --prediction-batch-size 8 --graph-batch-size 128
```

Expected terminal output:

```text
Wrote QSAR-104 MolE correlation analysis to: .../data/09.qsar_smiles/analysis
Main matrix: .../qsar104_mole_prediction_matrix.csv
Report: .../EXPERIMENT_REPORT.md
```

---

## Task 3: Verify Outputs

- [ ] **Step 1: Verify expected files exist**

Run:

```bash
find data/09.qsar_smiles/analysis -maxdepth 1 -type f | sort
```

Expected files:

```text
data/09.qsar_smiles/analysis/EXPERIMENT_REPORT.md
data/09.qsar_smiles/analysis/qsar104_analysis_manifest.json
data/09.qsar_smiles/analysis/qsar104_ecoli_correlations.csv
data/09.qsar_smiles/analysis/qsar104_feature_summary.csv
data/09.qsar_smiles/analysis/qsar104_key_scatterplots.png
data/09.qsar_smiles/analysis/qsar104_mole_feature_spearman_heatmap.png
data/09.qsar_smiles/analysis/qsar104_mole_prediction_matrix.csv
data/09.qsar_smiles/analysis/qsar104_mole_prediction_matrix_with_metadata.csv
data/09.qsar_smiles/analysis/qsar104_pmic_correlations.csv
data/09.qsar_smiles/analysis/qsar104_spearman_correlation_barplot.png
data/09.qsar_smiles/analysis/qsar104_trend_original_order_aligned_zscore.png
data/09.qsar_smiles/analysis/qsar104_trend_pmic_sorted_aligned_zscore.png
```

- [ ] **Step 2: Verify matrix shape**

Run:

```bash
pixi run python - <<'PY'
import pandas as pd
matrix = pd.read_csv("data/09.qsar_smiles/analysis/qsar104_mole_prediction_matrix.csv")
metadata = pd.read_csv("data/09.qsar_smiles/analysis/qsar104_mole_prediction_matrix_with_metadata.csv")
print("matrix", matrix.shape)
print("metadata", metadata.shape)
print("first columns", matrix.columns[:5].tolist())
print("last columns", matrix.columns[-3:].tolist())
assert matrix.shape == (104, 44)
assert matrix.columns[0] == "pMIC"
assert matrix.columns[-3:].tolist() == ["apscore_total", "apscore_gnegative", "apscore_gpositive"]
assert metadata.shape[0] == 104
PY
```

Expected:

```text
matrix (104, 44)
metadata (104, 48)
last columns ['apscore_total', 'apscore_gnegative', 'apscore_gpositive']
```

- [ ] **Step 3: Verify E. coli rows**

Run:

```bash
pixi run python - <<'PY'
import pandas as pd
ecoli = pd.read_csv("data/09.qsar_smiles/analysis/qsar104_ecoli_correlations.csv")
cols = ["feature", "spearman_rho", "spearman_q", "rank_by_directional_spearman"]
print(ecoli[cols].to_string(index=False))
assert len(ecoli) == 2
assert ecoli["feature"].str.contains("Escherichia coli", case=False, regex=False).all()
PY
```

Expected:

```text
Two rows: Escherichia coli ED1a and Escherichia coli IAI1.
```

- [ ] **Step 4: Verify correlation table**

Run:

```bash
pixi run python - <<'PY'
import pandas as pd
corr = pd.read_csv("data/09.qsar_smiles/analysis/qsar104_pmic_correlations.csv")
required = {
    "feature", "n", "pearson_r", "pearson_p", "pearson_q",
    "spearman_rho", "spearman_p", "spearman_q",
    "spearman_ci95_low", "spearman_ci95_high",
    "kendall_tau", "kendall_p", "kendall_q",
    "feature_type", "gram_group", "nt_number", "is_ecoli",
    "expected_correlation_with_pMIC",
    "directional_pearson_r", "directional_spearman_rho",
    "direction_matches_expected",
    "rank_by_directional_spearman", "rank_by_abs_spearman",
}
missing = required.difference(corr.columns)
print("shape", corr.shape)
print("missing", sorted(missing))
assert corr.shape[0] == 43
assert not missing
PY
```

Expected:

```text
shape (43, ...)
missing []
```

- [ ] **Step 5: Verify report sections**

Run:

```bash
for section in "实验目的" "实验过程" "实验方案" "结果解读" "E. coli 专项解读" "结论边界"; do
  grep -n "$section" data/09.qsar_smiles/analysis/EXPERIMENT_REPORT.md
done
```

Expected:

```text
Each section appears at least once.
```

- [ ] **Step 6: Capture final summary**

Run:

```bash
pixi run python - <<'PY'
import pandas as pd
corr = pd.read_csv("data/09.qsar_smiles/analysis/qsar104_pmic_correlations.csv")
cols = [
    "feature", "feature_type", "is_ecoli",
    "spearman_rho", "spearman_q",
    "directional_spearman_rho", "rank_by_directional_spearman",
]
print("Top directional Spearman:")
print(corr.sort_values("rank_by_directional_spearman")[cols].head(10).to_string(index=False))
print()
print("E. coli:")
print(corr[corr["is_ecoli"]][cols].to_string(index=False))
PY
```

Expected:

```text
Top 10 rows and the two E. coli rows are printed for final interpretation.
```

---

## Final Review Checklist

- [ ] `qsar104_mole_prediction_matrix.csv` has exactly 104 rows and 44 columns.
- [ ] Column 1 is `pMIC`.
- [ ] Columns 2-41 are the 40 Maier strain probability columns.
- [ ] Columns 42-44 are `apscore_total`, `apscore_gnegative`, `apscore_gpositive`.
- [ ] `qsar104_ecoli_correlations.csv` has exactly two rows.
- [ ] `qsar104_pmic_correlations.csv` includes Pearson, Spearman, Kendall, FDR q-values, bootstrap confidence intervals, expected direction, and directional rank.
- [ ] Trend plots use `-apscore_*` for direction alignment, but CSV keeps original `apscore_*`.
- [ ] `EXPERIMENT_REPORT.md` includes experiment purpose, process, plan, and interpretation.
- [ ] Final response reports the exact output folder and key correlation findings without overstating causality.

---

## Recommended Subagent Instruction

Use this instruction for a `gpt-5.5` agent with `xhigh` reasoning:

```text
Use the plan at docs/superpowers/plans/2026-04-24-qsar104-mole-pmic-correlation-analysis.md. Implement it exactly, using the current project pixi environment. Do not modify core predictor behavior. Generate all outputs under data/09.qsar_smiles/analysis/. After running, execute every verification command in Task 3 and report the output file paths, top directional Spearman correlations, and the two E. coli-specific correlation rows.
```
