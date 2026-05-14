# 还原胺化产物 MolE 与 pMIC_pred 相关性分析报告

## 实验背景

本分析针对泰乐菌素（Tylosin）和泰万菌素（Tylvalosin）两个16元大环内酯母核通过还原胺化反应拼接外源代谢物的产物，评估两种预测模型（MolE 和 2D-QSAR）输出之间的相关性。

## 数据来源

- **输入目录**: `data/11.macrolide_reductive_amination/`
- **分析目录**: `data/11.macrolide_reductive_amination/analysis/`

### 输入文件

| 文件名 | 母核 | 产物数量 |
|--------|------|----------|
| `tylosin_seed_amine_products_with_mole_predictions.csv` | 泰乐菌素 | 114 |
| `tylvalosin_seed_amine_products_with_mole_predictions.csv` | 泰万菌素 | 114 |

### 预测模型说明

| 模型 | 指标 | 含义 | 数值方向 |
|------|------|------|----------|
| 2D-QSAR | `pMIC_pred（μmol/ml）` | 预测的抗菌活性（pMIC） | 越大 = 活性越强 |
| MolE | `apscore_total` | 广谱抗菌潜力分数 | 越接近0 = 活性越强 |
| MolE | `ginhib_total` | 预测抑制菌株数 | 越大 = 活性越强 |

**注意**: MolE 的 apscore 是 `mean(log(probability))`，由于 probability 在 0-1 之间，log 后为负值。因此：
- apscore 接近 0 → 几何平均概率高 → 预测活性强
- apscore 很负 → 几何平均概率低 → 预测活性弱

## 相关性分析结果

### 泰乐菌素（Tylosin）

| 特征 | Pearson r | Pearson p | Spearman rho | Spearman p | 显著性 |
|------|-----------|-----------|--------------|------------|--------|
| apscore_total | +0.053 | 0.577 | -0.156 | 0.098 | ns |
| apscore_gnegative | +0.066 | 0.487 | -0.141 | 0.135 | ns |
| apscore_gpositive | +0.037 | 0.695 | -0.173 | 0.065 | ns |
| ginhib_total | +0.111 | 0.241 | -0.010 | 0.915 | ns |
| ginhib_gnegative | +0.140 | 0.138 | +0.030 | 0.750 | ns |
| ginhib_gpositive | +0.034 | 0.717 | -0.089 | 0.344 | ns |
| broad_spectrum | NaN | NaN | NaN | NaN | 常量 |

**结论**: 泰乐菌素产物中，MolE 预测与 pMIC_pred **无显著相关性**（所有 p > 0.05）。

### 泰万菌素（Tylvalosin）

| 特征 | Pearson r | Pearson p | Spearman rho | Spearman p | 显著性 |
|------|-----------|-----------|--------------|------------|--------|
| apscore_total | -0.516 | 4.3e-09 | -0.452 | 4.4e-07 | *** |
| apscore_gnegative | -0.441 | 8.9e-07 | -0.394 | 1.5e-05 | *** |
| apscore_gpositive | -0.568 | 4.5e-11 | -0.481 | 6.0e-08 | *** |
| ginhib_total | -0.399 | 1.1e-05 | -0.362 | 7.5e-05 | *** |
| ginhib_gnegative | -0.300 | 1.2e-03 | -0.268 | 4.0e-03 | ** |
| ginhib_gpositive | -0.382 | 2.7e-05 | -0.352 | 1.2e-04 | *** |
| broad_spectrum | NaN | NaN | NaN | NaN | 常量 |

显著性标记: *** p < 0.001, ** p < 0.01, * p < 0.05, ns 不显著

**结论**: 泰万菌素产物中，MolE 预测与 pMIC_pred 呈**显著负相关**（Spearman rho ≈ -0.45, p < 0.001）。

## 负相关的解释

### 分数方向对比

| 模型 | pMIC_pred高（活性强）时 | pMIC_pred低（活性弱）时 |
|------|------------------------|------------------------|
| 2D-QSAR | pMIC_pred 值大（接近0） | pMIC_pred 值小（很负） |
| MolE | apscore 更负（活性弱） | apscore 接近0（活性强） |

**两个模型的预测方向相反。**

### 数据验证

以泰万菌素为例：

| 产物 | pMIC_pred | apscore_total | 解释 |
|------|-----------|---------------|------|
| 0070 (pMIC最高) | **-0.21** (活性强) | **-3.12** (活性弱) | 模型预测相反 |
| 0034 (pMIC最低) | **-4.07** (活性弱) | **-1.88** (活性强) | 模型预测相反 |

### 可能原因

1. **训练数据不同**
   - 2D-QSAR: 基于104个化合物的实验MIC数据
   - MolE: 基于Maier等人40株细菌的筛选数据

2. **预测目标不同**
   - pMIC_pred: 针对特定细菌（如E. coli）的MIC
   - MolE: 广谱抗菌潜力（40株细菌的几何平均）

3. **模型架构不同**
   - 2D-QSAR: MOE分子描述符 + 回归模型
   - MolE: 图神经网络 + XGBoost分类器

## 区分度分析

| 指标 | 泰乐菌素 | 泰万菌素 |
|------|----------|----------|
| pMIC_pred变异系数 | **99.6%** | **48.6%** |
| MolE apscore变异系数 | 36.6% | **11.0%** |

**pMIC_pred 对化合物的区分度更高**，能更好地区分化合物优劣。

## 交集分析

### Top 10 交集

| 母核 | pMIC_pred Top 10 | MolE Top 10 | 交集数量 | 交集化合物 |
|------|------------------|-------------|----------|-----------|
| 泰乐菌素 | ✓ | ✓ | 1 | tylosin_seed_amine_0001 (ADC) |
| 泰万菌素 | ✓ | ✓ | 1 | tylvalosin_seed_amine_0001 (ADC) |

### Top 20 交集

| 母核 | 交集数量 | 交集化合物 |
|------|----------|-----------|
| 泰乐菌素 | 2 | 0001 (ADC), 0046 (Pretyrosine) |
| 泰万菌素 | 2 | 0001 (ADC), 0046 (Pretyrosine) |

### 交集化合物详情

#### tylosin_seed_amine_0001 (ADC)
- **Seed代谢物**: ADC
- **pMIC_pred**: +0.1305 (Top 1%)
- **apscore_total**: -1.038 (Top 9%)
- **ginhib_total**: 39/40

#### tylvalosin_seed_amine_0001 (ADC)
- **Seed代谢物**: ADC
- **pMIC_pred**: -0.5187 (Top 1%)
- **apscore_total**: -1.877 (Top 9%)
- **ginhib_total**: 39/40

#### tylosin_seed_amine_0046 (Pretyrosine)
- **Seed代谢物**: Pretyrosine
- **pMIC_pred**: +0.056 (Top 5%)
- **apscore_total**: -1.284 (Top 20%)

#### tylvalosin_seed_amine_0046 (Pretyrosine)
- **Seed代谢物**: Pretyrosine
- **pMIC_pred**: -0.593 (Top 5%)
- **apscore_total**: -2.441 (非Top 20%)

## 实验筛选建议

### 模型选择建议

**推荐以 pMIC_pred 为主，MolE 为辅**

理由：
1. pMIC_pred 区分度更强（变异系数更大）
2. MolE 区分度较弱（泰万菌素变异系数仅11%）
3. 训练数据相关性：若2D-QSAR包含大环内酯类化合物，预测更可靠

### 筛选方案

#### 方案A：保守策略（推荐）
选择 pMIC_pred Top 10-20 中的化合物，覆盖不同结构类型的 seed 代谢物。

**泰乐菌素推荐**:
1. tylosin_seed_amine_0070 - (S)-alpha-Amino-beta
2. tylosin_seed_amine_0001 - ADC
3. tylosin_seed_amine_0086 - Gly
4. tylosin_seed_amine_0046 - Pretyrosine
5. tylosin_seed_amine_0088 - (S)-2,3-Diaminopropa

**泰万菌素推荐**:
1. tylvalosin_seed_amine_0070 - (S)-alpha-Amino-beta
2. tylvalosin_seed_amine_0001 - ADC
3. tylvalosin_seed_amine_0086 - Gly
4. tylvalosin_seed_amine_0046 - Pretyrosine
5. tylvalosin_seed_amine_0088 - (S)-2,3-Diaminopropa

#### 方案B：交集策略
选择两个模型都排名靠前的化合物，即交集化合物。

**必选**:
- tylosin_seed_amine_0001 (ADC) - 两个母核都入选
- tylvalosin_seed_amine_0001 (ADC) - 两个母核都入选

#### 方案C：多样性策略（推荐用于初步探索）
- pMIC_pred Top 3
- MolE Top 3
- 交集化合物
- 随机抽样 2-3 个作为对照

### 需要确认的问题

1. **2D-QSAR训练数据**是否包含大环内酯类化合物？
2. **实验目标菌**是什么？如果是E. coli，pMIC_pred可能更相关
3. **实验预算**能测试多少个化合物？

## 输出文件

| 文件 | 说明 |
|------|------|
| `combined_correlations.csv` | 所有相关性统计结果 |
| `tylosin_seed_amine_products_correlation_heatmap.png` | 泰乐菌素相关性热图 |
| `tylosin_seed_amine_products_scatter_plots.png` | 泰乐菌素散点图 |
| `tylvalosin_seed_amine_products_correlation_heatmap.png` | 泰万菌素相关性热图 |
| `tylvalosin_seed_amine_products_scatter_plots.png` | 泰万菌素散点图 |
| `ANALYSIS_REPORT.md` | 本报告 |

## 脚本

分析脚本: `scripts/analyze_reductive_amination_correlation.py`

运行命令:
```bash
pixi run python scripts/analyze_reductive_amination_correlation.py
```

---

**报告生成日期**: 2026-05-14
**分析工具**: Python 3.10, scipy, pandas, matplotlib, seaborn
