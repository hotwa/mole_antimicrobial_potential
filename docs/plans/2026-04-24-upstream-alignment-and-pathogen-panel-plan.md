# Upstream 对齐与致病菌定向打分升级计划

**日期**: 2026-04-24  
**范围**: `mole_antimicrobial_potential` 主仓库（CLI/API/MCP + batch screen）  
**目标**:  
1) 明确是否与上游对齐；  
2) 落地“按指定致病菌集合打分”；  
3) 保留“少伤共生菌”的选择性评价能力。

---

## 0. 现状结论（已核对）

与上游 `rolayoalarcon/mole_antimicrobial_potential` 对比：

- `app_threshold` 默认值一致：`0.04374140128493309`
- `min_nkill` 默认值一致：`10`
- `growth_inhibition` 判定一致：`p >= app_threshold`
- `ginhib_total` 计算一致：对菌株二值抑制求和
- `broad_spectrum` 判定一致：`ginhib_total >= min_nkill`
- 唯一差异：`apscore_*` 对数底  
  - 上游：`log2(gmean(p))`
  - 当前：`log(gmean(p))`

---

## 1. 决策建议

### 建议 1：与上游对齐（建议执行）

建议将 `apscore_*` 的主口径对齐到上游 `log2`，原因：

- 便于与上游结果和文献口径直接比较；
- 不影响排序单调性（`log` 与 `log2` 是线性换底）；
- 不影响 `ginhib_*` / `broad_spectrum` 判定。

### 建议 2：致病菌集合打分优先于单纯提高全局阈值（建议执行）

当前筛选结果中 `ginhib_total` 分布已偏高，全局阈值调整容易出现“筛不掉”或“一刀切过猛”。  
优先引入“致病菌面板 + 共生菌惩罚”的选择性打分，能更稳定地区分候选优先级。

---

## 2. 面板定义（v1）

## 2.1 主要靶菌（Primary Pathogen Panel）

- Clostridium difficile (NT5083)
- Clostridium perfringens (NT5032)
- Bacteroides fragilis (ET) (NT5033)
- Fusobacterium nucleatum (NT5025)

## 2.2 扩展靶菌（Secondary Pathogen Panel）

- Bacteroides fragilis (NT) (NT5003)
- Bilophila wadsworthia (NT5036)
- Ruminococcus gnavus (NT5046)
- Escherichia coli ED1a (NT5078) [待注释确认]
- Escherichia coli IAI1 (NT5077) [待注释确认]

## 2.3 共生保留对照（Commensal Sparing Panel）

- Akkermansia muciniphila (NT5021)
- Bifidobacterium longum (NT5028)
- Bifidobacterium adolescentis (NT5022)
- Lactobacillus paracasei (NT5042)
- Eubacterium rectale (NT5009)
- Roseburia intestinalis (NT5011)
- Roseburia hominis (NT5079)
- Ruminococcus bromii (NT5045)

## 2.4 先决校验

所有菌名必须与 `workflows/reinvent4/inputs/strain_index.tsv` 完全匹配（含 NT 编号）。

---

## 3. 评分定义（建议）

设：

- `T` = 致病菌集合
- `C` = 共生菌集合
- `t` = `app_threshold`
- `tau` = 软阈值温度（默认 `0.02`）
- `sigmoid(x) = 1 / (1 + exp(-x))`

定义：

- `pathogen_soft = mean_{s in T} sigmoid((p_s - t)/tau)`
- `pathogen_hard = sum_{s in T} 1[p_s >= t]`
- `commensal_soft = mean_{s in C} sigmoid((p_s - t)/tau)`
- `selectivity_score = pathogen_soft - lambda * commensal_soft` （建议 `lambda=0.5~1.0` 可配）

排序优先级建议：

1. `selectivity_score` 降序  
2. `pathogen_hard` 降序  
3. `apscore_total` 升序（作为次级 tie-break）

---

## 4. 分阶段实施计划

### Phase A（快速可用，1-2 天）

- [ ] 新增面板配置文件  
  - Create: `workflows/reinvent4/inputs/objectives/pathogen_selective.panel.v1.json`
- [ ] 新增离线重打分脚本（读取 per-strain 概率表）  
  - Create: `scripts/pathogen_panel_rescore.py`
- [ ] 输出字段  
  - `pathogen_soft`, `pathogen_hard`, `commensal_soft`, `selectivity_score`
- [ ] 文档  
  - Modify: `docs/cli_reference.md`（新增“panel-rescore”示例）

验收：

- 输入含 `pred_id` 与 `antimicrobial_predictive_probability` 的 per-strain 结果；
- 输出可直接用于排序和阈值筛选；
- 与 `/score` 的 `single_strain`/`broad_spectrum_soft` 结果在子集场景下数值一致。

### Phase B（集成到主流程，3-5 天）

- [ ] 在 `mole screen` 增加可选参数  
  - `--panel-file`
  - `--panel-lambda`
  - `--panel-min-pathogen-hard`
  - `--panel-min-selectivity`
- [ ] process 模式支持“仅保留命中 + panel 指标列”
- [ ] manifest 写入 panel 配置快照（可复现）

验收：

- process 模式可直接输出 `hits`（带 panel 指标）；
- 重跑同配置可复现相同结果；
- 不启用 panel 时行为与当前完全一致。

### Phase C（上游口径对齐，1-2 天）

- [x] 将 `apscore_*` 主口径对齐到 `log2`
  - Modify: `src/predictor.py` (2026-05-05)
- [ ] 在 `manifest` 中记录 `apscore_log_base`
- [ ] 提供兼容脚本（旧结果 `ln -> log2` 换底）
- [ ] 更新文档说明”仅 apscore 数值尺度变化，不影响 broad-spectrum 判定”

验收：

- `ginhib_*` 与 `broad_spectrum` 在对齐前后逐条一致；
- `apscore_log2 == apscore_ln / ln(2)` 数值校验通过。

---

## 5. 测试计划

- [ ] Unit: `pathogen_panel_rescore` 数学正确性（含空值/缺菌报错）
- [ ] Unit: panel 配置校验（菌名存在性、重复名、空列表）
- [ ] Integration: `mole screen + panel` 输出字段与过滤逻辑
- [ ] Regression: 关闭 panel 时与当前版本输出一致
- [ ] Compatibility: 上游口径换底前后 `broad_spectrum` 全量一致

---

## 6. 风险与回滚

风险：

- 菌名不一致导致 panel 无法命中；
- 仅有 aggregate 结果时无法补算 panel（必须有 per-strain 概率）；
- `apscore` 换底影响历史阈值经验（但不影响排名与 broad-spectrum 判定）。

回滚：

- 保留当前默认行为（panel 功能显式开启）；
- 若换底上线后影响使用，临时切回 `log` 并通过 manifest 标记版本。

---

## 7. 交付物清单

- `workflows/reinvent4/inputs/objectives/pathogen_selective.panel.v1.json`
- `scripts/pathogen_panel_rescore.py`
- （可选）`src/panel_scoring.py`
- 文档更新：`docs/cli_reference.md`、`README_API.md`
- 回归报告：`docs/changes/<date>-pathogen-panel-upgrade.md`

