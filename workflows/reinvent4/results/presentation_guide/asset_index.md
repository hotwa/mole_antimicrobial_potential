# 汇报材料证据索引

这个文件列出本次汇报最值得直接引用的结果文件，以及每个文件适合说明什么问题。

## 一、输入与配置

### 0. 分数系统对齐说明

- `reward_alignment.md`

适合说明：

- 原始 MolE 分析分数与 REINVENT4 reward 的区别
- 如何从逐菌株概率重新计算单菌、多菌、广谱 RL 分数
- 为什么 `apscore_total` 不直接等于 RL reward

### 0.1 菌株选择与编号表

- `strain_selection.md`
- `../../inputs/strain_index.tsv`

适合说明：

- 当前 objective 为什么按菌株全名传参
- 单菌、多菌目标如何实际写入 objective JSON
- 1 到 40 的菌株编号与名称对应关系

### 0.2 长程 RL 策略

- `long_run_rl_strategy.md`

适合说明：

- 为什么应该按 A/B/C/D 组逐步扩大 panel
- 停止长程 RL 时不该只看单个高分分子
- `0.8` 阈值对哪些任务合理，对哪些任务不合理
- `stage1.chkpt`、`prior_file`、`agent_file` 的关系

### 1. 真实母核 scaffold

- `../../inputs/scaffolds/mother_scaffold.template.smi`

适合说明：

- 我们不是从头生成任意分子，而是固定一个真实母核
- REINVENT4 只在四个 attachment points 上做取代基装饰

### 2. 广谱 objective

- `../../inputs/objectives/broad_spectrum.real_scaffold.json`

适合说明：

- 当前广谱优化不是硬阈值分类，而是平滑的 `broad_spectrum_soft`
- `app_threshold` 和 `tau` 是 reward 设计里的关键参数

### 3. 单菌 objective

- `../../inputs/objectives/single_strain.akkermansia.real_scaffold.json`

适合说明：

- 单菌优化时，reward 直接对应目标菌株概率

## 二、先验采样结果

### 1. prior-only sampling 配置

- `../tests/real_scaffold_sampling/sampling.toml`

适合说明：

- 在进入 RL 之前，先用 prior-only sampling 验证母核是否可被装饰

### 2. prior-only sampling 主表

- `../tests/real_scaffold_sampling/sampling.csv`

适合说明：

- `libinvent.prior` 已经能在这个母核上生成新分子

### 3. prior-only sampling 结构图

- `../tests/real_scaffold_sampling/sampling_grid.png`

适合说明：

- 最直观展示“固定母核 + 不同取代基”的生成结果

### 4. prior-only sampling 汇总表

- `../tests/real_scaffold_sampling/sampling_prediction_summary.tsv`

适合说明：

- 先验生成分子已经被 MolE/XGBoost 成功打分
- 可以用这个表挑几条分子做定性展示

当前可直接提到的几个例子：

- `sample2`
  `broad_spectrum_soft_score = 0.5593713187599084`
  `ginhib_total = 19`
- `sample3`
  `broad_spectrum_soft_score = 0.33980859205800173`
  `ginhib_total = 11`

这两个分子适合放在“先验证母核可行性”的页里。

## 三、广谱 RL 结果

### 1. 广谱宏环模板实际配置

- `../tests/real_scaffold_broad_macrocycle_smoke/smoke_broad_macro/reinvent.toml`

适合说明：

- 当前 run 使用的是 macrocycle 模板
- scoring 由 `Antimicrobial Reward` 和 `Unwanted SMARTS` 两部分组成

### 2. 广谱 RL 主表

- `../tests/real_scaffold_broad_macrocycle_smoke/smoke_broad_macro/rl_run_1.csv`

适合说明：

- REINVENT4 每步实际优化的分子、`Score`、`Antimicrobial Reward`
- 当前 smoke run 已经出现稳定非零广谱 reward

### 3. 广谱 RL 分子图

- `../tests/real_scaffold_broad_macrocycle_smoke/smoke_broad_macro/rl_grid.png`

适合说明：

- RL 输出的代表性分子结构
- 适合配合 `rl_run_1.csv` 讲图中分子对应的 reward

### 4. 广谱 bridge 审计日志

- `../tests/real_scaffold_broad_macrocycle_smoke/smoke_broad_macro/bridge_audit.jsonl`

适合说明：

- 每批分子确实被送到了本地 `/score`
- 不是假数据，也不是静态模拟

## 四、单菌 RL 结果

### 1. 单菌宏环模板实际配置

- `../tests/real_scaffold_single_strain_macrocycle_smoke/smoke_single_macro_retry/reinvent.toml`

适合说明：

- 单菌任务也走同一条 REINVENT4 + `/score` 闭环

### 2. 单菌 RL 主表

- `../tests/real_scaffold_single_strain_macrocycle_smoke/smoke_single_macro_retry/rl_run_1.csv`

适合说明：

- 单菌 reward 已经打通，但当前数值仍然偏弱

### 3. 单菌 bridge 审计日志

- `../tests/real_scaffold_single_strain_macrocycle_smoke/smoke_single_macro_retry/bridge_audit.jsonl`

适合说明：

- 单菌 reward 也是通过真实 `/score` 回传

## 五、解释 Score 与 Antimicrobial Reward 时要用的证据

### 1. 广谱解释

- `../tests/real_scaffold_broad_macrocycle_smoke/smoke_broad_macro/rl_run_1.csv`

当前前几行中：

- `Score = Antimicrobial Reward`
- `Unwanted SMARTS = 1.0000000`

这说明当前这批分子没有被化学告警项惩罚，因此：

`Score = geometric_mean(Antimicrobial Reward, 1.0) = Antimicrobial Reward`

### 2. 单菌解释

- `../tests/real_scaffold_single_strain_macrocycle_smoke/smoke_single_macro_retry/rl_run_1.csv`

这里可以额外展示一个现象：

- 有些行 `Unwanted SMARTS = 0.0000000`
- 这时总 `Score` 会被压低甚至归零

这能帮助解释为什么：

- `Antimicrobial Reward` 代表纯抗菌目标
- `Score` 代表 REINVENT4 最终实际用于更新的总分

## 六、不要优先拿来当主证据的目录

下面这些目录更适合作为排错记录，不建议放到正式汇报主线里：

- `../tests/real_scaffold_single_strain_smoke/`
- `../tests/real_scaffold_broad_smoke/`
- `../tests/real_scaffold_single_strain_macrocycle_smoke/smoke_single_macro/`

原因：

- 它们包含早期模板或桥接问题
- 更适合解释“为什么后来要换 macrocycle 模板”
- 不适合作为“当前正式结果”的主要证据
