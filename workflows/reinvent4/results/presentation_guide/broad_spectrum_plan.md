# 固定真实母核的广谱抗菌优化讲解方案

这份材料对应的是：

- 固定真实母核
- 在四个 attachment points 上生成取代基
- 用 `broad_spectrum_soft` 作为 RL reward

建议把广谱部分讲成“为什么这样设计 reward，以及当前 smoke run 已经证明什么”。

## 一、讲解目标

这一部分最推荐的讲法是：

“我们不是直接让模型优化一个硬阈值的广谱标签，而是先把每个菌株的预测概率转成平滑抑菌分数，再平均成一个连续 reward，让 REINVENT4 在固定母核上持续朝广谱抗菌方向搜索取代基组合。”

这个表述有三个优点：

- 把 RL 设计说清楚了
- 把当前模型与 REINVENT4 的接口说清楚了
- 不会误导听众以为我们只是简单做了一个二分类筛选

## 二、建议的讲解顺序

### 1. 从问题定义开始

可以先讲：

- 已知一个真实母核
- 不希望打破母核骨架
- 只希望在四个指定取代位点上优化
- 目标是让分子对更多菌株更可能产生抑菌活性

这一页建议配：

- 母核 scaffold
- 四个位点编号说明

### 2. 说明为什么要先做 prior-only sampling

这里的核心信息是：

- 在做 RL 之前，必须先确认 `libinvent.prior` 能对这个母核生成有效分子
- 否则 RL 就没有可用搜索空间

建议使用：

- `../tests/real_scaffold_sampling/sampling_grid.png`
- `../tests/real_scaffold_sampling/sampling_prediction_summary.tsv`

这一页适合讲：

- 真实母核已经能产生一批带相同母核的新分子
- 这些分子已经能被 MolE/XGBoost 正常预测

### 3. 解释为什么广谱 reward 不直接用 `ginhib_total`

这里要明确说明：

- `ginhib_total` 是阈值计数，信号不够平滑
- `broad_spectrum` 是二值标签，更不适合 RL
- RL 需要连续 reward

因此当前设计使用：

`soft_i = sigmoid((p_i - app_threshold) / tau)`

`score = (1 / N) * Σ soft_i`

其中：

- `p_i` 是分子对第 `i` 个菌株的预测抑菌概率
- `app_threshold` 是当前系统已有的抑菌阈值
- `tau` 控制阈值附近的平滑程度

## 三、如何解释当前广谱结果

当前主要结果目录：

- `../tests/real_scaffold_broad_macrocycle_smoke/smoke_broad_macro/`

最值得展示的文件：

- `../tests/real_scaffold_broad_macrocycle_smoke/smoke_broad_macro/rl_grid.png`
- `../tests/real_scaffold_broad_macrocycle_smoke/smoke_broad_macro/rl_run_1.csv`

你可以这样解释：

### 1. 当前 run 的性质

这是一个 2-step smoke run，不是正式长程 RL。

它的意义是：

- 验证固定母核的广谱 reward 是否能真实闭环
- 验证 REINVENT4 是否能接收来自当前 `/score` 的连续反馈
- 验证 macrocycle 模板是否解决了默认模板误伤问题

### 2. `Score` 与 `Antimicrobial Reward` 如何讲

在当前广谱宏环 run 中，这两个值可以近似直接等同来讲。

原因是：

- `reinvent.toml` 里总分聚合方式是 `geometric_mean`
- 当前 component 只有：
  - `Antimicrobial Reward`
  - `Unwanted SMARTS`
- 本批结果里 `Unwanted SMARTS = 1.0`

所以：

`Score = geometric_mean(Antimicrobial Reward, 1.0) = Antimicrobial Reward`

你可以把它简化成一句话：

“在当前广谱宏环模板下，图里的 `Score` 就是 REINVENT4 实际看到的广谱抗菌 reward，因为这一批分子没有再被化学告警项额外惩罚。”

### 3. 当前数值应该怎么说

可以说：

- 当前 smoke run 里，广谱 reward 已经稳定出现大于 0 的连续值
- 目前示例值大约在 `0.10` 到 `0.16` 这一档
- 这说明模型已经能区分“更接近广谱抑菌目标”的分子

这里不要直接说：

- “已经找到最终最优广谱分子”
- “已经完成收敛”

更准确的表述是：

- “已经证明广谱 reward 设计是可用的”
- “已经证明固定母核的广谱 RL 闭环可以跑通”

## 四、为什么要强调 macrocycle 模板

这里建议单独放一页。

要讲清楚：

- 默认模板里有 `[*;r{8-17}]`
- 对宏环 scaffold，这条规则会把总分压成 0
- 这不是 REINVENT4 不支持宏环
- 而是默认模板把这类结构当作 unwanted SMARTS

因此当前广谱实验必须使用：

- `../../configs/templates/broad_spectrum_rl_macrocycle.toml.tpl`

这一页的结论可以是：

“如果不换 macrocycle 模板，就可能出现抗菌 reward 实际为正、但总 `Score` 被错误归零的情况，从而让 RL 朝错误方向更新，甚至让人误以为系统不支持这个母核。”

## 五、最适合放在 PPT 里的结论

可以直接写成下面这几条：

- 真实四位点母核已经可以被 LibInvent 装饰并生成有效分子
- 这些分子已经能被当前 MolE/XGBoost 体系成功打分
- 广谱 `broad_spectrum_soft` reward 已经能作为 REINVENT4 的连续优化信号
- 在宏环模板修正后，当前 smoke run 已经出现稳定非零广谱 reward
- 下一步应把当前 smoke run 扩展为更长步数的正式 RL 实验

## 六、汇报时建议避免的表达

不要直接说：

- “我们已经证明这个母核一定能优化成临床候选”
- “当前 2-step 结果已经足够说明最终活性”
- “Score 越高就等于一定广谱越强”

更准确的说法是：

- “Score 是当前 reward 设计下的优化目标”
- “它和真实实验活性之间还需要后续实验验证”
- “当前结果先证明计算闭环与 reward 设计可行”
