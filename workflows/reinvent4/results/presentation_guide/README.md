# 本次 REINVENT4 结果汇报材料总览

这个目录不是新的实验结果目录，而是为后续组会或汇报准备的讲解材料目录。  
目标是把“固定真实母核，在四个取代位点上做 REINVENT4 分子生成，再用当前 MolE/XGBoost 体系做抗菌优化”的过程拆成可直接用于 PPT 的叙事结构。

建议把汇报分成 4 个部分：

1. 研究目标与技术路线
2. 广谱抗菌优化方案与结果
3. 单菌活性优化方案与结果
4. 当前结论、限制与下一步实验

本目录下文件分工如下：

- `asset_index.md`
  汇报时要用到的原始结果文件索引，告诉你每个图表、CSV、JSON 在哪里
- `broad_spectrum_plan.md`
  如何讲“固定母核做广谱抗菌优化”的完整逻辑
- `reward_alignment.md`
  如何把原始 MolE 抗菌分数与 REINVENT4 的 RL reward 对齐，适合做公式页
- `site_reward_design.md`
  如果后续需要从“单个位点主导扩展”推进到“四个位点共同参与”，这里给出第一版弱结构辅助分数设计
- `four_site_cooptimization.md`
  如何解释“四个位点是否能同时优化”、当前为什么没有自然发生，以及后续如何推进
- `single_strain_plan.md`
  如何讲“固定母核做单菌活性优化”的完整逻辑
- `strain_selection.md`
  解释单菌、多菌目标菌株应如何传参，以及编号表如何使用
- `long_run_rl_strategy.md`
  长程 RL 的 panel 设计、终止条件、过拟合风险与 checkpoint 解释
- `README.md`
  本目录总说明

## 建议的汇报主线

最稳的讲法不是直接讲“我们已经做出高活性新分子”，而是按下面的顺序讲：

1. 先证明固定母核 scaffold 输入是合法的，REINVENT4 可以处理
2. 再证明 prior-only sampling 可以生成带这个母核的新分子
3. 再证明这些生成分子可以被当前仓库真实送入 MolE/XGBoost 做预测
4. 最后再讲 RL 如何把 reward 回传给 REINVENT4，用于继续朝广谱或单菌方向优化

这样讲的好处是：

- 不会把 smoke run 夸大成正式收敛结果
- 证据链完整
- 能清楚区分“已经打通闭环”和“已经完成优化收敛”

## 当前最稳的结论口径

当前已经可以稳定支持的说法：

- 真实母核经过 attachment point 重新编号后，已经能被 LibInvent 读取
- `libinvent.prior` 已经能在这个母核上生成新分子
- 这些分子已经能被当前 MolE/XGBoost 体系预测抗菌活性
- 本地 `/score` 已经能作为 REINVENT4 的 reward 信号工作
- 对广谱优化，当前 reward 信号明显比单菌更强
- 对单菌优化，闭环已打通，但当前 prior 下 reward 仍然偏弱

当前不建议直接下结论的说法：

- “已经收敛得到最终最优分子”
- “当前单菌策略已经优于广谱策略”
- “当前 smoke run 已经足够证明大规模 RL 有效”

## PPT 建议页序

可以按下面 8 页左右组织：

1. 研究目标：固定母核、四个位点、希望优化抗菌活性
2. 技术路线：REINVENT4 生成，MolE/XGBoost 打分，`/score` 回传 RL
3. 母核与输入规范：为什么 attachment points 必须从 `0` 连续编号
4. prior-only sampling 结果：先证明能生成带该母核的分子
5. 广谱优化 reward 设计：为什么用 `broad_spectrum_soft`
6. 广谱 smoke 结果：`Score` / `Antimicrobial Reward` 如何解释
7. 单菌优化 reward 设计与当前结果
8. 下一步：延长训练步数、调节 sigma / prior / objective

## 使用建议

真正做 PPT 时，先从 `asset_index.md` 找原始文件，再从：

- `broad_spectrum_plan.md`
- `reward_alignment.md`
- `four_site_cooptimization.md`
- `single_strain_plan.md`
- `strain_selection.md`
- `long_run_rl_strategy.md`

中摘取对应讲法。

如果你只准备讲一条主线，优先讲广谱，不要优先讲单菌。  
原因很简单：当前广谱 smoke run 已经出现了更清晰的连续 reward，而单菌 reward 还太弱，更适合作为“下一步重点优化方向”来讲。
