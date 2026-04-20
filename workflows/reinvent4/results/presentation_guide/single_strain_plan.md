# 固定真实母核的单菌活性优化讲解方案

这份材料对应的是：

- 固定真实母核
- 在四个 attachment points 上生成取代基
- 用目标菌株概率作为 RL reward

与广谱方案相比，单菌方案更强调“靶向某一菌株或少数菌株”的定向优化能力。

## 一、建议的核心表述

最稳的讲法是：

“在单菌任务里，我们让 REINVENT4 不再关注全菌株平均表现，而是直接接收目标菌株的预测抑菌概率作为 reward，从而把取代基搜索定向到某一特定菌株的活性提升上。”

## 二、单菌 reward 设计怎么讲

当前 objective 文件：

- `../../inputs/objectives/single_strain.akkermansia.real_scaffold.json`

如果只针对一个菌株，则：

`score = p_target`

其中 `p_target` 是目标菌株的抗菌预测概率。

如果后续要讲多菌加权版，可以在这一页补充：

`score = Σ(w_i' * p_i)`

其中 `w_i'` 是归一化后的权重。

这样可以自然引出：

- 单菌优化
- 多菌加权优化
- 广谱优化

三者其实只是 reward 设计不同，而生成框架是同一套。

## 三、如何介绍当前单菌实验结果

当前主要目录：

- `../tests/real_scaffold_single_strain_macrocycle_smoke/smoke_single_macro_retry/`

最值得展示的文件：

- `../tests/real_scaffold_single_strain_macrocycle_smoke/smoke_single_macro_retry/rl_run_1.csv`
- `../tests/real_scaffold_single_strain_macrocycle_smoke/smoke_single_macro_retry/bridge_audit.jsonl`

当前应该怎么讲：

- 单菌 RL 闭环已经真实跑通
- REINVENT4 已经能接收目标菌株的 reward
- 当前 reward 数值仍明显偏小，说明在当前 prior 和当前目标菌株组合下，信号还比较弱

这时你可以把它讲成：

“单菌方向已经不是接口问题，而是优化效率和 reward 强度问题。”

## 四、为什么当前单菌信号比广谱弱

这里要讲得准确一些，不要只说“结果不好”。

更合理的解释是：

- 广谱 reward 由多菌株平滑平均构成，更容易得到非零连续信号
- 单菌 reward 完全由一个目标菌株概率决定，稀疏性更高
- 当前 prior 在这个母核上的初始分子分布，未必天然偏向该目标菌株

因此当前看到：

- 广谱 reward 大约能到 `0.10` 到 `0.16`
- 单菌 reward 在当前 smoke 里仍处在很低量级

这不表示单菌路线错误，而表示后续需要：

- 更长步数训练
- 更合适的目标菌株选择
- 更针对性的 prior / checkpoint 初始化
- 可能的 reward rescale

## 五、`Score` 与 `Antimicrobial Reward` 在单菌里怎么讲

原则与广谱相同：

- `Antimicrobial Reward` 代表纯抗菌目标分数
- `Score` 代表 REINVENT4 最终用于训练的总分

不同之处在于，单菌 run 里更容易出现：

- `Antimicrobial Reward` 很小
- 甚至某些分子还会被 `Unwanted SMARTS` 项进一步压低

因此讲解时要提醒听众：

- 不要只看单菌 `Score` 的绝对大小
- 要先看 reward 是否已经连通、是否能区分不同分子

当前阶段最重要的结论不是“已经得到高分子”，而是：

- `single_strain` reward 已经能够被当前工作流真实调用

## 六、为什么单菌部分更适合作为“下一步重点”

你在 PPT 里可以把单菌部分放在广谱后面，作为“更有挑战但更有价值”的下一阶段。

推荐口径：

- 广谱优化先证明了宏环 scaffold + RL + `/score` 闭环是成立的
- 单菌优化则是下一步更细粒度、更高难度的方向
- 一旦单菌 reward 稳定起来，就能支持更精准的菌株定向分子设计

## 七、适合放在 PPT 里的结论

可以直接写成下面几条：

- 单菌优化与广谱优化使用同一套生成框架，只是 reward 定义不同
- 当前固定母核的单菌 RL 闭环已经跑通
- 当前 reward 偏弱，说明后续重点应放在训练步数、初始化策略和目标菌株选择上
- 单菌路线的价值在于支持菌株定向优化，而不是简单提升全局平均活性

## 八、汇报时建议避免的表达

不要直接说：

- “单菌路线已经失败”
- “当前单菌分数低，所以这个母核不适合定向优化”
- “当前单菌结果弱，说明模型没用”

更准确的说法是：

- “当前单菌 reward 已经打通，但还没有进入高效优化区间”
- “当前结果说明 reward 设计可接入，下一步需要强化信号和训练强度”
