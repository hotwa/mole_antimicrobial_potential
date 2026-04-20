# 多位点 `site_reward` 实现草案

## 目标

这个草案用于解决当前 REINVENT4 + LibInvent + MolE 组合里的一个实际问题：

- `MolE` / XGBoost 给出的 reward 是整分子抗菌分数
- 当前 RL 不知道四个位点各自的真实活性贡献
- 纯 MolE reward 会自然收敛到“单个位点主导扩展，其余位点最小 decoration”

因此，这里的目标不是去估计“位点 1/2/3/4 各自贡献多少活性”，而是定义一个**额外的结构辅助分数**，用来鼓励四个位点都真正参与修饰。

## 设计原则

本草案固定遵守下面 4 条原则：

1. 不改 `MolE` 原始定义  
   `MolE` 仍然只负责输出整分子的抗菌 reward。

2. `site_reward` 只做辅助，不喧宾夺主  
   RL 的主目标仍然是抗菌能力，不是单纯追求“大侧链”。

3. 第一版只用最稳定、最容易解释的结构量  
   首选每个位点 `R-group` 的 heavy atom count。

4. 第一版使用软分数，不使用硬阈值 gate  
   避免早期 RL 大量零分，导致训练不稳定。

## 为什么先用重原子数

当前最稳定的位点级结构信息，就是每个位点 decoration 的 heavy atom count。

原因：

- 可以直接从 `R-groups` 字段中稳定提取
- 与你手头真实化合物的经验分布一致
- 你已经观察到真实位点侧链大致落在 `4-27` 个重原子
- 这个特征不依赖额外实验真值

这意味着第一版 `site_reward` 不做位点活性归因，只回答两个问题：

- 四个位点是否都真的展开了
- 四个位点是否出现极端失衡

## 符号定义

设四个位点的 heavy atom count 分别为：

```text
h_1, h_2, h_3, h_4
```

其中：

- `h_j`：第 `j` 个位点对应 `R-group` 的重原子数
- `j ∈ {1,2,3,4}`
- heavy atom 不包含氢

设：

- 下限：`L = 4`
- 上限：`U = 27`

## 第一层：单个位点的区间分数

对每个位点定义一个软区间分数：

```text
r_j = sigmoid((h_j - L) / alpha) * sigmoid((U - h_j) / beta)
```

其中：

- `sigmoid(x) = 1 / (1 + e^(-x))`
- `alpha`：靠近下限时的平滑参数
- `beta`：靠近上限时的平滑参数

建议初始参数：

- `alpha = 1.5`
- `beta = 2.5`

这个设计的含义是：

- 当 `h_j < 4`，`r_j` 会下降，表示这个位点没有真正展开
- 当 `4 <= h_j <= 27`，`r_j` 接近 1，表示这个位点尺寸合理
- 当 `h_j > 27`，`r_j` 会再次下降，表示该位点可能过大

## 第二层：四个位点覆盖分数

定义四个位点覆盖分数：

```text
coverage = (1/4) * Σ r_j
```

解释：

- `coverage` 越高，表示越多位点同时落在合理区间
- 如果只有 1 个位点很大，其余 3 个都很小，`coverage` 不会高

这是第一版 `site_reward` 的主项。

## 第三层：位点平衡分数

仅有 `coverage` 还不够，因为：

- `[25, 4, 4, 4]`
- `[40, 1, 1, 1]`

都可能在某种程度上给出“至少有一个位点展开”的信号，但它们并不等价。

因此引入一个轻量平衡项。

先对每个位点做裁剪：

```text
h'_j = min(max(h_j, L), U)
```

然后定义：

```text
balance = exp(- std(h'_1, h'_2, h'_3, h'_4) / (mean(h'_1, h'_2, h'_3, h'_4) + eps))
```

其中：

- `std(...)`：四个位点的标准差
- `mean(...)`：四个位点的平均值
- `eps = 1e-6`

解释：

- 四个位点越均衡，`balance` 越接近 1
- 一个位点特别大、其余三个特别小，`balance` 会下降

注意：这里的平衡项只是弱辅助，不应压过抗菌主目标。

## 第一版 `site_reward`

综合上面两部分，第一版建议定义为：

```text
site_reward = 0.7 * coverage + 0.3 * balance
```

解释：

- `coverage` 是主项，优先鼓励“四个位点都参与”
- `balance` 是辅项，用来抑制“单个位点独大”的局部最优

在这个定义下，`site_reward` 的取值近似落在 `[0,1]`。

## 与 MolE reward 的组合方式

当前推荐的组合方式是线性加权：

```text
Final RL Score = lambda * MolE_reward + (1 - lambda) * site_reward
```

建议初始值：

```text
lambda = 0.85
```

即：

```text
Final RL Score = 0.85 * MolE_reward + 0.15 * site_reward
```

这样做的含义是：

- `MolE_reward` 仍然是主导目标
- `site_reward` 只是弱引导
- 不会改写 `MolE` 原始分数定义

## 为什么不直接用硬规则

不推荐一开始就使用下面这种硬 gate：

```text
site_reward = 1 if all(4 <= h_j <= 27) else 0
```

原因：

- RL 早期很容易产生大量零分
- 训练会变得非常不稳定
- 很难区分“完全失败”和“接近成功”

因此第一版必须使用软分数。

## 为什么不只用总重原子数

不推荐使用：

```text
f(h_1 + h_2 + h_3 + h_4)
```

原因：

- `[40,1,1,1]` 的总和并不低
- 但这不是你想要的四个位点共同展开

所以必须同时考虑：

- 每个位点 individually 是否展开
- 四个位点是否失衡

## 建议的实现位置

第一版实现应放在当前 REINVENT4 bridge / scoring 层，而不是改 `MolE` 模型本身。

推荐路径：

1. 训练时保留当前 `MolE_reward`
2. 从 `R-groups` 提取四个位点 decoration
3. 计算 `h_1` 到 `h_4`
4. 计算 `coverage`
5. 计算 `balance`
6. 组合得到 `site_reward`
7. 再与 `MolE_reward` 合成最终 `Final RL Score`

这样：

- `MolE` 原始打分仍可单独保存和汇报
- `site_reward` 只影响 RL 的训练方向
- 后续分析可以明确拆分“抗菌提升”和“多位点展开”

## 推荐的输出字段

如果后续进入实现，建议在 scoring 输出里额外记录：

- `mole_reward`
- `site_reward`
- `site_coverage`
- `site_balance`
- `site_heavy_atoms`
- `final_rl_score`

其中：

- `site_heavy_atoms` 建议保存成四元组，如 `[h_1, h_2, h_3, h_4]`

这样后面分析非常直接。

## 第一版验证标准

引入 `site_reward` 之后，不应只看最高分分子，还应同时比较：

1. `all_four_sites_ge_min_fraction`
2. `all_four_sites_in_reference_range_fraction`
3. top20 分子的四个位点 heavy atom 分布
4. `MolE_reward` 是否显著下降
5. `unique_smiles_ratio` 是否明显恶化

目标不是简单把侧链都做大，而是：

- 保持抗菌主目标仍然有效
- 同时让四个位点真正参与修饰

## 建议的推进顺序

建议分 3 步推进：

1. 先用当前纯 `MolE_reward` 作为 baseline  
   作为后续比较基准。

2. 引入当前草案中的弱 `site_reward`  
   先用 `lambda = 0.85` 做小规模验证。

3. 再根据结果调整  
   如果四个位点仍然明显不展开，再逐步提高 `site_reward` 权重，或进一步细化平衡项。

## 当前结论

当前最稳的工程判断是：

- 不需要改写 `MolE` 原始分数
- 但如果目标是“四个位点同时优化”，则需要额外结构辅助目标
- 第一版 `site_reward` 最适合基于 heavy atom count 实现
- 第一版应使用软分数，而不是硬 gate

这份草案的定位是“可落地的第一版设计”，不是最终最优形式。  
后续如果真实运行结果表明：

- `site_reward` 太弱
- 或者 `MolE_reward` 被干扰过多

再基于该版本做二次调整即可。
