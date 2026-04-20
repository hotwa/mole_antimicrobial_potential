# 基于真实母核的长程 RL 优化策略

这份材料回答下面几个关键问题：

1. 应该如何按“明确致病菌 -> 扩展致病菌 -> 全菌株广谱”逐步做长程 RL
2. 多菌联合优化的 panel 应该如何分组
3. 长程 RL 的终止指标应该怎么设计
4. `0.8` 这个阈值是否合理
5. 是否会过拟合
6. `stage1.chkpt` 是什么，权重从哪里来
7. `libinvent.prior` 到底是什么模型

这份策略不是证明“已经得到抗生素”，而是为后续**计算设计和候选分子优选**提供一个稳健的 RL 路线图。

## 一、先说结论

### 1. 当前最适合的路线

不建议一开始就直接跑“所有 40 个菌株的广谱长程 RL”。  
更稳的路线是分层推进：

- A 组：明确致病菌小面板
- B 组：扩大到更完整的致病 / pathobiont 面板
- C 组：加入“希望尽量保留的共生菌”做选择性优化
- D 组：最后再做全菌株 broad-spectrum soft 优化

### 2. 你的“连续三轮出现 0.8 以上就终止”思路是对的，但还不够稳

只看“是否出现一个 0.8 分子”太容易被偶然高分样本触发。  
更稳的做法是：

- 必须设置一个**硬上限 max_steps**
- 必须设置一个**最短训练长度 min_steps**
- 停止条件不要只看单个分子是否过线，还要看连续窗口内是否稳定重复出现

### 3. `0.8` 不是所有任务都合理

对于：

- 单菌优化
- 小规模多菌加权优化

`0.8` 是可以考虑的高阈值。

但对于：

- 全 40 菌株广谱优化

`0.8` 往往过高，容易让训练永远达不到，或者逼模型走向极端、化学空间塌缩。

### 4. `stage1.chkpt` 确实是 RL 微调后的新权重

是的。  
RL 过程中 agent 会被更新，阶段结束时写出 checkpoint。这个 checkpoint 就是你后续继续 RL 或采样的微调后权重。

### 5. 当前默认起点权重确实是 `libinvent.prior`

当前运行文件里 `prior_file` 和 `agent_file` 都指向：

- `/home/lingyuzeng/project/reinvent4_assets/priors/libinvent.prior`

所以第一轮 RL 默认是：

- `prior = libinvent.prior`
- `agent = libinvent.prior`

然后 agent 在 RL 中更新，stage 结束后保存为新的 `stage1.chkpt`。

### 6. `libinvent.prior` 是 LibInvent 的 RNN decorator 模型

不是 transformer 版。  
当前你用的是 LibInvent 的 `Decorator` 模型，底层是 LSTM。代码来源：

- [decorator.py](/home/lingyuzeng/project/REINVENT4-main/reinvent/models/libinvent/models/decorator.py)
- [model.py](/home/lingyuzeng/project/REINVENT4-main/reinvent/models/libinvent/models/model.py)

你之前运行日志里也已经显示了 `Decorator` 和 `LSTM`。

## 二、为什么要分组而不是一步到位

如果直接把所有菌株都压进同一个长程 RL 目标，会有三个问题：

1. reward 语义太宽，优化方向不够聚焦
2. 很难判断模型到底在学“抑制致病菌”还是“广泛误伤所有菌”
3. 很难做阶段性汇报和诊断

因此最合理的办法是从窄目标开始，再逐步扩大。

## 三、推荐的长程 RL 分组

下面分组是基于你给出的那份优先级建议，同时对当前 40 菌株面板做工程化落地。

### A 组：核心明确致病菌面板

这是第一组，最适合先做长程多菌联合 RL。

建议菌株：

- `Clostridium difficile (NT5083)`，编号 15
- `Clostridium perfringens (NT5032)`，编号 16
- `Bacteroides fragilis (ET) (NT5033)`，编号 3
- `Fusobacterium nucleatum (NT5025)`，编号 27

推荐用途：

- 第一轮“致病菌定向抑制”长程 RL

推荐 reward：

- 多菌加权 `single_strain + strain_indices + weights`

推荐初始权重：

```text
[0.30, 0.30, 0.20, 0.20]
```

理由：

- 两个 Clostridium 给更高权重
- ETBF 和 F. nucleatum 保留但略低

### B 组：扩展致病 / pathobiont 面板

在 A 组基础上扩大范围。

建议菌株：

- `Clostridium difficile (NT5083)`，15
- `Clostridium perfringens (NT5032)`，16
- `Bacteroides fragilis (ET) (NT5033)`，3
- `Fusobacterium nucleatum (NT5025)`，27
- `Bacteroides fragilis (NT) (NT5003)`，4
- `Bilophila wadsworthia (NT5036)`，12
- `Ruminococcus gnavus (NT5046)`，36

推荐用途：

- 第二轮“扩大病原谱”长程 RL

推荐 reward：

- 多菌加权 `single_strain + strain_indices + weights`

推荐初始权重：

```text
[0.22, 0.22, 0.16, 0.14, 0.10, 0.08, 0.08]
```

理由：

- 仍然让两株 Clostridium 和 ETBF 保持核心地位
- 其余作为扩展致病背景

### C 组：选择性优化面板

这个阶段不再只是“杀病原”，而是要考虑“少伤好菌”。

病原目标：

- `Clostridium difficile (NT5083)`，15
- `Clostridium perfringens (NT5032)`，16
- `Bacteroides fragilis (ET) (NT5033)`，3
- `Fusobacterium nucleatum (NT5025)`，27
- `Bilophila wadsworthia (NT5036)`，12

保留型对照菌：

- `Akkermansia muciniphila (NT5021)`，1
- `Bifidobacterium longum (NT5028)`，11
- `Bifidobacterium adolescentis (NT5022)`，10
- `Lactobacillus paracasei (NT5042)`，28
- `Eubacterium rectale (NT5009)`，26
- `Roseburia intestinalis (NT5011)`，34

注意：

当前仓库的 `/score` 还没有“病原加分、共生菌减分”的选择性复合 reward。  
如果要真正做选择性优化，需要新增一个 reward mode，例如：

```text
score_selective = alpha * pathogen_weighted_probability
                  - beta  * commensal_weighted_probability
```

如果暂时不改代码，C 组可以先作为：

- 分别跑病原多菌优化
- 再用共生菌 panel 做后评估过滤

### D 组：全菌株 broad-spectrum soft

最后才做：

- `mode = "broad_spectrum_soft"`
- 不传 `strains`

这时面板默认覆盖全部 40 个菌株。

推荐用途：

- 验证“固定母核是否能被推进到更全局的广谱方向”

但不建议把它作为第一组长程 RL。

## 四、我推荐的长程推进顺序

建议按下面顺序：

1. A 组长程 RL
2. B 组长程 RL
3. D 组广谱 RL
4. 如需选择性，再补 C 组 reward 扩展

原因：

- A 组最容易得到清晰信号
- B 组验证病原谱是否可扩展
- D 组验证是否能进一步转向全局广谱
- C 组最有科研价值，但需要额外 reward 设计，不适合先上

## 五、终止指标应该怎么设计

### 1. 不建议只用“出现一个 0.8 分子就停”

原因：

- REINVENT4 生成是随机的
- 单个高分子可能只是偶然样本
- 只看 top1 很容易过早停止

### 2. 也不建议“否则一直跑到无穷”

必须始终有：

- `min_steps`
- `max_steps`

因为当前 REINVENT4 默认 `simple` termination 的内置逻辑是：

- 至少跑到 `min_steps`
- 如果**平均 batch score** 超过 `max_score` 才停
- 否则一直跑到 `max_steps`

代码见：

- [terminators.py](/home/lingyuzeng/project/REINVENT4-main/reinvent/runmodes/RL/setup/terminators.py)
- [rl.md](/home/lingyuzeng/project/REINVENT4-main/contrib/reinvent-doc/tutorials/rl.md)

这里很重要：

**REINVENT4 内置终止看的是 mean batch score，不是“是否出现某个高分分子”。**

### 3. 我推荐的终止标准

对于长程 RL，我建议用“双层终止”：

#### 层 1：内置硬终止

- `min_steps = 100`
- `max_steps = 400` 或 `500`
- `max_score` 保守设置，不直接用 0.8 控死所有任务

#### 层 2：外部监控终止

每 20 或 25 steps 检查一次最近一个窗口，判断是否满足：

```text
连续 3 个窗口都满足：
1. 至少有 1 个分子 reward >= T_high
2. top5 平均 reward >= T_mean
```

这个条件比“连续三轮有一个 0.8 分子”稳得多。

### 4. 各任务推荐阈值

#### A 组核心致病菌面板

推荐：

- `T_high = 0.80`
- `T_mean = 0.60`

#### B 组扩展致病面板

推荐：

- `T_high = 0.75`
- `T_mean = 0.55`

#### D 组全菌株广谱

不建议一开始就用 `0.8`。

推荐：

- `T_high = 0.55` 到 `0.65`
- `T_mean = 0.40` 到 `0.50`

原因很直接：

- 全 40 菌株平均 soft reward 要到 `0.8` 非常苛刻
- 当前 broad macrocycle smoke 的高分也远没到这个量级
- 如果直接设 `0.8`，训练大概率只会顶到 `max_steps`

### 5. 你的“连续三轮 0.8”是否可保留

可以保留，但我建议改成：

```text
连续 3 个检查窗口中：
- 至少各出现 1 个 reward >= 0.8 的分子
- 且 top5 平均 reward >= 0.6
```

如果只保留第一条，过于脆弱。

## 六、`0.8` 合不合理

### 对单菌和小面板多菌

合理，但属于高阈值。

我们刚刚做过一个按编号别名的 3 菌多菌 RL 验证：

- 目标菌：16、34、37
- 即：
  - `Clostridium perfringens (NT5032)`
  - `Roseburia intestinalis (NT5011)`
  - `Ruminococcus torques (NT5047)`

实际 8 step 里已经出现：

- step 4：max reward `0.8802`
- step 8：max reward `0.8809`

但并没有满足“连续三轮都有 >= 0.8”的条件。  
这恰好说明：

- `0.8` 对小面板是可达的
- 但拿它做停机条件时，不能只靠一次命中

### 对全菌株 broad-spectrum

不合理，至少不适合作为首轮阈值。

## 七、是否会过拟合

会，而且长程 RL 必须默认考虑这个风险。

主要表现为：

1. reward hacking
   模型学会钻当前 reward 的空子，而不是真正变成可开发分子
2. scaffold 内部结构塌缩
   只剩少数重复取代模式
3. 误伤共生菌而不自知
   如果 reward 只奖励病原抑制，不惩罚共生菌损伤
4. 过度偏向模型盲点
   只在当前 XGBoost 判别边界上拿高分

### 如何降低过拟合风险

建议至少做这几件事：

- 保留 diversity filter
- 每 25 steps 导出一次 top 分子并复评估
- 同时看：
  - reward
  - `apscore_total`
  - `ginhib_total`
  - 共生菌单独概率
- 长程 RL 分成多个阶段，而不是一口气到底
- 后续用 checkpoint 继续跑时，不要只看 reward 最大值，还要看分子多样性

## 八、`stage1.chkpt` 是什么

它是当前 stage 结束时保存的 agent 权重文件。

来源说明：

- `prior_file` 是固定参考模型，不更新
- `agent_file` 是训练起点
- RL 更新的是 agent
- stage 结束后把 agent 写到 `chkpt_file`

官方说明见：

- [rl.md](/home/lingyuzeng/project/REINVENT4-main/contrib/reinvent-doc/tutorials/rl.md)

当前你仓库里的例子，例如：

- [stage1.chkpt](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/results/tests/real_scaffold_broad_macrocycle_smoke/smoke_broad_macro/stage1.chkpt)

就是广谱 smoke run 结束后得到的微调 agent 权重。

所以答案是：

- 是的，RL 强化过程会生成新的微调权重文件

## 九、这些微调权重从哪里来

当前默认起点是：

- `prior_file = /home/lingyuzeng/project/reinvent4_assets/priors/libinvent.prior`
- `agent_file = /home/lingyuzeng/project/reinvent4_assets/priors/libinvent.prior`

所以第一轮：

- prior 不变
- agent 从 `libinvent.prior` 初始化

RL 跑完后：

- agent 被更新并写成 `stage1.chkpt`

后续如果要继续第二轮或第二阶段，你完全可以把：

- `REINVENT4_AGENT_FILE=上一次的 stage1.chkpt`

而 `prior_file` 仍保持原始 `libinvent.prior`。

这也是 REINVENT4 官方推荐的继续训练方式。

## 十、`libinvent.prior` 是什么模型

当前使用的：

- `libinvent.prior`

对应的是 LibInvent 的 `Decorator` 模型，不是 transformer 版。

当前架构证据：

- 运行日志里显示 `Decorator`
- 日志里显示 encoder/decoder 都是 `LSTM`
- 代码实现见：
  - [decorator.py](/home/lingyuzeng/project/REINVENT4-main/reinvent/models/libinvent/models/decorator.py)
  - [model.py](/home/lingyuzeng/project/REINVENT4-main/reinvent/models/libinvent/models/model.py)

因此当前结论是：

- `libinvent.prior` 是基于 LSTM 的 LibInvent RNN decorator prior

不是：

- `libinvent_transformer_pubchem.prior`

后者是另一条路线，不是你当前默认在用的那一个。

## 十一、我推荐你现在采用的正式长程设置

### A 组核心致病菌长程 RL

- panel：15, 16, 3, 27
- weights：`[0.30, 0.30, 0.20, 0.20]`
- `min_steps = 100`
- `max_steps = 400`
- `sigma = 64` 起步
- 外部终止：
  - 连续 3 个检查窗口
  - 每个窗口至少 1 个分子 `reward >= 0.8`
  - 且 top5 平均 `>= 0.6`

### B 组扩展致病菌长程 RL

- panel：15, 16, 3, 27, 4, 12, 36
- `min_steps = 150`
- `max_steps = 500`
- `sigma = 64`
- 外部终止：
  - 连续 3 个检查窗口
  - 每个窗口至少 1 个分子 `reward >= 0.75`
  - 且 top5 平均 `>= 0.55`

### D 组全菌株 broad-spectrum 长程 RL

- mode：`broad_spectrum_soft`
- `min_steps = 150`
- `max_steps = 600`
- `sigma = 64` 或 `96`
- 外部终止：
  - 连续 3 个检查窗口
  - 每个窗口至少 1 个分子 `reward >= 0.60`
  - 且 top5 平均 `>= 0.45`

## 十二、PPT 最稳的表述

你可以直接这样讲：

“长程 RL 我们不采用一步到位的全广谱策略，而是从明确致病菌小面板开始，逐步扩展到更大的致病菌面板，最后再做全菌株 broad-spectrum soft 优化。终止条件也不只看单个高分分子是否出现，而是要求在连续多个检查窗口中稳定出现高分分子，并且 top 分子的平均 reward 维持在设定阈值以上。这样做可以降低偶然高分和 reward hacking 带来的过早停止风险。” 
