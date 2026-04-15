# 原始 MolE 抗菌分数与 REINVENT4 RL reward 的对齐说明

这份文件专门回答一个问题：

`当前仓库原始的 MolE / XGBoost 抗菌输出，如何换算成与 REINVENT4 一致的 RL reward？`

这份材料适合在 PPT 里单独做一页“分数系统对齐”说明。  
它的目标不是重复接口说明，而是明确：

1. 原始预测值是什么
2. 聚合指标是什么
3. RL 真正使用的 reward 又是什么
4. 如果想从原始 MolE 输出重新计算一个与 REINVENT4 一致的活性分数，应该怎么做

## 一、先把三层数值区分开

当前系统里至少有三层不同性质的数值：

### 1. 每个菌株的原始预测概率

记为：

```text
p_i = antimicrobial_predictive_probability
```

这是当前模型最原始、信息量最高的输出。  
来源是 XGBoost `predict_proba` 的正类概率列，代码见：

- [predictor.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/predictor.py#L235)
- [predictor.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/predictor.py#L255)

范围：

```text
0 <= p_i <= 1
```

解释：

- `p_i` 越大，表示该分子越可能抑制第 `i` 个菌株

### 2. 原始分析用聚合指标

这类指标主要包括：

- `growth_inhibition`
- `ginhib_total`
- `broad_spectrum`
- `apscore_total`
- `apscore_gnegative`
- `apscore_gpositive`

来源代码：

- [predictor.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/predictor.py#L177)
- [predictor.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/predictor.py#L242)

它们适合：

- 结果分析
- 结果展示
- 后处理筛选

但它们不是当前 RL 直接优化的目标。

### 3. REINVENT4 使用的 RL reward

这部分由 `POST /score` 定义，入口在：

- [api_server.py](/home/lingyuzeng/project/mole_antimicrobial_potential/api_server.py#L43)

真正计算逻辑在：

- [reinvent_scoring.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/reinvent_scoring.py#L61)
- [reinvent_scoring.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/reinvent_scoring.py#L86)

它是为 REINVENT4 专门做的 reward shaping，要求满足：

- 连续
- 方向统一为越大越好
- 支持批量计算

## 二、原始 MolE 分析指标是怎么计算出来的

### 1. 二值抑菌判断

先对每个菌株概率 `p_i` 做阈值化：

```text
growth_inhibition_i = 1[p_i >= app_threshold]
```

代码见：

- [predictor.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/predictor.py#L242)

默认阈值：

```text
app_threshold = 0.04374140128493309
```

定义见：

- [models.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/models.py#L60)

### 2. 广谱统计量

总抑菌菌株数：

```text
ginhib_total = Σ 1[p_i >= app_threshold]
```

广谱标签：

```text
broad_spectrum = 1[ginhib_total >= min_nkill]
```

代码见：

- [predictor.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/predictor.py#L197)
- [predictor.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/predictor.py#L250)

默认：

```text
min_nkill = 10
```

定义见：

- [models.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/models.py#L66)

### 3. `apscore_total`

当前仓库实现是：

```text
gmean_prob = geometric_mean(p_1, p_2, ..., p_N)
apscore_total = ln(gmean_prob)
```

代码见：

- [predictor.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/predictor.py#L183)
- [predictor.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/predictor.py#L186)

同样逻辑也在：

- [predict_api.py](/home/lingyuzeng/project/mole_antimicrobial_potential/predict_api.py#L352)

这里要特别强调方向：

因为：

```text
0 < gmean_prob <= 1
```

所以：

```text
apscore_total <= 0
```

并且：

```text
gmean_prob 越大  ->  apscore_total 越接近 0
gmean_prob 越小  ->  apscore_total 越负
```

因此在当前实现下：

- `apscore_total` 不是越低越强
- 更准确地说，是越接近 `0` 越强

## 三、为什么 RL 不直接用 `apscore_total`

这是后续汇报里必须讲清楚的点。

`apscore_total` 不直接用于 RL，原因有三个：

1. 它是 log 变换后的聚合指标，不是天然 `[0,1]` reward
2. 它会受 `log` 或 `log2` 实现口径影响
3. 它把所有菌株压成一个数以后，已经丢失了很多可用于定向优化的信息

所以当前 REINVENT4 一律不直接吃 `apscore_total`。  
真正进入 RL 的是基于每个菌株概率 `p_i` 重新计算的 reward。

## 四、如何把原始 MolE 输出重新计算成与 RL 一致的分数

### 方案 A：单菌优化

如果目标只是某一个菌株，最一致、最直接的 RL 分数就是：

```text
score_single = p_target
```

也就是目标菌株本身的预测概率。

来源代码：

- [reinvent_scoring.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/reinvent_scoring.py#L61)

当前实现里：

```text
score_single = weighted_probability
```

当 panel 里只有一个菌株时：

```text
weighted_probability = p_target
```

这也是为什么单菌 RL 分数的方向是：

- 越高越强

### 方案 B：多菌加权优化

如果要同时优化多个指定菌株，则先对权重归一化：

```text
w_i' = w_i / Σw_i
```

再计算：

```text
score_multi = Σ(w_i' * p_i)
```

来源代码：

- [models.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/models.py#L243)
- [reinvent_scoring.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/reinvent_scoring.py#L67)

它的方向同样是：

- 越高越强

### 方案 C：广谱优化

广谱模式不是直接用：

```text
ginhib_total = Σ 1[p_i >= app_threshold]
```

而是先做平滑变换：

```text
soft_i = sigmoid((p_i - app_threshold) / tau)
```

再对所有菌株求平均：

```text
score_broad = (1 / N) * Σ soft_i
```

代码见：

- [reinvent_scoring.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/reinvent_scoring.py#L86)

文档来源：

- [docs/reinvent4/README.md](/home/lingyuzeng/project/mole_antimicrobial_potential/docs/reinvent4/README.md#L158)

它的方向也统一为：

- 越高越强

## 五、如果手里只有 `apscore_total`，能否重算出 RL reward

这里要非常谨慎。

### 1. 能精确恢复什么

在当前仓库实现中：

```text
apscore_total = ln(gmean_prob)
```

所以可以精确恢复：

```text
gmean_prob = exp(apscore_total)
```

这给了你一个“方向与 RL 一致”的连续值：

- `gmean_prob` 越大越强

这个值虽然不是当前 RL 实际使用的 reward，但它比原始 `apscore_total` 更适合在 PPT 里和 RL 奖励做方向对齐。

### 2. 不能精确恢复什么

如果你手里只有：

- `apscore_total`
- `ginhib_total`
- `broad_spectrum`

而没有每个菌株的 `p_i`，那么你**不能精确重算**：

- 单菌 reward
- 多菌加权 reward
- `broad_spectrum_soft`

原因很直接：

- 单菌和多菌 reward 需要知道具体菌株的具体概率
- `broad_spectrum_soft` 也要求每个菌株的 `p_i`
- 聚合后的 `apscore_total` 已经丢掉了这些细节

因此：

**想与 REINVENT4 的 reward 完全一致，必须保留或重新请求 per-strain probabilities。**

最稳的做法是：

1. 用 `/predict` 获取 `aggregate_scores=false` 的逐菌株输出
2. 用同样公式重新计算单菌、多菌或广谱 reward
3. 或者直接调用 `/score`

## 六、推荐的对齐口径

如果你在 PPT 里想把“原始 MolE 分数”和“REINVENT4 reward”讲成一套清楚的体系，建议这样讲：

### 1. 原始模型基础层

最底层输出是每个菌株的预测概率：

```text
p_i = antimicrobial_predictive_probability
```

### 2. 原始分析展示层

从 `p_i` 可以计算：

- `growth_inhibition`
- `ginhib_total`
- `broad_spectrum`
- `apscore_total`

这一层主要用于：

- 结果分析
- 结果汇报
- 后筛选

### 3. RL 优化层

为了让 REINVENT4 更稳定优化，需要从同一批逐菌株概率 `p_i` 重新构造连续 reward。

这里：

- `p_i`：分子对第 `i` 个菌株的 `antimicrobial_predictive_probability`
- `i`：菌株索引
- `score`：最终返回给 REINVENT4 的 RL reward
- 所有 reward 的方向都统一为：越大越好

#### 单菌优化

如果只优化一个目标菌株，则：

```text
score = p_target
```

其中：

- `p_target`：目标菌株的预测抑菌概率

这表示 REINVENT4 直接优化该目标菌株的预测活性。

#### 多菌加权优化

如果同时优化多个指定菌株，则先对用户给定权重做归一化：

```text
w_i' = w_i / Σw_i
```

然后计算：

```text
score = Σ(w_i' * p_i)
```

其中：

- `p_i`：第 `i` 个目标菌株的预测抑菌概率
- `w_i`：用户指定的原始权重
- `w_i'`：归一化后的权重
- `Σ`：对所有目标菌株求和

这个 reward 的含义是：对指定病原菌集合做加权平均优化。  
权重越大，说明该菌株对最终 reward 的影响越大。

#### 广谱优化

广谱模式不直接使用二值抑菌结果，而是先对每个菌株做平滑变换：

```text
soft_i = sigmoid((p_i - app_threshold) / tau)
```

其中：

- `app_threshold`：当前系统用于定义抑菌的阈值
- `tau`：平滑温度参数，控制阈值附近 reward 变化的陡峭程度
- `sigmoid(x) = 1 / (1 + e^(-x))`

然后对整个 panel 求平均：

```text
score = (1 / N) * Σ soft_i
```

其中：

- `N`：参与广谱评估的菌株总数
- `soft_i`：第 `i` 个菌株的平滑抑菌值

这个 reward 的含义是：优化“更多菌株接近或超过抑菌阈值”的程度，而不是只优化某几个指定菌株的高概率。

因此：

- 单菌 reward：优化单个菌株
- 多菌 reward：优化指定菌株集合的加权平均概率
- 广谱 reward：优化整个菌株 panel 的平滑广谱抑菌能力

## 七、如果你希望原始展示分数也和 RL 方向一致，推荐怎么做

建议不要再直接把 `apscore_total` 当作“越低越好”的总活性分数来讲。

更稳的方式是：

### 1. 保留原始 `apscore_total`

因为它和历史实现一致，便于和旧结果对照。

### 2. 额外增加一个方向统一的展示量

例如：

```text
gmean_prob = exp(apscore_total)
```

在当前实现下，`gmean_prob` 与 `apscore_total` 完全等价，但方向更适合和 RL 对齐：

- `gmean_prob` 越高越强
- `apscore_total` 越接近 0 越强

### 3. 广谱页优先展示 `score_broad`

因为它就是当前 REINVENT4 实际优化的 reward。

### 4. 单菌页优先展示目标菌株概率

因为它就是当前单菌 REINVENT4 的 reward。

## 八、推荐你在 PPT 中直接使用的结论句

可以直接用下面这段：

“当前 MolE/XGBoost 的原始输出是每个菌株的抑菌预测概率。传统的 `apscore_total`、`ginhib_total` 和 `broad_spectrum` 只是基于这些概率做的分析型聚合指标；而 REINVENT4 并不直接优化这些聚合指标，而是从同一批逐菌株概率重新构造连续 reward。单菌任务直接使用目标菌株概率，多菌任务使用归一化加权概率，广谱任务使用基于阈值和平滑参数的 soft inhibition 平均值。因此，REINVENT4 优化的是与原始模型同源、但方向统一且更适合强化学习的活性分数。” 

## 九、代码与文档来源总表

### 原始 MolE 聚合逻辑

- [predictor.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/predictor.py#L177)
- [predict_api.py](/home/lingyuzeng/project/mole_antimicrobial_potential/predict_api.py#L352)

### `/score` 接口入口

- [api_server.py](/home/lingyuzeng/project/mole_antimicrobial_potential/api_server.py#L43)

### RL reward 计算逻辑

- [reinvent_scoring.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/reinvent_scoring.py#L61)
- [reinvent_scoring.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/reinvent_scoring.py#L86)

### 文档版公式说明

- [docs/reinvent4/README.md](/home/lingyuzeng/project/mole_antimicrobial_potential/docs/reinvent4/README.md#L101)
