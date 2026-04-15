# REINVENT4 `/score` API 设计说明

本目录说明当前仓库新增的 `POST /score` 接口。这个接口不是替代原有的 `POST /predict`，而是专门为 REINVENT4 强化学习设计的“连续 reward 层”。

## 1. 为什么单独做 `/score`

真正影响 REINVENT4 强化学习成败的，不是 `log` 的底，而是 reward 是否满足下面三个条件：

- 连续：reward 不能过于离散，否则 RL 很难学习。
- 方向正确：分数必须满足“越大越好”，不能让模型朝相反方向优化。
- 批量化：一次请求必须能同时给一批分子打分，否则 RL 每个 epoch 会很慢。

当前仓库原本的 `/predict` 更适合预测与分析：

- `aggregate_scores=false` 时返回每个 `compound:strain` 的概率。
- `aggregate_scores=true` 时返回 `apscore_*`、`ginhib_*`、`broad_spectrum`。

这些字段对科研分析很有价值，但不适合直接原样作为 REINVENT4 reward：

- `growth_inhibition` 和 `broad_spectrum` 是阈值后二值结果，reward 太稀疏。
- `ginhib_total` 是整数计数，比二值更好，但依然是阈值驱动，不够平滑。
- `apscore_*` 是对概率做聚合后再取对数，更适合汇总与排序，不适合直接拿原值做 RL reward。

因此这里新增 `/score`，直接从每个菌株的原始预测概率构造 REINVENT4 需要的连续标量分数 `score ∈ [0, 1]`。

## 2. 一个关键原则：`/score` 不直接使用 `apscore_*`

`/score` 的 reward 直接基于每个菌株的 `antimicrobial_predictive_probability` 计算，而不是基于 `apscore_total` 反推。

这样做有三个原因：

- 避免 `log` 与 `log2` 口径差异影响 RL。
- 保证 reward 天然处于 `[0,1]`，更符合 REINVENT4 的 scoring 习惯。
- 保留最细粒度信息，方便做单菌、多菌、全菌株三类优化。

也就是说，即使未来你把上游的 `np.log` 改成 `np.log2` 合进来，`/score` 的 reward 定义也不需要改变。

## 3. `/score` 的接口设计

### 3.1 请求地址

- `POST /score`

### 3.2 请求目标

当前版本支持两种 reward 模式：

- `single_strain`
- `broad_spectrum_soft`

### 3.3 请求体

```json
{
  "smiles": ["CCO", "CCN"],
  "chem_id": ["mol1", "mol2"],
  "app_threshold": 0.04374140128493309,
  "min_nkill": 10,
  "objective": {
    "mode": "single_strain",
    "strain": "Akkermansia muciniphila (NT5021)"
  }
}
```

核心字段：

- `smiles` / `chem_id` / `molecules`
  与 `/predict` 相同，支持单个或批量输入。
- `objective.mode`
  reward 计算模式。
- `objective.strain`
  单个目标菌株的快捷写法。
- `objective.strains`
  多个目标菌株列表。
- `objective.weights`
  多菌株优化时的加权系数，长度必须与 `strains` 一致。
- `objective.tau`
  `broad_spectrum_soft` 的平滑温度参数。
- `app_threshold`
  当前模型定义抑菌的概率阈值，默认仍然是原始论文阈值。
- `min_nkill`
  仅作为响应中的参考元数据返回，不参与 `broad_spectrum_soft` 的 reward 计算。

## 4. 为什么 reward 要这样设计

这里重复强调一次，因为这决定了 REINVENT4 后续是否好用：

- 不是 `log` 底数决定 RL 成败。
- 真正关键的是 reward 是否连续、方向是否正确、接口是否批量化。

`/score` 的设计就是围绕这三件事展开的：

- 连续：不用二值 `growth_inhibition` 直接做 reward。
- 方向正确：统一返回“越大越好”的 `score`。
- 批量化：一次请求返回整批分子的 reward。

## 5. 三种优化场景

虽然当前只有两个 `mode`，但实际上已经覆盖了三种 REINVENT4 常见优化场景。

### 5.1 单菌优化

使用：

- `objective.mode = "single_strain"`
- `objective.strain = "<目标菌株>"`

公式：

```text
score = p_target
```

其中 `p_target` 是目标菌株的 `antimicrobial_predictive_probability`。

这样设计的原因：

- 概率本身就是连续值。
- 天然位于 `[0,1]`。
- 方向明确，越大表示越可能抑制该菌株。

### 5.2 多菌株优化

使用：

- `objective.mode = "single_strain"`
- `objective.strains = ["strain_a", "strain_b", ...]`
- `objective.weights = [w1, w2, ...]`

如果不提供 `weights`，默认等权。

实际使用前，接口会先把权重归一化：

```text
ŵ_i = w_i / Σw_i
```

然后 reward 为：

```text
score = Σ(ŵ_i * p_i)
```

其中 `p_i` 是第 `i` 个目标菌株的预测概率。

这样设计的原因：

- 允许你显式偏向某几个更重要的菌株。
- reward 保持连续。
- 分数仍然在 `[0,1]`。

这是当前仓库中“多菌株定向优化”最直接、最稳定的方式。

### 5.3 所有菌株优化

使用：

- `objective.mode = "broad_spectrum_soft"`
- 不传 `objective.strains`

这时接口默认把全部菌株作为 panel。

每个菌株先从概率 `p_i` 变成一个平滑抑菌值：

```text
soft_i = sigmoid((p_i - t) / tau)
```

其中：

- `t = app_threshold`
- `tau = objective.tau`

然后对整个 panel 取平均：

```text
score = (1 / N) * Σ soft_i
```

同时还会返回：

- `soft_inhibition_count = Σ soft_i`
- `soft_inhibition_ratio = score`

这样设计的原因：

- 它保留了“接近是否抑菌阈值”的信息，而不是一刀切。
- 比 `ginhib_total` 更平滑，更适合 RL。
- 比 `broad_spectrum` 二值标签信息量更高。

## 6. 为什么 `broad_spectrum_soft` 不直接用 `ginhib_total`

`ginhib_total` 的定义是：

```text
ginhib_total = Σ 1[p_i >= app_threshold]
```

它的问题是：

- 只要某个菌株概率还没过阈值，再接近也得不到额外 reward。
- 对 RL 来说，很多“快成功”的分子与“明显失败”的分子会被压成同一档。

`broad_spectrum_soft` 把这个硬阈值换成 sigmoid 平滑形式：

```text
soft_i = sigmoid((p_i - t) / tau)
```

这样有两个直接好处：

- 分子刚开始接近阈值时就能获得增量 reward。
- RL 更容易沿着“越来越广谱”的方向稳定爬升。

这也是为什么这里强调：真正影响 RL 成败的，不是 `log` 的底，而是你给 REINVENT4 的 reward 是否连续、方向是否正确、接口是否批量化。

## 7. 响应格式

```json
{
  "mode": "reinvent_score",
  "objective": {
    "mode": "single_strain",
    "strains": ["Akkermansia muciniphila (NT5021)"],
    "app_threshold": 0.04374140128493309,
    "min_nkill": 10,
    "tau": 0.02
  },
  "items": [
    {
      "chem_id": "mol1",
      "score": 0.81,
      "objective_mode": "single_strain",
      "selected_strains": ["Akkermansia muciniphila (NT5021)"],
      "selected_probabilities": {
        "Akkermansia muciniphila (NT5021)": 0.81
      },
      "weighted_probability": 0.81,
      "panel_mean_probability": 0.81,
      "panel_gmean_probability": 0.81
    }
  ]
}
```

`items[*].score` 才是推荐直接喂给 REINVENT4 的 reward。

其他字段主要用于：

- 调试 reward 是否按预期工作
- 后处理分析
- 检查多菌株加权是否正确

## 8. 三个调用示例

### 8.1 单菌优化

```json
{
  "smiles": ["CCO"],
  "chem_id": ["mol1"],
  "objective": {
    "mode": "single_strain",
    "strain": "Akkermansia muciniphila (NT5021)"
  }
}
```

### 8.2 多菌株加权优化

```json
{
  "smiles": ["CCO"],
  "chem_id": ["mol1"],
  "objective": {
    "mode": "single_strain",
    "strains": [
      "Akkermansia muciniphila (NT5021)",
      "Bacteroides fragilis (NT5004)"
    ],
    "weights": [0.7, 0.3]
  }
}
```

### 8.3 全菌株广谱优化

```json
{
  "smiles": ["CCO"],
  "chem_id": ["mol1"],
  "app_threshold": 0.04374140128493309,
  "objective": {
    "mode": "broad_spectrum_soft",
    "tau": 0.02
  }
}
```

## 9. 与原有 `apscore_*` 的关系

当前新增的 `/score` 不修改原有预测数学定义：

- `/predict` 仍然照常返回 `apscore_*`、`ginhib_*`、`broad_spectrum`
- `/score` 只是额外新增一层面向 REINVENT4 的 reward 计算

也就是说：

- 原来的科研分析流程不受影响
- 新的 REINVENT4 强化学习接口可以单独演进

如果未来需要合并上游 `log2` 提交，也只会影响 `/predict` 里 `apscore_*` 的展示尺度，不会影响 `/score` 的 reward 方向与定义。

## 10. 推荐组织方式

当前采用：

- [api_server.py](/home/lingyuzeng/project/mole_antimicrobial_potential/api_server.py)
  暴露 `POST /score`
- [src/reinvent_scoring.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/reinvent_scoring.py)
  放置纯 reward 逻辑
- [src/models.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/models.py)
  放置 `/score` 请求模型
- [test/test_score_api.py](/home/lingyuzeng/project/mole_antimicrobial_potential/test/test_score_api.py)
  放置接口测试

这个组织方式比把所有逻辑都塞进 `api_server.py` 更好，原因是：

- reward 公式可以单独测试
- 以后增加新的 REINVENT4 reward mode 更容易
- FastAPI 层和分数计算层职责分离更清楚
