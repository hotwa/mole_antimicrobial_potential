# REINVENT4 真实母核实验结果说明

这个目录保存的是本次针对真实四取代位点母核所做的全部 REINVENT4 集成验证结果。  
目标不是只做静态配置检查，而是确认下面这条链真正能跑通：

`真实母核 scaffold -> LibInvent 生成分子 -> 当前仓库 MolE 表征 -> XGBoost 抗菌预测 -> /score reward -> REINVENT4 RL`

这份说明重点回答四类问题：

1. `results/` 里目前有哪些目录，哪些可用，哪些只是排错中间产物
2. 每个目录和文件分别代表什么
3. 为什么新增了两个 macrocycle 模板
4. 为什么默认模板会让宏环母核的结果看起来像“不支持”或者“全部 0 分”

## 一、当前目录结构总览

本目录现在按用途拆成了两层：

- `presentation_guide/`
- `tests/`
- `runs/`

其中：

- `tests/` 保存已经完成的 smoke / sampling / validation 产物归档
- `runs/` 预留给后续正式长程 RL 输出
- `presentation_guide/` 保存汇报材料，不是运行产物

它们的定位不同，不要混着用。

## 二、哪些结果应该看，哪些不要看

### 建议作为正式参考的目录

请优先查看这三个目录：

- `tests/real_scaffold_sampling/`
- `tests/real_scaffold_single_strain_macrocycle_smoke/smoke_single_macro_retry/`
- `tests/real_scaffold_broad_macrocycle_smoke/smoke_broad_macro/`

原因：

- 这些目录是在修复兼容性问题之后得到的结果
- 它们已经避开了默认模板对宏环母核的错误惩罚
- 它们代表当前这套 workflow 对真实母核的有效参考结果

### 不建议作为科学结论使用的目录

下面这些目录属于排错过程：

- `tests/real_scaffold_single_strain_smoke/`
- `tests/real_scaffold_broad_smoke/`
- `tests/real_scaffold_single_strain_macrocycle_smoke/smoke_single_macro/`

这些目录有历史价值，但不要拿来解释“这个母核不行”或者“reward 没信号”，因为它们主要反映的是配置问题，不是化学问题。

## 三、每个目录的作用

### 1. `tests/smoke_rl/`

这是短程 smoke test 使用的运行时环境配置目录。

当前文件：

- `local.env`

这个文件的作用是给 smoke test 提供一个尽量小、尽量快、但仍然是完整闭环的运行环境。  
它不是正式生产配置，而是为了回答两个问题：

1. 这条链能不能跑通
2. 跑通后 reward 会不会被模板或桥接过程错误压掉

这个文件里比较重要的参数有：

- `REINVENT_DEVICE=cpu`
  强制用 CPU 做 smoke test，避免先被 GPU 环境问题阻塞
- `REINVENT_BATCH_SIZE=8`
  降低一次生成和一次打分的分子数，便于快速验证
- `REINVENT_MIN_STEPS=1`
- `REINVENT_MAX_STEPS=2`
  把 RL 限制在 2 step，只验证闭环，不追求真正收敛
- `REINVENT_MAX_SCORE=0.999`
  让 REINVENT4 保持合法配置，同时不至于因为阈值太低而立即终止
- `MOLE_REQUEST_TIMEOUT_SECONDS=300`
  给第一次加载 MolE/XGBoost 足够时间

### 2. `tests/real_scaffold_sampling/`

这个目录是最基础、最重要的先验检查目录。

它回答的问题是：

`libinvent.prior` 能不能基于这个真实母核生成有效分子？

如果这个目录都跑不出来，就没有必要继续做 RL。

当前文件：

- `sampling.toml`
- `sampling.csv`
- `sampling_single_strain_score.json`
- `sampling_broad_spectrum_score.json`
- `sampling_aggregate_prediction.json`
- `sampling_prediction_summary.tsv`
- `sampling_grid.png`

文件含义：

- `sampling.toml`
  本次 prior-only sampling 的实际配置文件
- `sampling.csv`
  LibInvent 原始采样结果，每行一个生成分子
- `sampling_single_strain_score.json`
  这些采样分子在单菌 objective 下的 `/score` 输出
- `sampling_broad_spectrum_score.json`
  这些采样分子在广谱 objective 下的 `/score` 输出
- `sampling_aggregate_prediction.json`
  这些采样分子的 `/predict aggregate` 输出
- `sampling_prediction_summary.tsv`
  本目录最值得直接看的汇总表，已经合并了：
  - 生成 SMILES
  - sampling NLL
  - 单菌 reward
  - 广谱 soft reward
  - `apscore_total`
  - `ginhib_total`
  - `broad_spectrum`
- `sampling_grid.png`
  从 `sampling.csv` 画出来的分子网格图，便于人工看结构

这次的实际意义是：

- 你的真实母核经过重新编号后，`libinvent.prior` 已经可以生成有效分子
- 这些生成分子已经被当前仓库成功送入 MolE 做抗菌预测
- 因此“生成 -> 表征 -> 打分”这条链已经被真实验证

### 3. `tests/real_scaffold_single_strain_smoke/`

这是最早的单菌 RL smoke 目录，使用的是通用模板：

- `single_strain_rl.toml.tpl`

这个目录存在的意义主要是暴露问题，不是给你做科学分析。

这里出现的问题包括：

- 初始时的桥接配置错误
- 模板中默认 `Unwanted SMARTS` 对宏环 scaffold 的惩罚

因此这里的结果不能用来证明“单菌 reward 不可用”。

### 4. `tests/real_scaffold_single_strain_macrocycle_smoke/`

这是修正模板之后的单菌 RL 目录。

你真正应该看的子目录是：

- `tests/real_scaffold_single_strain_macrocycle_smoke/smoke_single_macro_retry/`

当前文件通常包括：

- `objective.normalized.json`
- `reinvent.toml`
- `rl_run_1.csv`
- `run.summary.json`
- `scaffold.validation.json`
- `bridge_audit.jsonl`
- `stage1.chkpt`

含义：

- `objective.normalized.json`
  本次单菌目标的规范化配置
- `reinvent.toml`
  实际执行的 REINVENT4 配置
- `rl_run_1.csv`
  第 1 stage 的主要结果表
- `run.summary.json`
  本次运行的命令、输出目录、health check 等摘要
- `scaffold.validation.json`
  scaffold 校验记录
- `bridge_audit.jsonl`
  每批 `ExternalProcess -> /score` 的审计日志
- `stage1.chkpt`
  本轮 smoke run 的 checkpoint

这个目录的意义是：

- 单菌 reward 已经真正接入 REINVENT4
- 当前 prior 下，reward 数值偏小，但不是 0，也不是断开的

### 5. `tests/real_scaffold_broad_smoke/`

这是最早的广谱 RL smoke 目录，使用的是通用模板：

- `broad_spectrum_rl.toml.tpl`

和前面的 generic single-strain 一样，这个目录主要用于暴露默认模板不适合宏环 scaffold 的问题。

### 6. `tests/real_scaffold_broad_macrocycle_smoke/`

这是修正模板之后的广谱 RL 目录。

你真正应该看的子目录是：

- `tests/real_scaffold_broad_macrocycle_smoke/smoke_broad_macro/`

当前文件通常包括：

- `objective.normalized.json`
- `reinvent.toml`
- `rl_run_1.csv`
- `run.summary.json`
- `scaffold.validation.json`
- `bridge_audit.jsonl`
- `stage1.chkpt`
- `rl_grid.png`

这个目录的实际意义非常重要：

- 广谱 reward 已经在真实母核上返回连续值
- 它不是 dry-run，而是跑过真实 2-step smoke 的结果
- 从目前结果看，广谱 reward 比单菌 reward 更有信号

## 四、每个结果文件具体表示什么

下面按“最常见文件名”解释它们的含义。

### `reinvent.toml`

这是最终渲染出来、真正交给 REINVENT4 执行的配置文件。

它不是模板，而是已经把路径、objective、输出目录都填好的最终版本。

如果你想排查：

- 本次到底用了哪个 prior
- 指向了哪个 scaffold
- 用的是哪种 scoring template
- batch size / steps / device 是什么

先看这个文件。

### `objective.normalized.json`

这是本次 objective 的规范化版本。

它的作用是确保：

- `single_strain` / `broad_spectrum_soft` 被明确记录
- `tau`、`app_threshold`、`min_nkill` 被固定下来
- 如果输入里有缩写或默认值，最终都展开成稳定形式

### `scaffold.validation.json`

这是对 scaffold 的校验结果。

这里会告诉你：

- 输入 scaffold 在第几行
- attachment points 是哪些编号
- attachment point 个数是多少

对于你的母核，这个文件能证明：

- 重新编号之后确实是 `[*:0]` 到 `[*:3]`
- workflow 把它视为 4 个合法位点

### `run.summary.json`

这是运行摘要文件。

里面记录：

- 本次运行目录
- REINVENT4 根目录
- 实际执行命令
- 渲染配置路径
- scaffold 文件路径
- objective 文件路径
- 运行前 API health 检查结果

这个文件最适合回答：

`这次到底是怎么启动起来的？`

### `bridge_audit.jsonl`

这个文件很关键。

它记录的是每一批 REINVENT4 调用 `ExternalProcess` 时，桥接脚本对 `/score` 的真实交互审计。

适合排查的问题：

- `/score` 是否真的被调用了
- 返回的 reward 是多少
- 返回顺序是否和输入 SMILES 一致
- 为什么 `rl_run_1.csv` 中总分不对

如果你看到：

- `Antimicrobial Reward` 明明应该有值
- 但 `Score` 还是 0

先看这个文件，再看 `rl_run_1.csv`。

### `rl_run_1.csv`

这是 REINVENT4 最核心的结果表。

最重要的列：

- `SMILES`
  生成分子
- `Input_Scaffold`
  本轮输入给 LibInvent 的 scaffold 表示
- `R-groups`
  本轮生成的 decoration 片段
- `Score`
  REINVENT 的最终 stage score
- `Antimicrobial Reward`
  本地 `/score` 组件返回的 reward
- `Antimicrobial Reward (raw)`
  原始值
- `Unwanted SMARTS`
  自定义 alert 组件输出
- `Unwanted SMARTS (raw)`
  原始值
- `step`
  RL 的 step

### `stage1.chkpt`

这是第 1 stage 写出的 checkpoint。

用途：

- 后续继续训练时作为 agent 初始化
- 保存本次 smoke 的状态快照

### `sampling.csv`

这是 prior 直接采样的输出，不是 RL。

用途：

- 看 prior 在你的 scaffold 上到底会生成什么类型的化学空间
- 看是否生成有效分子
- 看是否明显偏向某些装饰模式

### `sampling_prediction_summary.tsv`

这是目前采样目录里最重要的文件。

因为它已经把以下信息合并起来了：

- 生成分子
- sampling NLL
- 单菌 reward
- 广谱 reward
- `apscore_total`
- `ginhib_total`
- `broad_spectrum`

如果你现在只想看：

`这批生成分子里，哪些最值得继续追踪？`

先看这个文件。

### `sampling_grid.png` / `rl_grid.png`

这是结构可视化图片。

用途：

- 快速人工看结构
- 观察是否出现明显重复 scaffold decoration 模式
- 观察 top 分子是否有共同取代特征

## 五、为什么新增了两个 macrocycle 模板

新增的两个模板是：

- `single_strain_rl_macrocycle.toml.tpl`
- `broad_spectrum_rl_macrocycle.toml.tpl`

### 原因不是 REINVENT4 不支持宏环

这一点必须说清楚：

REINVENT4 不是“不支持大环母核”。

这次实际已经验证：

- `libinvent.prior` 能读入这个宏环母核
- 能基于它生成有效分子
- 生成分子能被 MolE 成功表征
- 结果能被 `/score` 成功用于 RL

所以不是引擎不支持，而是默认模板里的过滤规则不适合这个母核系列。

### 问题出在默认模板的 `Unwanted SMARTS`

默认模板里有这一条：

- `[*;r{8-17}]`

它出现在：

- 上游 REINVENT4 默认示例 [staged_learning.toml](/home/lingyuzeng/project/REINVENT4-main/configs/staged_learning.toml)
- 本仓库原始模板 [single_strain_rl.toml.tpl](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/configs/templates/single_strain_rl.toml.tpl)
- 本仓库原始模板 [broad_spectrum_rl.toml.tpl](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/configs/templates/broad_spectrum_rl.toml.tpl)

这条规则本质上是在说：

`带有 attachment point 的 8~17 元环结构，默认视为 unwanted SMARTS`

这是一种比较保守的默认药化过滤，不是语法限制。

### 为什么你的母核会命中它

因为你的真实母核本身就是 macrocycle，而且 attachment point 就落在这个大环体系上。

所以生成分子非常容易继续保持：

- 宏环
- 环上的 attachment / decoration 连接特征

于是默认模板会频繁命中 `[*;r{8-17}]`。

### 为什么命中后总分会看起来像“不支持”

因为模板使用的是：

- `[stage.scoring] type = "geometric_mean"`

同时打分组件里包含：

- `custom_alerts`
- `ExternalProcess`

几何平均有一个重要性质：

- 只要有一个 component 是 `0`
- 最终总 `Score` 就会变成 `0`

于是会出现一种非常迷惑的现象：

- `Antimicrobial Reward` 明明不是 0
- 但 `Score` 却是 0

这会让人误以为：

- REINVENT4 不支持这个母核
- `/score` 没被调用
- 宏环生成失败

实际上都不是。  
真正的原因只是：

- `custom_alerts` 把这类宏环分子打成了 0
- 然后几何平均把总分压成了 0

## 六、为什么要继续使用新的两个模板

后续你应该优先用这两个模板：

- `single_strain_rl_macrocycle.toml.tpl`
- `broad_spectrum_rl_macrocycle.toml.tpl`

原因很简单：

- 它们保留了整体 RL 结构不变
- 保留了 `ExternalProcess` + `/score` 的打分方式
- 保留了通用的 alert 思路
- 但去掉了对当前宏环母核明显不合适的 `[*;r{8-17}]`

因此这两个模板的作用不是“放宽所有过滤”，而是：

`去掉一个对本项目母核族明显会误伤的默认过滤规则`

## 七、这两个新模板的区别是什么

这两个模板在结构上几乎一样，真正的区别不在 TOML 主体，而在你配给它们的 objective。

### `single_strain_rl_macrocycle.toml.tpl`

用途：

- 针对单个目标菌株优化
- 或者通过 `single_strain + strains + weights` 做多菌加权优化

它通常配合：

- `single_strain.akkermansia.real_scaffold.json`

或者未来你自己的多菌 objective JSON。

这个模板中的 `ExternalProcess` 最终会调用 `/score`，而 `/score` 会按：

- `single_strain`

逻辑去返回 reward。

也就是：

- 单菌时直接用目标菌概率
- 多菌时用加权平均概率

### `broad_spectrum_rl_macrocycle.toml.tpl`

用途：

- 针对广谱抗菌能力优化

它通常配合：

- `broad_spectrum.real_scaffold.json`

这个模板中的 `ExternalProcess` 也调用 `/score`，但 objective 是：

- `broad_spectrum_soft`

返回的是连续化后的广谱 reward，而不是某个单菌概率。

### 换句话说

这两个模板的主要区别是：

- `single_strain_rl_macrocycle.toml.tpl`
  用于优化特定菌或指定菌群
- `broad_spectrum_rl_macrocycle.toml.tpl`
  用于优化整体广谱能力

它们都适用于你的宏环母核。  
它们和旧模板的核心差异是：

- 去掉了 `[*;r{8-17}]`

而不是更改 REINVENT4 算法本身。

## 八、以后看结果文件时，如何判断是不是又遇到了“默认模板不支持”的假象

如果你以后看到下面这种情况：

- `Antimicrobial Reward` 有值
- 但 `Score` 长时间接近 0
- `Unwanted SMARTS` 是 0
- `Unwanted SMARTS (raw)` 也是 0
- `matchting_patterns (Unwanted SMARTS)` 里出现某条 SMARTS 命中

那首先不要下结论说：

- 模型不支持这个母核
- scaffold 不合法
- `/score` 出错

而应该优先怀疑：

- 某条 alert 规则和当前化学空间不兼容

特别是宏环项目里，要重点排查：

- `[*;r{8-17}]`

### 最典型的排查顺序

1. 看 `rl_run_1.csv`
   重点看：
   - `Score`
   - `Antimicrobial Reward`
   - `Unwanted SMARTS`
   - `matchting_patterns (Unwanted SMARTS)`
2. 看 `bridge_audit.jsonl`
   确认 `/score` 的原始 reward 是否正常
3. 看使用的是不是 generic 模板
4. 如果是宏环 scaffold，优先换成 macrocycle 模板

## 九、这次实验里“generic 模板”和“macrocycle 模板”的实际区别

### generic 模板

文件：

- `configs/templates/single_strain_rl.toml.tpl`
- `configs/templates/broad_spectrum_rl.toml.tpl`

包含：

```toml
params.smarts = [
    "[*;r{8-17}]",
    "[#8][#8]",
    "[#6;+]",
    "[#16][#16]"
]
```

### macrocycle 模板

文件：

- `configs/templates/single_strain_rl_macrocycle.toml.tpl`
- `configs/templates/broad_spectrum_rl_macrocycle.toml.tpl`

只保留：

```toml
params.smarts = [
    "[#8][#8]",
    "[#6;+]",
    "[#16][#16]"
]
```

所以本质差异只有一条：

- macrocycle 模板移除了 `[*;r{8-17}]`

但这条差异对你的母核是决定性的。

## 十、如何重新生成结构图

当前可视化脚本是：

- `scripts/visualize_smiles_grid.py`

例子：

```bash
pixi run python workflows/reinvent4/scripts/visualize_smiles_grid.py \
  --input-file workflows/reinvent4/results/tests/real_scaffold_sampling/sampling.csv \
  --output-file workflows/reinvent4/results/tests/real_scaffold_sampling/sampling_grid.png \
  --legend-columns NLL
```

```bash
pixi run python workflows/reinvent4/scripts/visualize_smiles_grid.py \
  --input-file workflows/reinvent4/results/tests/real_scaffold_broad_macrocycle_smoke/smoke_broad_macro/rl_run_1.csv \
  --output-file workflows/reinvent4/results/tests/real_scaffold_broad_macrocycle_smoke/smoke_broad_macro/rl_grid.png \
  --legend-columns Score "Antimicrobial Reward"
```

## 十一、当前最重要的结论

这次实验最重要的不是“已经收敛”，而是下面这些基础判断已经被验证：

1. 真实母核在重新编号后可以进入 LibInvent
2. `libinvent.prior` 可以基于它生成有效分子
3. 生成分子可以被当前仓库 MolE + XGBoost 成功预测
4. `/score` 已经能给 REINVENT4 返回连续 reward
5. generic 模板不适合这个宏环母核
6. macrocycle 模板是当前应该继续使用的模板
7. 目前广谱 reward 比单菌 reward 更有信号

后续如果你继续跑正式实验：

- 宏环母核优先使用：
  - `single_strain_rl_macrocycle.toml.tpl`
  - `broad_spectrum_rl_macrocycle.toml.tpl`
- 先从广谱优化开始更稳
- 单菌优化仍然可以做，但需要更长训练或更强初始化
