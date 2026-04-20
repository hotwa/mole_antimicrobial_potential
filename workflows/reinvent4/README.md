# 本仓库中的 REINVENT4 工作流说明

这个目录把当前仓库本地提供的 `POST /score` 打分接口，组织成一套可直接运行的 `REINVENT4 + LibInvent scaffold decoration` 工作流。

目标不是单纯渲染一份 REINVENT4 配置，而是让下面这条链稳定可重复：

`母核 scaffold -> LibInvent 生成分子 -> 当前仓库 MolE 表征 -> XGBoost 抗菌预测 -> /score reward -> REINVENT4 RL`

这份说明重点回答下面几类问题：

1. 这个工作流目录里有哪些文件和脚本，各自负责什么
2. 为什么 RL 必须通过 `/score`，而不是直接使用 `apscore_*`
3. 单菌、多菌、广谱三类优化目标分别怎么配置
4. 母核 scaffold 应该怎么准备
5. 为什么宏环母核要用 `*_macrocycle.toml.tpl`
6. 出现“结果像是不支持”时，应该先检查什么

## 一、核心原则

### 1. RL 真实优化目标不是 `apscore_*`

后续你用 REINVENT4 做强化学习时，真正影响 RL 成败的不是 `log` 的底，也不是 `apscore_total` 这个字段本身，而是 reward 是否满足三件事：

- 连续
- 方向正确
- 支持批量打分

当前工作流里，REINVENT4 一律通过本地 `POST /score` 获取 reward，原因就是这三个条件。

`apscore_*` 仍然保留，但它们在这里的定位是：

- 用于后分析
- 用于结果排序和人工筛选
- 用于和历史结果对比

而不是直接作为 RL 的原始 reward。

### 2. 为什么不直接把 `apscore_*` 喂给 REINVENT4

原因很直接：

- `apscore_*` 是聚合后的 log 变换分数，不是天然适合 RL 的 `[0,1]` 连续奖励
- 它的数值尺度会受 `log` 或 `log2` 的实现影响
- 它更适合做分析，不适合直接做训练驱动信号

当前 `/score` 已经把这些问题隔离掉了。  
因此 REINVENT4 端只需要消费一个稳定、连续、批量化的 `score`。

### 3. 当前支持的 RL 模式

这个工作流目前围绕 3 类目标组织：

- `single_strain`
- 多菌加权优化
- `broad_spectrum_soft`

其中多菌加权优化在接口层仍然使用 `single_strain` 模式，只是把 `strains + weights` 一起传入。

## 二、你需要准备什么

这个仓库不会 vendor REINVENT4 本体，而是接入一个外部 clone。

你需要准备：

- 一个本地可运行的 REINVENT4 目录
- 一个 `LibInvent prior` 权重文件
- 一个带 attachment points 的母核 scaffold 文件

当前机器上已经按推荐方式放好了这些路径：

- REINVENT4 外部目录：
  [REINVENT4-main](/home/lingyuzeng/project/REINVENT4-main)
- prior 建议目录：
  [reinvent4_assets/priors](/home/lingyuzeng/project/reinvent4_assets/priors)
- MolE 预训练 checkpoint 来源说明：
  [ckpt/README.md](/home/lingyuzeng/project/mole_antimicrobial_potential/ckpt/README.md)
- 当前真实母核 scaffold：
  [mother_scaffold.template.smi](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/inputs/scaffolds/mother_scaffold.template.smi)

推荐的本地环境模板在：

- [local.env.recommended.example](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/configs/local.env.recommended.example)

## 三、目录结构

### 1. `configs/`

这个目录放运行时环境模板和 REINVENT4 模板。

- [local.env.example](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/configs/local.env.example)
  通用示例环境文件
- [local.env.recommended.example](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/configs/local.env.recommended.example)
  按当前机器路径准备好的推荐版本
- [templates](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/configs/templates)
  各种 REINVENT4 staged-learning 模板
- [long_runs](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/configs/long_runs)
  长程 RL experiment spec，包含单菌基线和 A/B/C/D 分组

### 2. `inputs/objectives/`

这个目录放 `/score` 的目标配置示例。

- 单菌 objective 示例
- 多菌加权 objective 示例
- 广谱 objective 示例
- 真实母核用的单菌和广谱配置

### 3. `inputs/scaffolds/`

这个目录放 LibInvent scaffold 输入文件。

当前重要文件包括：

- [four_attachment.example.smi](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/inputs/scaffolds/four_attachment.example.smi)
- [mother_scaffold.template.smi](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/inputs/scaffolds/mother_scaffold.template.smi)

### 4. `scripts/`

关键脚本有：

- [validate_scaffold.py](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/validate_scaffold.py)
  校验 scaffold 是否满足 attachment point 规则
- [score_bridge.py](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/score_bridge.py)
  REINVENT4 `ExternalProcess` 桥接脚本，负责把分子批次送到本地 `/score`
- [run_reinvent4.py](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/run_reinvent4.py)
  渲染 TOML、校验输入、检查 API 健康状态、再启动 REINVENT4
- [run_long_rl.py](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/run_long_rl.py)
  按 chunk 方式串联多轮 RL，复用每一轮 `stage1.chkpt`，并执行外部早停监控
- [evaluate_results.py](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/evaluate_results.py)
  对 RL 输出的 top 分子做进一步 `/score` 和 `/predict` 复评估
- [visualize_smiles_grid.py](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/visualize_smiles_grid.py)
  生成分子网格图，方便人工查看结构

### 5. `results/`

这个目录是默认输出目录，保存：

- 渲染后的 `reinvent.toml`
- scaffold 校验记录
- objective 规范化结果
- bridge 审计日志
- RL 结果 CSV
- checkpoint
- 采样结果
- 可视化 PNG

具体到每个结果文件是什么意思，已经在：

- [results/README.md](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/results/README.md)

里做了详细中文说明。

## 四、Objective 配置格式

`objective` 文件本质上就是 `/score` 请求体里的 `objective` 部分。

### 1. 单菌优化

```json
{
  "mode": "single_strain",
  "strain": "Akkermansia muciniphila (NT5021)",
  "app_threshold": 0.04374140128493309,
  "min_nkill": 10,
  "tau": 0.02
}
```

这时 reward 直接来自目标菌株的预测概率。

### 2. 多菌加权优化

```json
{
  "mode": "single_strain",
  "strains": [
    "Akkermansia muciniphila (NT5021)",
    "Bacteroides thetaiotaomicron (NT5004)"
  ],
  "weights": [0.6, 0.4],
  "app_threshold": 0.04374140128493309,
  "min_nkill": 10,
  "tau": 0.02
}
```

这里的 reward 不是简单拼接，而是按归一化权重求加权平均。  
如果原始权重是 `w_i`，归一化后为：

`w_i' = w_i / Σw_i`

则多菌 reward 为：

`score = Σ(w_i' * p_i)`

其中 `p_i` 是对应菌株的抗菌预测概率。

### 3. 广谱优化

```json
{
  "mode": "broad_spectrum_soft",
  "app_threshold": 0.04374140128493309,
  "min_nkill": 10,
  "tau": 0.02
}
```

广谱模式不直接用硬阈值计数做 RL reward，而是使用平滑版的 soft inhibition：

`soft_i = sigmoid((p_i - app_threshold) / tau)`

再对 panel 求平均：

`score = (1 / N) * Σ soft_i`

这样做的原因是：

- 比 `broad_spectrum` 二值标签更连续
- 比 `ginhib_total` 更平滑
- 更适合 RL 在早期训练阶段接收信号

## 五、为什么必须先做 sampling

上游 REINVENT4 的 LibInvent 示例大多是 2 个 attachment points。  
本地校验器支持 4 个位点，但这不代表任意 prior 都一定擅长装饰四位点 scaffold。

因此，对四取代位点真实母核，sampling 不是可选项，而是前置检查。

应该先回答这个问题：

`prior 能不能基于这个 scaffold 生成足够多的有效分子？`

如果 prior-only sampling 都失败，后面的 RL 不值得继续。

当前真实母核已经通过了这一步验证，结果目录在：

- [real_scaffold_sampling](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/results/tests/real_scaffold_sampling)

## 六、Scaffold 输入规则

每个 scaffold 文件是一行一个 LibInvent scaffold。

当前工作流强制要求：

- 每个非注释行必须包含编号 attachment points，例如 `[*:0]`
- attachment point 编号必须唯一
- 编号必须从 `0` 开始连续
- 如果 RDKit 可用，attachment point 必须是 terminal 的，不允许多个 attachment point 指向同一原子

这意味着四取代位点母核是支持的，但必须写成：

- `[*:0]`
- `[*:1]`
- `[*:2]`
- `[*:3]`

你当前真实母核已经做过这个归一化处理，文件在：

- [mother_scaffold.template.smi](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/inputs/scaffolds/mother_scaffold.template.smi)

## 七、为什么新增了三个 macrocycle 模板

这个项目里新增了这三个模板：

- [single_strain_rl_macrocycle.toml.tpl](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/configs/templates/single_strain_rl_macrocycle.toml.tpl)
- [multi_strain_rl_macrocycle.toml.tpl](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/configs/templates/multi_strain_rl_macrocycle.toml.tpl)
- [broad_spectrum_rl_macrocycle.toml.tpl](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/configs/templates/broad_spectrum_rl_macrocycle.toml.tpl)

新增原因不是因为 `REINVENT4` 本身不支持宏环母核。

真实原因是：默认 generic 模板和上游示例模板里，`Unwanted SMARTS` 包含这条规则：

`[*;r{8-17}]`

它出现在这些文件里：

- [single_strain_rl.toml.tpl](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/configs/templates/single_strain_rl.toml.tpl)
- [broad_spectrum_rl.toml.tpl](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/configs/templates/broad_spectrum_rl.toml.tpl)
- [staged_learning.toml](/home/lingyuzeng/project/REINVENT4-main/configs/staged_learning.toml)

这条规则不是程序“不支持宏环”的硬编码，而是上游示例里默认加入的一条药化过滤规则。  
它会把带有较大环 attachment point 的结构当作 unwanted SMARTS。

你的真实母核正好就是宏环 scaffold，因此 generic 模板会把它误伤。

## 八、为什么 generic 模板会让结果看起来像“不支持”

这里的关键不是“生成失败”，而是“总分被打成 0”。

generic 模板和 macrocycle 模板都使用：

`type = "geometric_mean"`

这意味着总分是各个 scoring component 的几何平均。  
一旦 `custom_alerts` 命中并给出 `0`，总 `Score` 就会被压成 `0`。

因此你会看到一种容易误解的现象：

- `ExternalProcess` 返回的 `Antimicrobial Reward` 实际上是正的
- 但最终 `Score` 仍然是 `0`

这并不表示：

- REINVENT4 不支持这个母核
- LibInvent 无法生成分子
- 当前 MolE 无法预测这些分子

真正表示的是：

- 你用了不适合宏环项目的默认 template
- `custom_alerts` 在打分聚合时把总分归零了

## 九、generic 模板和 macrocycle 模板的区别

### 1. `single_strain_rl.toml.tpl` 和 `single_strain_rl_macrocycle.toml.tpl`

这两个模板的结构、学习策略、diversity filter、bridge 接法都一样。  
真正的差异只有一处：

- generic 版本包含 `[*;r{8-17}]`
- macrocycle 版本去掉了这条规则

也就是说：

- generic 版本适合默认药化过滤更保守的项目
- macrocycle 版本适合当前这种母核本身就是大环的项目

两者的 reward 目标相同，都是：

- `params.property = ${score_property}`
- 通常对应单菌 reward

### 2. `broad_spectrum_rl.toml.tpl` 和 `broad_spectrum_rl_macrocycle.toml.tpl`

这两个模板和上面的关系完全一样：

- learning strategy 一样
- diversity filter 一样
- `ExternalProcess` 桥接方式一样
- 唯一实际差异仍然是删除了 `[*;r{8-17}]`

区别在于它们使用的 `/score` 属性不同：

- `single_strain` 模板用于单菌或多菌加权 reward
- `broad_spectrum` 模板用于 `broad_spectrum_soft` reward

### 3. `multi_strain_rl.toml.tpl` 和 `multi_strain_rl_macrocycle.toml.tpl`

这组模板专门用于 A/B/C 这类多菌联合长程 RL。

它们和单菌模板的关系相同：

- generic 版本包含 `[*;r{8-17}]`
- macrocycle 版本移除了这条不适合当前大环母核的默认 alert

这意味着对当前真实宏环母核：

- 多菌联合 RL 也不应继续使用 generic `multi_strain_rl.toml.tpl`
- 应改用 `multi_strain_rl_macrocycle.toml.tpl`

因此这两个 macrocycle 模板的“功能区别”不在于是否支持宏环，而在于优化目标不同：

- `single_strain_rl_macrocycle.toml.tpl`
  面向单菌或多菌加权优化
- `broad_spectrum_rl_macrocycle.toml.tpl`
  面向广谱 soft reward 优化

## 十、什么时候应该优先用哪一个模板

如果你的母核是常规 scaffold，不是宏环，也不命中这类 ring attachment alert，可以先尝试 generic 模板：

- [single_strain_rl.toml.tpl](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/configs/templates/single_strain_rl.toml.tpl)
- [broad_spectrum_rl.toml.tpl](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/configs/templates/broad_spectrum_rl.toml.tpl)

但对你当前这个真实母核，应优先使用：

- [single_strain_rl_macrocycle.toml.tpl](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/configs/templates/single_strain_rl_macrocycle.toml.tpl)
- [multi_strain_rl_macrocycle.toml.tpl](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/configs/templates/multi_strain_rl_macrocycle.toml.tpl)
- [broad_spectrum_rl_macrocycle.toml.tpl](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/configs/templates/broad_spectrum_rl_macrocycle.toml.tpl)

原因不是“迁就模型”，而是避免一个与当前化学问题不匹配的默认 alert 把有效 reward 全部压掉。

## 十一、快速开始

### 1. 准备环境文件

优先从下面这个文件复制：

- [local.env.recommended.example](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/configs/local.env.recommended.example)

### 2. 启动本地 API

确保当前仓库的 FastAPI 服务已经启动，并且可访问：

- `POST /score`
- `POST /predict`
- `GET /health`

### 3. 校验 scaffold

```bash
pixi run python workflows/reinvent4/scripts/validate_scaffold.py \
  workflows/reinvent4/inputs/scaffolds/mother_scaffold.template.smi
```

### 4. 先做 dry-run

单菌示例：

```bash
pixi run python workflows/reinvent4/scripts/run_reinvent4.py \
  --template workflows/reinvent4/configs/templates/single_strain_rl_macrocycle.toml.tpl \
  --env-file workflows/reinvent4/configs/local.env.recommended.example \
  --objective-file workflows/reinvent4/inputs/objectives/single_strain.akkermansia.real_scaffold.json \
  --scaffold-file workflows/reinvent4/inputs/scaffolds/mother_scaffold.template.smi \
  --experiment-name real_scaffold_single_strain_demo \
  --dry-run
```

广谱示例：

```bash
pixi run python workflows/reinvent4/scripts/run_reinvent4.py \
  --template workflows/reinvent4/configs/templates/broad_spectrum_rl_macrocycle.toml.tpl \
  --env-file workflows/reinvent4/configs/local.env.recommended.example \
  --objective-file workflows/reinvent4/inputs/objectives/broad_spectrum.real_scaffold.json \
  --scaffold-file workflows/reinvent4/inputs/scaffolds/mother_scaffold.template.smi \
  --experiment-name real_scaffold_broad_demo \
  --dry-run
```

确认 dry-run 没问题后，再移除 `--dry-run` 正式启动。

## 十二、结果评估

RL 跑完后，建议再用评估脚本把 top 分子重新送进 `/score` 和 `/predict`：

```bash
pixi run python workflows/reinvent4/scripts/evaluate_results.py \
  --env-file workflows/reinvent4/configs/local.env.recommended.example \
  --objective-file workflows/reinvent4/inputs/objectives/broad_spectrum.real_scaffold.json \
  --reinvent-csv workflows/reinvent4/results/tests/real_scaffold_broad_macrocycle_smoke/smoke_broad_macro/rl_run_1.csv
```

这一步会额外产出：

- `score_response.json`
- `predict_per_strain.json`
- `predict_aggregate.json`
- `summary.tsv`

这样你就能明确区分：

- REINVENT4 训练时实际优化的是 `/score`
- 科学分析阶段再看 `apscore_*`、`ginhib_total`、`broad_spectrum`

## 十三、出现异常时优先排查什么

如果结果看起来像“不支持”或者“全 0 分”，优先按这个顺序排查：

1. 先看 `sampling.csv`
   是否 prior 本身就没有生成有效分子
2. 再看 `bridge_audit.jsonl`
   是否 `/score` 已经返回了正的 reward
3. 再看 `rl_run_1.csv`
   是否 `Antimicrobial Reward` 为正，但总 `Score` 为 0
4. 最后看 `reinvent.toml`
   是否仍然使用了 generic 模板，且保留了 `[*;r{8-17}]`

只要出现下面这种情况：

- `Antimicrobial Reward > 0`
- 但 `Score = 0`

就优先怀疑是 template 的 `custom_alerts` 问题，而不是模型能力问题。

## 十四、为什么这里使用 ExternalProcess bridge

当前工作流没有直接让 REINVENT4 调一个裸 REST 端点，而是使用：

- [score_bridge.py](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/score_bridge.py)

原因有四个：

- 可以明确保持输入分子顺序
- 可以稳定写出每批审计日志
- 可以把 `/score` 作为唯一 reward 真值源
- 即使未来 `/score` 扩展新 reward mode，REINVENT4 模板也不用大改

这比把 reward 逻辑散落到 REINVENT4 配置里更稳。
