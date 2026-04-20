# 四个位点同时优化可行性说明

## 结论

可行，但**不会在当前纯 MolE reward 配置下自动发生**。

当前 REINVENT4 + LibInvent 在四个位点宏环母核上已经证明：

- 四个位点 scaffold 输入合法
- 可以正常生成分子
- 可以正常用 MolE/XGBoost 打抗菌分数
- 可以正常进行 RL 微调

但当前 reward 只优化抗菌能力，不优化“四个位点是否同时展开”。  
因此模型会自然学到一种更便宜的策略：

- 只在 1 个主导位点长出大侧链
- 其余 3 个位点保持极小 decoration

## 当前证据

四个位点 scaffold：

- [mother_scaffold.template.smi](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/inputs/scaffolds/mother_scaffold.template.smi)

真实 RL 运行的四个位点分析结果：

- [chunk_001 summary](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/results/runs/pathogen_group_a_chunk001/real_chunk001_cpu_fix/chunk_001/rgroup_site_summary.json)
- [chunk_002 summary](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/results/runs/pathogen_group_a_chunk002/real_chunk002_cpu_fix/chunk_001/rgroup_site_summary.json)

关键结果：

- `chunk_001`: `all_four_sites_ge_min_fraction = 0.0`
- `chunk_002`: `all_four_sites_ge_min_fraction = 0.0`

说明：

- 没有任何一个分子实现四个位点同时达到 `>= 4` 重原子

最优分子的四个位点重原子数：

- `chunk_001` 最优分子：`[41, 1, 1, 1]`
- `chunk_002` 最优分子：`[1, 1, 36, 1]`

这两个结果都说明：当前高分分子主要来自单个位点主导扩展，而不是四个位点协同展开。

## 为什么当前 reward 下不会自然出现四个位点同时展开

当前 RL score 的定义是抗菌目标：

- 单菌：`score = p_target`
- 多菌：`score = Σ(w_i' * p_i)`
- 广谱：`score = (1/N) * Σ sigmoid((p_i - app_threshold) / tau)`

这些分数都不包含“每个位点都要长出有效侧链”的信息。  
因此在当前目标函数下，REINVENT4 只会寻找“最容易提高抗菌分数”的结构策略，而不是“最均匀展开四个位点”的结构策略。

## 技术上是否可行

技术上是可行的。

实现路径有三类：

1. 仅做后筛选  
   训练时仍然用纯 MolE reward；训练后只保留四个位点都满足条件的分子。

2. 加弱结构先验  
   MolE reward 仍是主项，额外加一个较弱的“多位点展开”结构信号。

3. 从生成机制层面约束  
   改 decoration 生成或过滤逻辑，减少“小 decoration 占位多个位点”的情况。

## 当前最推荐路线

当前最推荐路线不是立刻引入硬约束，而是：

1. 继续保留纯 MolE reward baseline  
2. 用离线分析脚本持续统计四个位点展开情况  
3. 如果多轮 run 后仍然只出现单个位点主导扩展，再引入弱结构先验

原因：

- 先把抗菌主目标上限摸清楚
- 先确认当前 scaffold 上哪些位点天然更容易承担主导优化
- 避免过早用结构约束干扰抗菌主目标学习

## 当前已经实现的分析工具

用于统计四个位点侧链大小的脚本：

- [analyze_rgroup_sites.py](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/scripts/analyze_rgroup_sites.py)

它会输出：

- `rgroup_site_analysis.tsv`
- `rgroup_site_summary.json`
- `rgroup_site_top_molecules.tsv`

这些文件可以直接回答：

- 每个位点的重原子数是多少
- 有几个位点真正展开了
- 高分分子到底是单个位点主导还是多位点协同

## 汇报建议

组会里建议用下面这句话总结：

“当前 REINVENT4 在固定四位点宏环 scaffold 上能够稳定进行抗菌 RL 优化，但在纯抗菌 reward 下，优化会自然集中到单个主导位点。四个位点同时展开在当前配置下尚未出现，因此如果项目目标要求多位点协同修饰，后续需要额外的结构引导信号。”  
