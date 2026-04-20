# RL 目标菌株如何选择与传参

这份材料专门回答下面几个问题：

1. 单菌优化时，应该如何把目标菌株传给 REINVENT4 / `/score`
2. 多菌联合优化时，应该如何传多个菌株和权重
3. 当前系统是否支持用数字编号代替菌株名称
4. 当前 40 个菌株的编号与名称对应表在哪里

## 一、先说结论

当前实现里：

- `/score` 内部仍然按 **菌株全名字符串** 做校验和打分
- 但 REINVENT4 objective 解析层现在已经支持数字编号别名
- 如果名称拼错，或者编号超出范围，当前代码会直接报错

对应校验逻辑见：

- [models.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/models.py#L226)
- [reinvent_scoring.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/reinvent_scoring.py#L40)

尤其是这里：

- [reinvent_scoring.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/reinvent_scoring.py#L55)

如果请求的菌株名不在当前可用面板里，会报：

```text
Unknown strains requested: [...]
```

所以现阶段有两种稳妥做法：

- objective 里直接写完整菌株名
- objective 里写 `strain_index` / `strain_indices`，由 workflow 先映射成正式菌株名

## 二、当前支持多少个菌株

当前 REINVENT4 工作流使用的菌株面板来自：

- [maier_screening_results.tsv.gz](/home/lingyuzeng/project/mole_antimicrobial_potential/data/01.prepare_training_data/maier_screening_results.tsv.gz)

读取逻辑见：

- [mcp_server_enhanced.py](/home/lingyuzeng/project/mole_antimicrobial_potential/mcp_server_enhanced.py#L26)
- [reinvent4_workflow.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/reinvent4_workflow.py#L499)

当前面板总数是：

```text
40
```

因此你现在讲 PPT 时应该说：

- 当前本仓库 REINVENT4 工作流针对的是 40 个菌株的预测面板

不要再讲“100+ 菌株”作为这条 RL workflow 的当前实际面板数。  
`100+` 是项目更宽泛的历史描述，但当前这一套 RL 接入所用训练头是 40 个菌株。

## 三、当前 40 个菌株编号与名称对照表

我已经把编号与名称对照表整理成：

- [strain_index.tsv](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/inputs/strain_index.tsv)

格式是：

```text
index	strain_name
1	Akkermansia muciniphila (NT5021)
2	Bacteroides caccae (NT5050)
...
40	Veillonella parvula (NT5017)
```

这个编号表目前的作用是：

- 方便人工查阅
- 方便做 PPT
- 方便内部讨论“第几个菌株”
- 也作为当前 objective 编号别名输入的映射来源

## 四、单菌优化应该怎么传参

当前单菌模式写法是：

```json
{
  "mode": "single_strain",
  "strain": "Akkermansia muciniphila (NT5021)",
  "app_threshold": 0.04374140128493309,
  "min_nkill": 10,
  "tau": 0.02
}
```

当前真实示例在：

- [single_strain.akkermansia.real_scaffold.json](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/inputs/objectives/single_strain.akkermansia.real_scaffold.json)

这里的 `strain` 必须是完整菌株名。

例如，如果你想优化第 16 个菌株 `Clostridium perfringens (NT5032)`，应该写：

```json
{
  "mode": "single_strain",
  "strain": "Clostridium perfringens (NT5032)",
  "app_threshold": 0.04374140128493309,
  "min_nkill": 10,
  "tau": 0.02
}
```

现在也可以写成编号别名：

```json
{
  "mode": "single_strain",
  "strain_index": 16
}
```

当前 workflow 会先把它映射成正式菌株名，再继续运行。

## 五、多个菌株联合优化应该怎么传参

当前多菌联合优化仍然使用 `mode = "single_strain"`，只是把 `strain` 换成：

- `strains`
- `weights`

格式如下：

```json
{
  "mode": "single_strain",
  "strains": [
    "Akkermansia muciniphila (NT5021)",
    "Bacteroides thetaiotaomicron (NT5004)",
    "Clostridium perfringens (NT5032)"
  ],
  "weights": [0.5, 0.2, 0.3],
  "app_threshold": 0.04374140128493309,
  "min_nkill": 10,
  "tau": 0.02
}
```

当前示例文件：

- [multi_strain.example.json](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/inputs/objectives/multi_strain.example.json)
- [multi_strain.by_indices.example.json](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/inputs/objectives/multi_strain.by_indices.example.json)

权重会先归一化，代码见：

- [models.py](/home/lingyuzeng/project/mole_antimicrobial_potential/src/models.py#L243)

然后按下面公式计算 reward：

```text
w_i' = w_i / Σw_i
score = Σ(w_i' * p_i)
```

也可以使用编号别名：

```json
{
  "mode": "single_strain",
  "strain_indices": [1, 6, 16],
  "weights": [0.5, 0.2, 0.3]
}
```

当前 workflow 会先解析为：

- `Akkermansia muciniphila (NT5021)`
- `Bacteroides thetaiotaomicron (NT5004)`
- `Clostridium perfringens (NT5032)`

## 六、为什么底层仍然使用菌株名称而不是数字编号

因为当前实现设计上把“菌株名称”当作唯一标识，而不是“编号”。

原因有两个：

1. 名称更稳定、更可读
2. 编号本质上只是当前一个面板文件里的显示顺序，不是模型内部的全局 ID

当前数据头里的菌株顺序来自：

- [maier_screening_results.tsv.gz](/home/lingyuzeng/project/mole_antimicrobial_potential/data/01.prepare_training_data/maier_screening_results.tsv.gz)

如果未来数据头顺序调整，编号就可能变，但菌株名称本身不变。

所以从工程稳定性看：

- 名称传参比数字编号更安全

## 七、当前最推荐的实际使用方式

现在最推荐你这样做：

### 1. 内部讨论和 PPT

使用：

- `strain_index.tsv` 的编号表

便于说：

- “我们现在重点看第 1、6、16 号菌株”

### 2. 真正跑 REINVENT4 objective

可以用两种方式：

- 菌株全名字符串
- 编号别名

例如：

```json
{
  "mode": "single_strain",
  "strains": [
    "Akkermansia muciniphila (NT5021)",
    "Bacteroides thetaiotaomicron (NT5004)",
    "Clostridium perfringens (NT5032)"
  ],
  "weights": [0.5, 0.2, 0.3]
}
```

## 八、PPT 里最稳的表述

你可以直接这样讲：

“当前 REINVENT4 objective 在用户输入层已经支持按 1 到 40 的编号选择目标菌株，但在系统内部仍然统一映射成完整菌株名称进行校验和打分。这样既保留了汇报和配置时的简洁性，又避免了因编号顺序变化直接污染底层 scoring 逻辑。” 
