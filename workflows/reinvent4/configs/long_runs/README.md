# 长程 RL 实验 spec

这个目录保存给 `workflows/reinvent4/scripts/run_long_rl.py` 使用的长程实验 spec。

使用原则：

- 单菌基线使用 `single_strain` objective
- A/B/C 组使用多菌加权 `single_strain` objective
- D 组使用 `broad_spectrum_soft`
- 真实宏环母核请优先搭配：
  - `single_strain_rl_macrocycle.toml.tpl`
  - `multi_strain_rl_macrocycle.toml.tpl`
  - `broad_spectrum_rl_macrocycle.toml.tpl`

建议运行顺序：

1. 单菌基线
2. A 组核心致病菌
3. B 组扩展致病菌
4. C 组更大病原面板
5. D 组全 40 菌株 broad-spectrum
