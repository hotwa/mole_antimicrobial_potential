# REINVENT4 测试与验证产物归档

这个目录保存历史 smoke run、sampling、验证性 RL 运行和桥接调试产物。

目录用途：

- 保留已经产生过的真实证据链
- 避免把大体积生成物继续堆在 `results/` 根目录
- 给 `presentation_guide/` 和 `results/README.md` 提供可追溯的原始材料来源

约定：

- `results/tests/` 下的生成型文件默认不再跟踪到 Git
- 这里只保留说明文档被跟踪
- 后续正式长程 RL 结果统一输出到 `../runs/`

当前已归档的主要目录：

- `real_scaffold_sampling/`
- `real_scaffold_single_strain_smoke/`
- `real_scaffold_single_strain_macrocycle_smoke/`
- `real_scaffold_broad_smoke/`
- `real_scaffold_broad_macrocycle_smoke/`
- `smoke_rl/`
- `validation_alias_rl/`
- `validation_multi_strain_indices/`
