# 调研任务：优化 RDKit 分子组装吞吐量

## 背景

我们有一个高通量虚拟筛选系统，当前瓶颈在 RDKit 分子组装阶段（materialization）。

### 当前实现
```python
def _assemble_molecule(scaffold_smiles: str, fragments: Mapping[str, str]) -> str:
    scaffold = Chem.MolFromSmiles(scaffold_smiles)  # 11.7 μs
    current = scaffold
    for position_name in ("pos4", "pos3", "pos13", "pos12"):
        current = _attach_fragment(current, attachment_label, fragments[position_name])
        # 每次 _attach_fragment: 455 μs
        #   1. Chem.MolFromSmiles(fragment_smiles)
        #   2. _find_attachment_index()
        #   3. EditableMol 操作
        #   4. 返回新 Mol
    return Chem.MolToSmiles(current, canonical=True, isomericSmiles=True)  # 9.5 μs
```

### 性能 Profiling 数据（1000 次调用平均）
| 操作 | 耗时 (μs) | 调用次数/分子 | 总耗时 (μs) | 占比 |
|------|-----------|--------------|-------------|------|
| Chem.MolFromSmiles | 11.7 | 5 | 58.5 | 1.9% |
| Chem.MolToSmiles(canonical=True) | 9.5 | 1 | 9.5 | 0.3% |
| _attach_fragment | 455.1 | 4 | 1,820.4 | 57.6% |
| 其他开销 | - | - | ~1,268 | 40.2% |
| **总计** | - | - | **3,157** | 100% |

### 当前吞吐量
- 单进程: ~316 mol/s (3,157 μs/分子)
- 4 workers ProcessPoolExecutor: ~316 mol/s（未提升）
- 目标: 1000+ mol/s

### 约束
1. 必须生成 canonical SMILES（用于下游去重和预测）
2. 必须保持 global_combination_index 映射正确
3. 不能破坏分子结构正确性

## 核心瓶颈：_attach_fragment

```python
def _attach_fragment(scaffold: Chem.Mol, attachment_label: int, fragment_smiles: str) -> Chem.Mol:
    fragment = Chem.MolFromSmiles(fragment_smiles)  # ~12 μs
    scaffold_idx = _find_scaffold_attachment_index(scaffold, attachment_label)  # 遍历原子
    fragment_idx = _find_fragment_attachment_index(fragment)  # 遍历原子

    # 获取邻居原子
    scaffold_atom = scaffold.GetAtomWithIdx(scaffold_idx)
    fragment_atom = fragment.GetAtomWithIdx(fragment_idx)

    # 获取 bond type
    scaffold_bond = scaffold.GetBondBetweenAtoms(...)
    fragment_bond = fragment.GetBondBetweenAtoms(...)

    # 创建 EditableMol 进行编辑
    editable = Chem.EditableMol(scaffold)  # 复制整个分子
    editable.RemoveAtom(scaffold_idx)
    # ... 更多编辑操作

    new_mol = editable.GetMol()
    # ... 清理和验证
    return new_mol
```

## 调研问题

### 1. _attach_fragment 为什么这么慢？
- `Chem.MolFromSmiles(fragment_smiles)` 每次都重新解析 SMILES，能否预计算？
- `Chem.EditableMol(scaffold)` 每次都复制整个分子，能否避免？
- `_find_attachment_index()` 每次都遍历原子，能否缓存？

### 2. 替代实现方案
- **预计算 fragment Mol 对象**: 在初始化时预计算所有 fragment 的 Mol 对象和 attachment index
- **SMILES 模板组装**: 直接操作 SMILES 字符串而不是 RDKit Mol 对象
- **RGroupDecomposition**: 使用 RDKit 的 RGroup 工具
- **Chem.CombineMols + Chem.EditableMol**: 优化的合并方式
- **C++ 扩展**: 用 pybind11 包装自定义的 C++ 组装函数

### 3. 并行化分析
- 当前 4 workers 未提升的原因是什么？
- 如果增加 workers 数量（如 8、16），吞吐量会线性提升吗？
- 是否存在 RDKit 内部的 GIL 或锁？

### 4. 具体优化建议
请给出按优先级排序的优化方案，每个方案包含:
- 预期吞吐量提升倍数
- 实现复杂度（低/中/高）
- 是否改变输出格式
- 代码改动范围

## 参考资料

- RDKit 官方文档: https://www.rdkit.org/docs/
- RDKit 性能优化: https://rdkit.org/docs/GettingStartedInPython.html#advanced-features
- SMILES 规范: https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html
- RDKit GitHub Issues: https://github.com/rdkit/rdkit/issues

## 输出格式

请输出:
1. 瓶颈详细分析（为什么 _attach_fragment 需要 455 μs）
2. 优化方案列表（按优先级排序）
3. 每个方案的详细实现步骤和代码示例
4. 预期效果对比表
5. 推荐的下一步行动
