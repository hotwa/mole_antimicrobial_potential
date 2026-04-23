# Antimicrobial API (FastAPI) 接口说明文档

此文档描述 `MODE=api` 的 FastAPI 接口。MCP JSON-RPC 规范与示例请见 `README.md`。

---

## 1. API 地址

- 健康检查: `GET http://<YOUR_HOST>:8000/health`
- 预测接口: `POST http://<YOUR_HOST>:8000/predict`
- REINVENT4 打分接口: `POST http://<YOUR_HOST>:8000/score`

---

## 2. 请求参数格式（JSON）

请求体直接使用 `MoleculeInput` 格式：

```json
{
  "molecules": [
    {"smiles": "CCO", "chem_id": "ethanol"},
    {"smiles": "CCN", "chem_id": "ethylamine"}
  ],
  "smiles": "CCO",
  "chem_id": "ethanol",
  "aggregate_scores": true,
  "app_threshold": 0.04374,
  "min_nkill": 10
}
```

字段说明：

- `molecules`: 结构化输入（优先级最高）
- `smiles`: 单个或列表
- `chem_id`: 与 smiles 对齐的单个或列表
- `aggregate_scores`: 是否返回聚合结果
- `app_threshold`: growth_inhibition 的阈值
- `min_nkill`: broad_spectrum 判定阈值

`/score` 使用专门的请求体，详细设计与 reward 公式见
[docs/reinvent4/README.md](/home/lingyuzeng/project/mole_antimicrobial_potential/docs/reinvent4/README.md)。
REINVENT4 LibInvent 的运行脚本和模板见
[workflows/reinvent4/README.md](/home/lingyuzeng/project/mole_antimicrobial_potential/workflows/reinvent4/README.md)。
如果你在新的 NVIDIA 机器上部署，请先执行 `pixi install`
并运行 `pixi run install-cuda-torch`，必要时可通过
`MOLE_TORCH_CUDA_TAG` 和 `MOLE_TORCH_VERSION` 切换到对应 CUDA wheel。

如果你需要查找“功能 -> Python 文件 -> 推荐用法”的映射，请看：

- [docs/cli_reference.md](/home/lingyuzeng/project/mole_antimicrobial_potential/docs/cli_reference.md)
- [docs/repo_layout.md](/home/lingyuzeng/project/mole_antimicrobial_potential/docs/repo_layout.md)
- [docs/batch_screening_input_format.md](/home/lingyuzeng/project/mole_antimicrobial_potential/docs/batch_screening_input_format.md)

### API 与 CLI 调优参数的关系

FastAPI `/predict` 和 `/score` 不暴露 `mole predict` / `mole screen` 的 CLI
调优参数，比如：

- `--num-graph-workers`
- `--graph-batch-size`
- `--prefetch-batches`
- `--classifier-backend`
- `--max-batch-size`
- `--target-gpu-memory-fraction`
- `--profiling`

这些参数由服务启动时的运行环境和共享 service/scheduler 配置决定。
也就是说，API 请求体不会逐次覆盖这些 CLI tuning knobs，但服务会继承：

- `MOLE_CLASSIFIER_BACKEND`
- `MOLE_MOLE_MODEL_PATH`
- `MOLE_PICKLE_MODEL_PATH`
- `MOLE_TIMBER_MODEL_DIR`
- `CUDA_VISIBLE_DEVICES`

安装层面的 CUDA wheel 选择仍由以下环境变量控制：

- `MOLE_TORCH_VERSION`
- `MOLE_TORCH_CUDA_TAG`
- `MOLE_TORCH_INDEX_URL`

如果你需要完整的参数解释、哪些参数会影响结果内容、哪些只影响性能/调度，
请以 [README.md](/home/lingyuzeng/project/mole_antimicrobial_potential/README.md)
和
[docs/cli_reference.md](/home/lingyuzeng/project/mole_antimicrobial_potential/docs/cli_reference.md)
、
[docs/batch_screening_input_format.md](/home/lingyuzeng/project/mole_antimicrobial_potential/docs/batch_screening_input_format.md)
为准。

---

## 3. 返回格式

```json
{
  "mode": "aggregate" | "per_strain",
  "items": [ ... ]
}
```

- `aggregate` 模式: 含 `chem_id`, `apscore_total`, `ginhib_total`, `broad_spectrum` 等
- `per_strain` 模式: 含 `pred_id`, `antimicrobial_predictive_probability`, `growth_inhibition`
- `/score` 模式: 返回 REINVENT4 可直接消费的 `score ∈ [0,1]` 及调试细节

---

## 4. 示例

### curl

```bash
curl -sS http://localhost:8000/predict \
  -H "Content-Type: application/json" \
  -d '{"smiles":"CCO","aggregate_scores":true}'
```

### Python (httpx)

```python
import httpx

payload = {"smiles": "CCO", "aggregate_scores": True}
resp = httpx.post("http://localhost:8000/predict", json=payload, timeout=120.0)
print(resp.json())
```

### REINVENT4 score

```python
import httpx

payload = {
    "smiles": ["CCO", "CCN"],
    "chem_id": ["mol1", "mol2"],
    "objective": {
        "mode": "single_strain",
        "strain": "Akkermansia muciniphila (NT5021)"
    }
}
resp = httpx.post("http://localhost:8000/score", json=payload, timeout=120.0)
print(resp.json())
```

### 本地脚本

项目内已提供 `test/test_api_request.py`，可直接运行：

```bash
API_URL=http://localhost:8000 \
python test/test_api_request.py
```

---

## 5. 说明

- 该接口与 MCP 服务共享同一个 predictor 和 adaptive batching scheduler 单例，避免重复加载模型，并在后台自动进行 OOM 降级重试。
- 旧的 `/mcp` SSE 接口已废弃，建议使用 MCP JSON-RPC 或该 FastAPI 接口。
- 对 REINVENT4，推荐调用 `/score` 而不是直接拿 `apscore_*` 做 RL reward。
- 如果 `site_reward.enabled=true`，则 `/score` 请求里的 `objective.site_reward`
  必须提供 `scaffold_smiles`。仓库中的 `mole` CLI 和 REINVENT4 bridge 会在
  只提供 `--scaffold-file` 时自动补齐该字段，默认使用
  `workflows/reinvent4/inputs/scaffolds/mother_scaffold.template.smi`。
