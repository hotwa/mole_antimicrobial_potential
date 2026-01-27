# Antimicrobial API (FastAPI) 接口说明文档

此文档描述 `MODE=api` 的 FastAPI 接口。MCP JSON-RPC 规范与示例请见 `README.md`。

---

## 1. API 地址

- 健康检查: `GET http://<YOUR_HOST>:8000/health`
- 预测接口: `POST http://<YOUR_HOST>:8000/predict`

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

### 本地脚本

项目内已提供 `test/test_api_request.py`，可直接运行：

```bash
API_URL=http://localhost:8000 \
python test/test_api_request.py
```

---

## 5. 说明

- 该接口与 MCP 服务共享同一个 predictor 单例，避免重复加载模型。
- 旧的 `/mcp` SSE 接口已废弃，建议使用 MCP JSON-RPC 或该 FastAPI 接口。
