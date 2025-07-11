# Antimicrobial MCP API 接口说明文档

本接口用于通过分子的 SMILES 表达式预测其对微生物的抗菌活性，可以输出详细分布结果，或聚合（summary）分子抗菌谱等。接口支持单分子和多分子，并支持自定义 chem\_id。

---

## 1. API 地址

POST http\://\<YOUR\_HOST>:8000/mcp

---

## 2. 请求参数格式（JSON）

### 通用格式：

```json
{
  "id": "请求ID，可自定义",
  "method": "predict",
  "params": {
    // 下列参数
  }
}
```

### params 可选参数：

* **smiles**：单个 SMILES 字符串 或 SMILES 字符串列表
* **chem\_id**：可选，单个分子的 chem\_id 或 chem\_id 列表，和 smiles 一一对应
* **aggregate\_scores**：是否返回聚合（汇总）结果（true/false）
* **app\_threshold**：二值化判定 growth\_inhibition 的阈值，默认 0.04374
* **min\_nkill**：broad\_spectrum 的判定阈值，默认 10
* **molecules**：结构化方式（列表，每个元素有 smiles 和 chem\_id 字段），优先级高于 smiles/chem\_id

#### molecules 示例

```json
"molecules": [
  {"smiles": "CCO", "chem_id": "ethanol"},
  {"smiles": "CCN", "chem_id": "ethylamine"}
]
```

---

## 3. 返回格式

返回为 Server-Sent Events (SSE) 流，每一行 event 对应一次推理步骤：

* `event: start`：任务开始
* `event: data`：结果数据
* `event: end`：任务结束
* `event: error`：错误信息

### 结果字段解释：

#### 1) 聚合（aggregate\_scores=true）：

每个分子返回：

* `chem_id`: 分子ID
* `apscore_total`：全菌种抗菌得分（对数几何均值）
* `apscore_gnegative`/`apscore_gpositive`：革兰阴/阳菌群抗菌得分
* `ginhib_total`: 被预测抑制的菌株数量
* `ginhib_gnegative`/`ginhib_gpositive`：被抑制的阴性/阳性菌株数量
* `broad_spectrum`: 是否广谱抗生素（1/0）

#### 2) 详细分布（aggregate\_scores=false）：

每个分子-菌株对返回：

* `pred_id`: 复合ID（如 "ethanol\:Akkermansia muciniphila (NT5021)"）
* `antimicrobial_predictive_probability`: 预测概率
* `growth_inhibition`: 是否抑制（1/0）

### 错误返回（event: error）:

```json
{"id": "请求ID", "error": {"message": "错误信息"}}
```

---

## 4. 用法示例（Python requests）

#### (1) 单个分子，聚合结果

```python
data = {
    "id": "case1",
    "method": "predict",
    "params": {
        "smiles": "CCO",
        "aggregate_scores": True
    }
}
```

#### (2) 多个分子，详细分布

```python
data = {
    "id": "case2",
    "method": "predict",
    "params": {
        "smiles": ["CCO", "CCN"],
        "aggregate_scores": False
    }
}
```

#### (3) 多个分子+自定义 chem\_id

```python
data = {
    "id": "case3",
    "method": "predict",
    "params": {
        "smiles": ["CCO", "CCN"],
        "chem_id": ["ethanol", "ethylamine"],
        "aggregate_scores": True
    }
}
```

#### (4) molecules 结构化用法

```python
data = {
    "id": "case4",
    "method": "predict",
    "params": {
        "molecules": [
            {"smiles": "CCO", "chem_id": "ethanol"},
            {"smiles": "CCN", "chem_id": "ethylamine"}
        ],
        "aggregate_scores": True
    }
}
```

#### (5) 单分子，指定 chem\_id

```python
data = {
    "id": "case5",
    "method": "predict",
    "params": {
        "smiles": "CCO",
        "chem_id": "ethanol",
        "aggregate_scores": False
    }
}
```

#### (6) 自定义阈值

```python
data = {
    "id": "case6",
    "method": "predict",
    "params": {
        "smiles": "CCO",
        "aggregate_scores": True,
        "app_threshold": 0.01,
        "min_nkill": 3
    }
}
```

#### 发送请求：

```python
import requests
r = requests.post('http://YOUR_HOST:8000/mcp', json=data, stream=True)
for line in r.iter_lines():
    print(line.decode())
```

---

## 5. 结果解读

* `apscore_total`/`apscore_gnegative`/`apscore_gpositive`：值越高（数值为负，对比即可），表示抗菌活性越强
* `ginhib_total`/`ginhib_gnegative`/`ginhib_gpositive`：预测抑制的菌株个数
* `broad_spectrum`: 1 表示广谱抗生素
* `antimicrobial_predictive_probability`：分子对某菌株的抑菌概率，越高越可能抑制
* `growth_inhibition`：是否达到阈值被判为抑制

---

## 6. FAQ

* 支持单分子/多分子批量预测
* chem\_id 可自定义（建议保证唯一性）
* molecules 结构优先级最高，推荐批量用
* 默认阈值来自原始论文，如需自定义可传入

---

## 7. 对 LLM/Agent 调用说明

* 每个字段均有含义描述
* 可以按需聚合或详细查看每一分子-菌株对
* molecules 字段推荐用作批量结构化传参
* 出错直接返回错误描述
* 结果均为标准 JSON 格式，适合自动化处理

---

如有问题可联系接口维护者。
