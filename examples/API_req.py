import requests

test_cases = [
    {
        "desc": "单个分子，默认chem_id",
        "data": {
            "id": "case1",
            "method": "predict",
            "params": {"smiles": "CCO", "aggregate_scores": True}
        }
    },
    {
        "desc": "多个分子，默认chem_id",
        "data": {
            "id": "case2",
            "method": "predict",
            "params": {"smiles": ["CCO", "CCN"], "aggregate_scores": False}
        }
    },
    {
        "desc": "多个分子+自定义chem_id",
        "data": {
            "id": "case3",
            "method": "predict",
            "params": {"smiles": ["CCO", "CCN"], "chem_id": ["ethanol", "ethylamine"], "aggregate_scores": True}
        }
    },
    {
        "desc": "molecules结构化用法",
        "data": {
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
    },
    {
        "desc": "单分子，指定chem_id",
        "data": {
            "id": "case5",
            "method": "predict",
            "params": {"smiles": "CCO", "chem_id": "ethanol", "aggregate_scores": False}
        }
    },
]

for case in test_cases:
    print(f"\n==== {case['desc']} ====")
    r = requests.post("http://localhost:8000/mcp", json=case["data"], stream=True)
    for line in r.iter_lines():
        if line:
            print(line.decode())
