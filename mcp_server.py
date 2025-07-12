# mcp_server.py
import argparse
from fastmcp import FastMCP
from predict_api import AntimicrobialPredictor, MoleculeInput

predictor = AntimicrobialPredictor()
mcp = FastMCP("Antimicrobial MCP Server")

@mcp.tool
async def predict_antimicrobial_activity(input_data: MoleculeInput) -> list:
    """
    分子抗菌潜力预测（MCP标准工具，支持单分子、多分子、聚合/非聚合输出、结构化自定义ID等多种用法）。

    工具说明:
        该工具可用于预测一个或多个小分子的抗菌活性。支持自动匹配分子ID、批量预测、结果聚合统计，适合大模型Agent、Bio研发和自动化高通量筛选。

    参数:
        input_data: MoleculeInput
            molecules: List[dict], 可选。结构化分子输入（推荐复杂批量、带自定义ID场景）。每个dict支持:
                - smiles: str    必填，分子的SMILES字符串
                - chem_id: str   可选，分子自定义ID（如"ethanol"），不填则自动编号
            smiles: str | List[str]，可选。单个或多个SMILES（与chem_id一一对应，结构化用molecules优先）。
            chem_id: str | List[str]，可选。分子ID，推荐与smiles配套用。
            aggregate_scores: bool，是否聚合结果。True时返回汇总分数，False返回每个分子对每个菌株的概率。（默认False）
            app_threshold: float，可选。抑制判定阈值（概率>=此值为有效抑制，默认0.04374）。
            min_nkill: int，可选。定义广谱抗菌的“最小被抑制菌株数”下限（默认10）。

    输入用法示例（适配所有Agent/平台）:
        # 1. 单分子预测（详细菌株分布）
        {
            "smiles": "CCO"
        }
        # 2. 单分子聚合结果
        {
            "smiles": "CCO",
            "aggregate_scores": true
        }
        # 3. 多分子批量
        {
            "smiles": ["CCO", "CCN"]
        }
        # 4. 多分子聚合+自定义ID
        {
            "smiles": ["CCO", "CCN"],
            "chem_id": ["ethanol", "ethylamine"],
            "aggregate_scores": true
        }
        # 5. 结构化输入（推荐复杂批量/带自定义ID）
        {
            "molecules": [
                {"smiles": "CCO", "chem_id": "ethanol"},
                {"smiles": "CCN", "chem_id": "ethylamine"}
            ],
            "aggregate_scores": true
        }
        # 6. 自定义判定阈值和广谱下限
        {
            "smiles": "CCO",
            "aggregate_scores": true,
            "app_threshold": 0.1,
            "min_nkill": 5
        }

    返回:
        - 若 aggregate_scores=False:
            List[dict]: 每个分子的每个菌株的预测概率和抑制判定
            [
                {
                    "pred_id": "ethanol:Escherichia coli (NT5077)",
                    "antimicrobial_predictive_probability": 0.00005,
                    "growth_inhibition": 0
                },
                ...
            ]
        - 若 aggregate_scores=True:
            List[dict]: 每个分子的综合抗菌分数/抑制菌株数/广谱性
            [
                {
                    "chem_id": "ethanol",
                    "apscore_total": -11.7,
                    "apscore_gnegative": -11.6,
                    "apscore_gpositive": -11.8,
                    "ginhib_total": 0,
                    "ginhib_gnegative": 0,
                    "ginhib_gpositive": 0,
                    "broad_spectrum": 0
                },
                ...
            ]

    字段解释:
        - pred_id: 分子ID+菌株（仅非聚合结果），唯一定位某分子某菌株
        - antimicrobial_predictive_probability: 该分子对该菌株的抑菌概率（0~1）
        - growth_inhibition: 是否有效抑制（1是, 0否, 阈值由app_threshold控制）
        - chem_id: 分子的自定义ID（或自动编号，如mol1、ethanol等）
        - apscore_total: 全部菌株的几何均值分数（越低代表更强抗菌潜力）
        - apscore_gnegative / apscore_gpositive: G- / G+ 菌株的分数（对比细分抗菌能力）
        - ginhib_total: 总共被有效抑制的菌株数
        - ginhib_gnegative / ginhib_gpositive: G-/G+菌株被抑制数
        - broad_spectrum: 是否广谱抑制（1为广谱，0为否，判据由min_nkill控制）

    业务建议:
        - 高通量分析可批量传入smiles数组或molecules结构体
        - 推荐结果聚合（aggregate_scores=True）用于药物筛选
        - chem_id自定义便于和数据库对接

    LLM/Agent调用建议:
        - 可自动读取字段类型和示例，支持智能参数补全
        - 支持OpenAPI、Schema、Claude等主流工具发现

    文档:
        - 全功能及字段说明见项目README、/mcp/schema端点、llms.txt协议文档
    """
    return await predictor.predict(input_data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--transport", default="sse")  # 也可传 http、stdio
    args = parser.parse_args()
    mcp.run(host=args.host, port=args.port, transport=args.transport)
