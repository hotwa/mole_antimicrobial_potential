#!/usr/bin/env python3
"""
标准 FastMCP 服务器
提供与 predict_api.py 相同的预测功能，但使用标准 MCP 协议
"""

import os
import asyncio
from fastmcp import FastMCP

# 创建 FastMCP 实例
mcp = FastMCP("mole-antimicrobial-prediction")

@mcp.tool
async def predict(
    molecules=None,
    smiles=None,
    chem_id=None,
    aggregate_scores: bool = False,
    app_threshold: float = 0.04374140128493309,
    min_nkill: int = 10,
):
    """
    预测小分子的抗菌潜力
    
    参数:
        molecules: 结构化输入，每个dict含smiles/chem_id
        smiles: SMILES字符串或字符串列表
        chem_id: 化合物ID或ID列表
        aggregate_scores: 是否聚合输出
        app_threshold: 抑制阈值
        min_nkill: 广谱判断的最小菌株数
    
    返回:
        预测结果列表
    """
    # 延迟导入，避免启动时实例化
    from src.models import MoleculeInput
    from src.service import get_predictor

    predictor = get_predictor()
    input_data = MoleculeInput(
        molecules=molecules,
        smiles=smiles,
        chem_id=chem_id,
        aggregate_scores=aggregate_scores,
        app_threshold=app_threshold,
        min_nkill=min_nkill,
    ).normalize()
    await predictor.ensure_loaded()
    return await predictor.predict(input_data)

def main():
    """启动 FastMCP 服务器"""
    host = os.getenv("HOST", "0.0.0.0")
    port = int(os.getenv("PORT_MCP", "8001"))
    transport = os.getenv("MCP_TRANSPORT", "http")  # http recommended
    
    print(f"[INFO] Starting FastMCP server on {host}:{port} via {transport}")
    mcp.run(transport=transport, host=host, port=port)

if __name__ == "__main__":
    main()
