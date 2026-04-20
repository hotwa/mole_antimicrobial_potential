#!/usr/bin/env python3
"""
测试增强版MCP服务器的功能
"""

import pytest

pytest.skip("Legacy in-process MCP tests; use test_mcp_http_e2e.py", allow_module_level=True)

import asyncio
import json
import sys
import os

# 添加项目路径
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.models import MoleculeInput, ServiceStatus
from src.service import get_predictor


async def test_status():
    """测试status工具"""
    print("=" * 60)
    print("测试 1: status() 工具")
    print("=" * 60)
    
    predictor = get_predictor()
    status_info = predictor.get_status()
    
    result = {
        "service": "mole-antimicrobial-prediction",
        "version": "1.0.0",
        "loaded": status_info.loaded,
        "device": status_info.device,
        "model_path": status_info.model_path
    }
    
    print("✓ Status 工具返回:")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    print()


async def test_predict_single_aggregated():
    """测试单分子聚合预测"""
    print("=" * 60)
    print("测试 2: 单分子聚合预测 (aggregate_scores=True)")
    print("=" * 60)
    
    # 模拟输入
    input_data = MoleculeInput(
        smiles="CCO",  # 乙醇
        aggregate_scores=True,
        app_threshold=0.04374140128493309,
        min_nkill=10
    )
    
    print("输入参数:")
    print(f"  smiles: {input_data.smiles}")
    print(f"  aggregate_scores: {input_data.aggregate_scores}")
    print(f"  app_threshold: {input_data.app_threshold}")
    print(f"  min_nkill: {input_data.min_nkill}")
    print()
    
    # 执行预测
    predictor = get_predictor()
    try:
        result = await predictor.predict(input_data)
        print("✓ 预测结果:")
        print(f"  模式: {result['mode']}")
        print(f"  结果数量: {len(result['items'])}")
        if result['items']:
            item = result['items'][0]
            print(f"  化合物ID: {item.get('chem_id', 'N/A')}")
            print(f"  总分数: {item.get('apscore_total', 'N/A')}")
            print(f"  抑制菌株数: {item.get('ginhib_total', 'N/A')}")
            print(f"  广谱抗菌: {item.get('broad_spectrum', 'N/A')}")
    except Exception as e:
        print(f"✗ 预测失败: {e}")
    print()


async def test_predict_single_detailed():
    """测试单分子详细预测"""
    print("=" * 60)
    print("测试 3: 单分子详细预测 (aggregate_scores=False)")
    print("=" * 60)
    
    input_data = MoleculeInput(
        smiles="CCO",
        aggregate_scores=False,
        app_threshold=0.04374140128493309
    )
    
    print("输入参数:")
    print(f"  smiles: {input_data.smiles}")
    print(f"  aggregate_scores: {input_data.aggregate_scores}")
    print()
    
    predictor = get_predictor()
    try:
        result = await predictor.predict(input_data)
        print("✓ 预测结果:")
        print(f"  模式: {result['mode']}")
        print(f"  结果数量: {len(result['items'])}")
        if result['items']:
            # 显示前3个菌株的结果
            for i, item in enumerate(result['items'][:3]):
                print(f"  [{i+1}] {item.get('pred_id', 'N/A')}:")
                print(f"      概率: {item.get('antimicrobial_predictive_probability', 'N/A'):.4f}")
                print(f"      抑制: {item.get('growth_inhibition', 'N/A')}")
    except Exception as e:
        print(f"✗ 预测失败: {e}")
    print()


async def test_predict_batch():
    """测试批量预测"""
    print("=" * 60)
    print("测试 4: 批量预测 (molecules)")
    print("=" * 60)
    
    input_data = MoleculeInput(
        molecules=[
            {"smiles": "CCO", "chem_id": "ethanol"},
            {"smiles": "CCN", "chem_id": "ethylamine"}
        ],
        aggregate_scores=True,
        min_nkill=5  # 降低阈值以便测试
    )
    
    print("输入参数:")
    print(f"  分子数量: {len(input_data.molecules)}")
    for mol in input_data.molecules:
        print(f"    - {mol['chem_id']}: {mol['smiles']}")
    print(f"  aggregate_scores: {input_data.aggregate_scores}")
    print(f"  min_nkill: {input_data.min_nkill}")
    print()
    
    predictor = get_predictor()
    try:
        result = await predictor.predict(input_data)
        print("✓ 预测结果:")
        print(f"  模式: {result['mode']}")
        print(f"  结果数量: {len(result['items'])}")
        for item in result['items']:
            print(f"  {item.get('chem_id', 'N/A')}: "
                  f"总分数={item.get('apscore_total', 'N/A'):.4f}, "
                  f"抑制数={item.get('ginhib_total', 'N/A')}, "
                  f"广谱={item.get('broad_spectrum', 'N/A')}")
    except Exception as e:
        print(f"✗ 预测失败: {e}")
    print()


async def test_input_validation():
    """测试输入验证"""
    print("=" * 60)
    print("测试 5: 输入验证")
    print("=" * 60)
    
    # 测试1: 长度不匹配
    print("测试 5.1: chem_id 长度不匹配")
    try:
        MoleculeInput(
            smiles=["CCO", "CCN"],
            chem_id=["ethanol"],  # 只有1个，但smiles有2个
            aggregate_scores=True
        )
        print("✗ 应该抛出验证错误")
    except ValueError as e:
        print(f"✓ 正确捕获错误: {e}")
    
    # 测试2: 无效SMILES
    print("\n测试 5.2: 无效SMILES")
    try:
        MoleculeInput(
            smiles="invalid_smiles",
            aggregate_scores=True
        )
        print("✗ 应该抛出验证错误")
    except ValueError as e:
        print(f"✓ 正确捕获错误: {e}")
    
    # 测试3: 阈值超出范围
    print("\n测试 5.3: 阈值超出范围")
    try:
        MoleculeInput(
            smiles="CCO",
            app_threshold=1.5,
            aggregate_scores=True
        )
        print("✗ 应该抛出验证错误")
    except ValueError as e:
        print(f"✓ 正确捕获错误: {e}")
    
    print()


async def test_resources():
    """测试资源访问"""
    print("=" * 60)
    print("测试 6: Resources 访问")
    print("=" * 60)
    
    # 模拟资源访问
    print("资源列表:")
    print("  - resource://strains: 菌株信息")
    print("  - resource://schema: API Schema")
    print("  - resource://about: 服务信息")
    print()
    
    # 模拟about资源
    about = {
        "service": "mole-antimicrobial-prediction",
        "version": "1.0.0",
        "description": "AI-powered antimicrobial activity prediction using MolE + XGBoost",
        "paper": {
            "title": "Pre-trained molecular representations enable antimicrobial discovery",
            "journal": "Nature Communications",
            "year": 2025
        }
    }
    print("✓ resource://about 示例:")
    print(json.dumps(about, indent=2, ensure_ascii=False))
    print()


async def main():
    """运行所有测试"""
    print("\n" + "=" * 60)
    print("增强版MCP服务器功能测试")
    print("=" * 60)
    print()
    
    try:
        await test_status()
        await test_predict_single_aggregated()
        await test_predict_single_detailed()
        await test_predict_batch()
        await test_input_validation()
        await test_resources()
        
        print("=" * 60)
        print("✓ 所有测试完成")
        print("=" * 60)
        
    except Exception as e:
        print(f"\n✗ 测试过程中发生错误: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    exit_code = asyncio.run(main())
    sys.exit(exit_code)
