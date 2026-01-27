#!/usr/bin/env python3
"""
简化版增强MCP服务器测试 - 带超时和详细日志
"""

import pytest

pytest.skip("Legacy in-process MCP tests; use test_mcp_http_e2e.py", allow_module_level=True)

import asyncio
import json
import sys
import os
import time

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
    return True


async def test_model_loading():
    """测试模型加载"""
    print("=" * 60)
    print("测试 2: 模型加载")
    print("=" * 60)
    
    print("正在加载模型，这可能需要一些时间...")
    start_time = time.time()
    
    predictor = get_predictor()
    
    # 等待模型加载完成
    try:
        await asyncio.wait_for(predictor.ensure_loaded(), timeout=120)
        load_time = time.time() - start_time
        print(f"✓ 模型加载完成，耗时: {load_time:.2f}秒")
        
        status = predictor.get_status()
        print(f"  - 加载状态: {status.loaded}")
        print(f"  - 使用设备: {status.device}")
        return True
    except asyncio.TimeoutError:
        print("✗ 模型加载超时")
        return False
    except Exception as e:
        print(f"✗ 模型加载失败: {e}")
        return False


async def test_predict_simple():
    """测试简单预测"""
    print("=" * 60)
    print("测试 3: 简单预测")
    print("=" * 60)
    
    input_data = MoleculeInput(
        smiles="CCO",
        aggregate_scores=True,
        app_threshold=0.04374140128493309,
        min_nkill=10
    )
    
    print("输入参数:")
    print(f"  smiles: {input_data.smiles}")
    print(f"  aggregate_scores: {input_data.aggregate_scores}")
    print()
    
    predictor = get_predictor()
    
    try:
        print("正在执行预测...")
        start_time = time.time()
        
        result = await asyncio.wait_for(
            predictor.predict(input_data),
            timeout=60
        )
        
        predict_time = time.time() - start_time
        print(f"✓ 预测完成，耗时: {predict_time:.2f}秒")
        print(f"  模式: {result['mode']}")
        print(f"  结果数量: {len(result['items'])}")
        
        if result['items']:
            item = result['items'][0]
            print(f"  第一个结果:")
            for key, value in item.items():
                if isinstance(value, float):
                    print(f"    {key}: {value:.6f}")
                else:
                    print(f"    {key}: {value}")
        
        return True
        
    except asyncio.TimeoutError:
        print("✗ 预测超时")
        return False
    except Exception as e:
        print(f"✗ 预测失败: {e}")
        import traceback
        traceback.print_exc()
        return False


async def test_input_validation():
    """测试输入验证"""
    print("=" * 60)
    print("测试 4: 输入验证")
    print("=" * 60)
    
    tests = [
        {
            "name": "长度不匹配",
            "input": {"smiles": ["CCO", "CCN"], "chem_id": ["ethanol"]},
            "should_fail": True
        },
        {
            "name": "无效SMILES",
            "input": {"smiles": "invalid_smiles"},
            "should_fail": True
        },
        {
            "name": "有效输入",
            "input": {"smiles": "CCO", "aggregate_scores": True},
            "should_fail": False
        }
    ]
    
    passed = 0
    for test in tests:
        print(f"  测试: {test['name']}")
        try:
            MoleculeInput(**test['input'])
            if test['should_fail']:
                print("    ✗ 应该失败但通过了")
            else:
                print("    ✓ 正确通过")
                passed += 1
        except Exception as e:
            if test['should_fail']:
                print(f"    ✓ 正确捕获错误: {type(e).__name__}")
                passed += 1
            else:
                print(f"    ✗ 应该通过但失败: {e}")
    
    print(f"\n输入验证: {passed}/{len(tests)} 通过")
    return passed == len(tests)


async def main():
    """运行所有测试"""
    print("\n" + "=" * 60)
    print("增强版MCP服务器简化测试")
    print("=" * 60)
    print()
    
    results = []
    
    # 测试1: Status
    try:
        result = await test_status()
        results.append(("Status", result))
    except Exception as e:
        print(f"Status测试异常: {e}")
        results.append(("Status", False))
    
    # 测试2: 模型加载
    try:
        result = await test_model_loading()
        results.append(("Model Loading", result))
    except Exception as e:
        print(f"模型加载测试异常: {e}")
        results.append(("Model Loading", False))
    
    # 测试3: 预测
    if any(r[1] for r in results if r[0] == "Model Loading"):
        try:
            result = await test_predict_simple()
            results.append(("Prediction", result))
        except Exception as e:
            print(f"预测测试异常: {e}")
            results.append(("Prediction", False))
    else:
        print("跳过预测测试（模型未加载）")
        results.append(("Prediction", False))
    
    # 测试4: 输入验证
    try:
        result = await test_input_validation()
        results.append(("Input Validation", result))
    except Exception as e:
        print(f"输入验证测试异常: {e}")
        results.append(("Input Validation", False))
    
    # 总结
    print("\n" + "=" * 60)
    print("测试总结")
    print("=" * 60)
    
    for name, result in results:
        status = "✓ 通过" if result else "✗ 失败"
        print(f"  {name}: {status}")
    
    total_passed = sum(1 for _, r in results if r)
    total_tests = len(results)
    
    print(f"\n总成绩: {total_passed}/{total_tests} 通过")
    
    return 0 if total_passed == total_tests else 1


if __name__ == "__main__":
    exit_code = asyncio.run(main())
    sys.exit(exit_code)
