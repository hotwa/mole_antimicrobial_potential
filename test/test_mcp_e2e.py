#!/usr/bin/env python3
"""
MCP端到端测试 - 验证协议发现和调用
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

async def test_mcp_discovery():
    """测试MCP协议发现功能"""
    print("=" * 60)
    print("测试 1: MCP 协议发现")
    print("=" * 60)
    
    # 模拟MCP客户端的tools/list调用
    from src.service import get_predictor
    from mcp_server_enhanced import status, predict_antimicrobial_potential, get_strain_info, get_schema_info, get_about
    
    print("✓ 工具列表:")
    print("  - status")
    print("  - predict_antimicrobial_potential")
    print()
    
    print("✓ 资源列表:")
    print("  - resource://strains")
    print("  - resource://schema")
    print("  - resource://about")
    print()
    
    return True


async def test_status_tool():
    """测试status工具"""
    print("=" * 60)
    print("测试 2: status() 工具")
    print("=" * 60)
    
    from mcp_server_enhanced import status
    
    result = await status()
    print("✓ status() 返回:")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    print()
    
    # 验证返回结构
    required_keys = ["service", "version", "loaded", "device", "model_path"]
    for key in required_keys:
        if key not in result:
            print(f"✗ 缺少必需字段: {key}")
            return False
    
    print("✓ 所有必需字段存在")
    return True


async def test_resources():
    """测试资源访问"""
    print("=" * 60)
    print("测试 3: Resources 访问")
    print("=" * 60)
    
    from mcp_server_enhanced import get_strain_info, get_schema_info, get_about
    
    # 测试 strains 资源
    print("测试 3.1: resource://strains")
    strains = await get_strain_info()
    print(f"✓ 返回类型: {type(strains)}")
    print(f"✓ 内容预览: {strains[:100]}...")
    print()
    
    # 测试 schema 资源
    print("测试 3.2: resource://schema")
    schema = await get_schema_info()
    print(f"✓ 返回类型: {type(schema)}")
    print(f"✓ 工具数量: {len(schema.get('tools', {}))}")
    print()
    
    # 测试 about 资源
    print("测试 3.3: resource://about")
    about = await get_about()
    print(f"✓ 返回类型: {type(about)}")
    print(f"✓ 服务: {about.get('service')}")
    print()
    
    return True


async def test_predict_tool():
    """测试predict工具"""
    print("=" * 60)
    print("测试 4: predict_antimicrobial_potential() 工具")
    print("=" * 60)
    
    from mcp_server_enhanced import predict_antimicrobial_potential
    from src.service import get_predictor
    
    # 先确保模型加载
    print("确保模型加载...")
    predictor = get_predictor()
    await predictor.ensure_loaded()
    print("✓ 模型已加载")
    print()
    
    # 测试4.1: 单分子聚合预测
    print("测试 4.1: 单分子聚合预测")
    try:
        result = await predict_antimicrobial_potential(
            smiles="CCO",
            aggregate_scores=True,
            app_threshold=0.04374140128493309,
            min_nkill=10
        )
        print(f"✓ 返回模式: {result['mode']}")
        print(f"✓ 结果数量: {len(result['items'])}")
        if result['items']:
            item = result['items'][0]
            print(f"✓ 示例字段: chem_id={item.get('chem_id')}, apscore_total={item.get('apscore_total'):.4f}")
        print()
    except Exception as e:
        print(f"✗ 失败: {e}")
        return False
    
    # 测试4.2: 单分子详细预测
    print("测试 4.2: 单分子详细预测")
    try:
        result = await predict_antimicrobial_potential(
            smiles="CCO",
            aggregate_scores=False
        )
        print(f"✓ 返回模式: {result['mode']}")
        print(f"✓ 结果数量: {len(result['items'])}")
        if result['items']:
            item = result['items'][0]
            print(f"✓ 示例字段: pred_id={item.get('pred_id')}, probability={item.get('antimicrobial_predictive_probability'):.4f}")
        print()
    except Exception as e:
        print(f"✗ 失败: {e}")
        return False
    
    # 测试4.3: 批量预测
    print("测试 4.3: 批量预测")
    try:
        result = await predict_antimicrobial_potential(
            molecules=[
                {"smiles": "CCO", "chem_id": "ethanol"},
                {"smiles": "CCN", "chem_id": "ethylamine"}
            ],
            aggregate_scores=True,
            min_nkill=5
        )
        print(f"✓ 返回模式: {result['mode']}")
        print(f"✓ 结果数量: {len(result['items'])}")
        for item in result['items']:
            print(f"  - {item.get('chem_id')}: ginhib={item.get('ginhib_total')}, broad={item.get('broad_spectrum')}")
        print()
    except Exception as e:
        print(f"✗ 失败: {e}")
        return False
    
    return True


async def test_error_handling():
    """测试错误处理"""
    print("=" * 60)
    print("测试 5: 错误处理")
    print("=" * 60)
    
    from mcp_server_enhanced import predict_antimicrobial_potential
    
    # 测试5.1: 缺少必需参数
    print("测试 5.1: 缺少必需参数")
    try:
        await predict_antimicrobial_potential()
        print("✗ 应该抛出异常")
        return False
    except ValueError as e:
        print(f"✓ 正确捕获异常: {e}")
    print()
    
    # 测试5.2: 长度不匹配
    print("测试 5.2: chem_id 长度不匹配")
    try:
        await predict_antimicrobial_potential(
            smiles=["CCO", "CCN"],
            chem_id=["ethanol"],
            aggregate_scores=True
        )
        print("✗ 应该抛出异常")
        return False
    except ValueError as e:
        print(f"✓ 正确捕获异常: {e}")
    print()
    
    # 测试5.3: 无效SMILES
    print("测试 5.3: 无效SMILES")
    try:
        await predict_antimicrobial_potential(
            smiles="invalid_smiles",
            aggregate_scores=True
        )
        print("✗ 应该抛出异常")
        return False
    except Exception as e:
        print(f"✓ 正确捕获异常: {type(e).__name__}")
    print()
    
    return True


async def test_unified_response_structure():
    """测试统一响应结构"""
    print("=" * 60)
    print("测试 6: 统一响应结构验证")
    print("=" * 60)
    
    from mcp_server_enhanced import predict_antimicrobial_potential
    from src.service import get_predictor
    
    # 确保模型加载
    predictor = get_predictor()
    await predictor.ensure_loaded()
    
    # 测试聚合模式
    print("测试 6.1: 聚合模式响应结构")
    agg_result = await predict_antimicrobial_potential(
        smiles="CCO",
        aggregate_scores=True
    )
    
    if agg_result.get("mode") != "aggregate":
        print(f"✗ mode 应为 'aggregate', 得到: {agg_result.get('mode')}")
        return False
    
    if "items" not in agg_result:
        print("✗ 缺少 'items' 字段")
        return False
    
    if agg_result["items"]:
        item = agg_result["items"][0]
        required_fields = ["chem_id", "apscore_total", "ginhib_total", "broad_spectrum"]
        for field in required_fields:
            if field not in item:
                print(f"✗ 聚合结果缺少字段: {field}")
                return False
    
    print("✓ 聚合模式结构正确")
    print()
    
    # 测试详细模式
    print("测试 6.2: 详细模式响应结构")
    detail_result = await predict_antimicrobial_potential(
        smiles="CCO",
        aggregate_scores=False
    )
    
    if detail_result.get("mode") != "per_strain":
        print(f"✗ mode 应为 'per_strain', 得到: {detail_result.get('mode')}")
        return False
    
    if "items" not in detail_result:
        print("✗ 缺少 'items' 字段")
        return False
    
    if detail_result["items"]:
        item = detail_result["items"][0]
        required_fields = ["pred_id", "antimicrobial_predictive_probability", "growth_inhibition"]
        for field in required_fields:
            if field not in item:
                print(f"✗ 详细结果缺少字段: {field}")
                return False
    
    print("✓ 详细模式结构正确")
    print()
    
    return True


async def main():
    """运行所有E2E测试"""
    print("\n" + "=" * 60)
    print("MCP 端到端测试套件")
    print("=" * 60)
    print()
    
    tests = [
        ("MCP 协议发现", test_mcp_discovery),
        ("Status 工具", test_status_tool),
        ("Resources 访问", test_resources),
        ("Predict 工具", test_predict_tool),
        ("错误处理", test_error_handling),
        ("统一响应结构", test_unified_response_structure),
    ]
    
    results = []
    
    for name, test_func in tests:
        try:
            result = await test_func()
            results.append((name, result))
        except Exception as e:
            print(f"✗ {name} 测试异常: {e}")
            import traceback
            traceback.print_exc()
            results.append((name, False))
    
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
