#!/usr/bin/env python3
"""
MCP SSE服务连接测试脚本
测试MolE抗菌预测服务的MCP SSE连接
"""

import requests
import json

def test_mcp_connection():
    """测试MCP SSE服务连接"""
    print("=" * 60)
    print("MCP SSE服务连接测试")
    print("=" * 60)
    
    # 测试1: 检查服务健康状态
    print("\n1. 检查服务健康状态...")
    try:
        response = requests.get("http://localhost:8000/health", timeout=5)
        if response.status_code == 200:
            health_data = response.json()
            print(f"   ✅ 服务状态: {health_data['status']}")
            print(f"   ✅ 模型加载: {health_data['model_loaded']}")
            print(f"   ✅ 设备: {health_data['device']}")
        else:
            print(f"   ❌ 健康检查失败: {response.status_code}")
            return False
    except Exception as e:
        print(f"   ❌ 连接失败: {e}")
        return False
    
    # 测试2: 测试MCP端点
    print("\n2. 测试MCP端点...")
    try:
        # 发送MCP请求
        mcp_request = {
            "id": "test_connection",
            "method": "predict",
            "params": {
                "smiles": "CCO",  # 乙醇
                "aggregate_scores": True
            }
        }
        
        response = requests.post(
            "http://localhost:8000/mcp",
            json=mcp_request,
            headers={"Content-Type": "application/json"},
            timeout=30,
            stream=True
        )
        
        events = []
        for line in response.iter_lines():
            if line:
                line_str = line.decode()
                events.append(line_str)
        
        # 验证事件序列
        has_start = any("event: start" in e for e in events)
        has_data = any("event: data" in e for e in events)
        has_end = any("event: end" in e for e in events)
        
        if has_start and has_data and has_end:
            print("   ✅ SSE流式响应正常")
            print("   ✅ 事件序列完整")
            
            # 解析数据事件
            for event in events:
                if "event: data" in event:
                    # 获取下一行的数据
                    idx = events.index(event)
                    if idx + 1 < len(events):
                        data_line = events[idx + 1]
                        if data_line.startswith("data:"):
                            try:
                                data_json = json.loads(data_line[5:])
                                if "result" in data_json:
                                    result = data_json["result"]
                                    if result and len(result) > 0:
                                        print(f"   ✅ 预测结果: {len(result)} 个菌株预测")
                                        print(f"   ✅ 化合物: {result[0].get('chem_id', 'mol1')}")
                                        print(f"   ✅ 广谱分数: {result[0].get('apscore_total', 'N/A')}")
                            except:
                                pass
            return True
        else:
            print(f"   ❌ SSE事件缺失")
            print(f"   Events: {events}")
            return False
            
    except Exception as e:
        print(f"   ❌ MCP测试失败: {e}")
        return False
    
    # 测试3: 测试健康检查端点
    print("\n3. 测试健康检查端点...")
    endpoints = ["/health", "/health/detailed", "/health/ready", "/health/alive"]
    for endpoint in endpoints:
        try:
            response = requests.get(f"http://localhost:8000{endpoint}", timeout=5)
            if response.status_code == 200:
                print(f"   ✅ {endpoint}: 正常")
            else:
                print(f"   ❌ {endpoint}: {response.status_code}")
        except Exception as e:
            print(f"   ❌ {endpoint}: {e}")
    
    return True

def show_mcp_config():
    """显示MCP配置信息"""
    print("\n" + "=" * 60)
    print("MCP配置信息")
    print("=" * 60)
    print("""
服务地址: http://localhost:8000/mcp
传输协议: SSE (Server-Sent Events)

MCP配置文件位置:
~/.vscode-server/data/User/globalStorage/saoudrizwan.claude-dev/settings/cline_mcp_settings.json

配置内容:
{
  "mcpServers": {
    "mole-antimicrobial-prediction": {
      "command": "npx",
      "args": [
        "mcp-remote",
        "http://localhost:8000/mcp",
        "--transport",
        "sse"
      ],
      "description": "MolE Antimicrobial Prediction Service"
    }
  }
}

可用工具:
- predict_antimicrobial_activity: 预测小分子抗菌活性

测试用例:
1. 健康检查: GET /health
2. MCP预测: POST /mcp {"method": "predict", "params": {"smiles": "CCO"}}
3. 聚合结果: POST /mcp {"method": "predict", "params": {"smiles": "CCO", "aggregate_scores": true}}
""")

if __name__ == "__main__":
    print("开始MCP SSE服务连接测试...\n")
    
    if test_mcp_connection():
        print("\n" + "=" * 60)
        print("🎉 所有测试通过！MCP SSE服务配置成功！")
        print("=" * 60)
        show_mcp_config()
    else:
        print("\n" + "=" * 60)
        print("❌ 测试失败，请检查服务状态")
        print("=" * 60)
