#!/usr/bin/env bash
set -euo pipefail

echo "[TEST] MCP connect via mcp-remote"
echo "测试标准 MCP 服务器连接..."
echo ""

# 测试 MCP 服务器是否响应
if curl -s http://localhost:8001/mcp > /dev/null; then
    echo "✅ MCP 服务器正在运行"
else
    echo "❌ MCP 服务器未响应"
    exit 1
fi

# 使用 mcp-remote 连接
echo "正在启动 mcp-remote..."
npx mcp-remote http://localhost:8001/mcp --transport http --debug <<<"" 2>&1 | head -20

echo ""
echo "[TEST] 完成"
