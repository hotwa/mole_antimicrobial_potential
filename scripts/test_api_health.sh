#!/usr/bin/env bash
set -euo pipefail

echo "[TEST] FastAPI Health Check"
echo "测试 FastAPI 服务健康状态..."
echo ""

# 测试 FastAPI 健康检查
echo "测试 /health 端点..."
curl -fsS http://localhost:8000/health | jq . 2>/dev/null || curl -fsS http://localhost:8000/health

echo ""
echo "[TEST] 完成"
