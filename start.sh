#!/bin/bash
set -e

MODE="${MODE:-api}"
HOST="${HOST:-0.0.0.0}"

PORT_API="${PORT_API:-8000}"
WORKERS_API="${WORKERS_API:-2}"

PORT_MCP="${PORT_MCP:-8001}"
MCP_TRANSPORT="${MCP_TRANSPORT:-http}"

# GPU 配置：优先使用 GPU_IDS（避免与 nvidia runtime 冲突）
GPU_ID="${GPU_IDS:-0}"
export CUDA_VISIBLE_DEVICES="$GPU_ID"

echo "[INFO] MODE: $MODE"
echo "[INFO] HOST: $HOST"
echo "[INFO] CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"

if [ "$MODE" = "api" ]; then
    echo "[INFO] Starting FastAPI service"
    echo "[INFO] PORT_API: $PORT_API"
    echo "[INFO] WORKERS_API: $WORKERS_API"
    exec uvicorn api_server:app --host "$HOST" --port "$PORT_API" --workers "$WORKERS_API"

elif [ "$MODE" = "mcp_enhanced" ]; then
    echo "[INFO] Starting Enhanced FastMCP service"
    echo "[INFO] PORT_MCP: $PORT_MCP"
    echo "[INFO] MCP_TRANSPORT: $MCP_TRANSPORT"
    exec python mcp_server_enhanced.py
else
    echo "[ERROR] Unknown MODE: $MODE"
    echo "[ERROR] Supported MODE values: api, mcp_enhanced"
    exit 1
fi
