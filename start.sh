#!/bin/bash  
set -e  
  
# 支持常用参数，均可通过环境变量传递，未传递时使用默认值  
MODE="${MODE:-api}"         # 启动模式: api (FastAPI) 或 mcp (FastMCP标准MCP)  
HOST="${HOST:-0.0.0.0}"    # 监听地址  

# FastAPI 配置
PORT_API="${PORT_API:-8000}"       # API端口  
WORKERS_API="${WORKERS_API:-2}"    # API并发进程数  

# FastMCP 配置  
PORT_MCP="${PORT_MCP:-8001}"       # MCP端口  
MCP_TRANSPORT="${MCP_TRANSPORT:-http}"  # MCP传输协议: http, sse, stdio  

# 打印当前配置  
echo "[INFO] MODE: $MODE"  
echo "[INFO] HOST: $HOST"  

# 启动逻辑  
if [ "$MODE" = "api" ]; then  
    echo "[INFO] 启动 FastAPI 服务"  
    echo "[INFO] PORT_API: $PORT_API"  
    echo "[INFO] WORKERS_API: $WORKERS_API"  
    exec uvicorn predict_api:app --host "$HOST" --port "$PORT_API" --workers "$WORKERS_API"  
  
elif [ "$MODE" = "mcp" ]; then  
    echo "[INFO] 启动 FastMCP 标准 MCP 服务"  
    echo "[INFO] PORT_MCP: $PORT_MCP"  
    echo "[INFO] MCP_TRANSPORT: $MCP_TRANSPORT"  
    # 强制单进程：不要用 uvicorn workers  
    exec python mcp_server.py  
  
else  
    echo "[ERROR] 未知启动模式：$MODE"  
    echo "[ERROR] 支持的模式: api, mcp"  
    exit 1  
fi

# # API 模式使用 uvicorn 额外参数  
# MODE=api EXTRA_UVICORN_ARGS="--access-log --log-config logging.json" ./start.sh  
  
# # MCP 模式使用 FastMCP 额外参数（如日志级别）  
# MODE=mcp EXTRA_FASTMCP_ARGS="--log-level DEBUG" ./start.sh  
  
# # MCP 模式使用 STDIO 传输  
# MODE=mcp TRANSPORT=stdio ./start.sh
