#!/bin/bash  
set -e  
  
# 支持常用参数，均可通过环境变量传递，未传递时使用默认值  
MODE="${MODE:-mcp}"         # 启动模式: api（老接口） 或 mcp（MCP标准）  
HOST="${HOST:-0.0.0.0}"    # 监听地址  
PORT="${PORT:-8000}"       # 端口  
RELOAD="${RELOAD:-}"       # 是否开启 --reload，默认不开，可设为 "--reload"  
WORKERS="${WORKERS:-1}"    # Uvicorn并发进程数，默认1  
TRANSPORT="${TRANSPORT:-http}"  # MCP传输协议: stdio, http, sse  
  
# 额外参数（按模式分离）  
EXTRA_UVICORN_ARGS="${EXTRA_UVICORN_ARGS:-}"      # API模式专用的uvicorn参数  
EXTRA_FASTMCP_ARGS="${EXTRA_FASTMCP_ARGS:-}"      # MCP模式专用的fastmcp参数  
  
# 打印当前配置  
echo "[INFO] MODE: $MODE"  
echo "[INFO] HOST: $HOST"  
echo "[INFO] PORT: $PORT"  
if [ "$MODE" = "api" ]; then  
    echo "[INFO] RELOAD: $RELOAD"  
    echo "[INFO] WORKERS: $WORKERS"  
    echo "[INFO] EXTRA_UVICORN_ARGS: $EXTRA_UVICORN_ARGS"  
elif [ "$MODE" = "mcp" ]; then  
    echo "[INFO] TRANSPORT: $TRANSPORT"  
    echo "[INFO] EXTRA_FASTMCP_ARGS: $EXTRA_FASTMCP_ARGS"  
fi  
  
# 启动逻辑  
if [ "$MODE" = "api" ]; then  
    echo "[INFO] 启动老API predict_api.py"  
    exec uvicorn predict_api:app --host "$HOST" --port "$PORT" $RELOAD --workers "$WORKERS" $EXTRA_UVICORN_ARGS  
elif [ "$MODE" = "mcp" ]; then  
    echo "[INFO] 启动MCP标准接口 mcp_server.py"  
    exec fastmcp run mcp_server.py --transport "$TRANSPORT" --host "$HOST" --port "$PORT" $EXTRA_FASTMCP_ARGS  
else  
    echo "[ERROR] 未知启动模式：$MODE"  
    exit 1  
fi

# # API 模式使用 uvicorn 额外参数  
# MODE=api EXTRA_UVICORN_ARGS="--access-log --log-config logging.json" ./start.sh  
  
# # MCP 模式使用 FastMCP 额外参数（如日志级别）  
# MODE=mcp EXTRA_FASTMCP_ARGS="--log-level DEBUG" ./start.sh  
  
# # MCP 模式使用 STDIO 传输  
# MODE=mcp TRANSPORT=stdio ./start.sh