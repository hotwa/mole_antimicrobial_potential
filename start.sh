#!/bin/bash
# 设置默认值
HOST=${UVICORN_HOST:-0.0.0.0}
PORT=${UVICORN_PORT:-8000}
RELOAD=${UVICORN_RELOAD:-"--reload"}

exec uvicorn predict_api:app --host $HOST --port $PORT $RELOAD
