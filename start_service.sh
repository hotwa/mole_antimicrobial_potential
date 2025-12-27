#!/bin/bash
# 启动服务脚本

echo "启动MCP服务..."
echo "模式: http"
echo "主机: 0.0.0.0"
echo "端口: 8000"
echo ""

# 检查是否在正确的目录
if [ ! -f "mcp_server.py" ]; then
    echo "错误: 请在项目根目录运行此脚本"
    exit 1
fi

# 检查模型文件是否存在
if [ ! -f "pretrained_model/model_ginconcat_btwin_100k_d8000_l0.0001/config.yaml" ]; then
    echo "警告: 模型配置文件不存在"
    echo "请确保模型文件已下载到 pretrained_model/model_ginconcat_btwin_100k_d8000_l0.0001/"
    echo "继续启动，但预测功能将不可用..."
    echo ""
fi

# 启动服务
python mcp_server.py --host 0.0.0.0 --port 8000 --transport http
