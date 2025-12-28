"""
Shared predictor service for both FastAPI and FastMCP
避免重复加载模型，提供统一的预测服务
"""

_predictor = None

def get_predictor():
    """
    获取单例 predictor 实例
    避免 FastAPI 和 FastMCP 重复加载模型
    """
    global _predictor
    if _predictor is None:
        # 延迟导入，避免循环依赖和启动时实例化
        from predict_api import AntimicrobialPredictor
        _predictor = AntimicrobialPredictor()
    return _predictor
