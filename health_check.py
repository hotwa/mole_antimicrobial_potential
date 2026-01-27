"""
健康检查端点 - 快速诊断服务状态
"""
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import torch

class HealthStatus(BaseModel):
    status: str
    model_loaded: bool
    device: str
    gpu_available: bool
    error: str = None

def add_health_check(app: FastAPI, predictor):
    """
    为FastAPI应用添加健康检查端点
    
    Args:
        app: FastAPI应用实例
        predictor: AntimicrobialPredictor实例
    """
    
    @app.get("/health", response_model=HealthStatus, tags=["Monitoring"])
    async def health_check():
        """
        健康检查端点 - 检查模型加载状态和系统资源
        
        Returns:
            {
                "status": "ready" | "loading" | "error",
                "model_loaded": true | false,
                "device": "cuda:0" | "cpu",
                "gpu_available": true | false,
                "error": "错误信息（如有）"
            }
        """
        try:
            # 检查GPU可用性
            gpu_available = torch.cuda.is_available()
            
            # 检查预测器状态
            if hasattr(predictor, '_model_loaded'):
                model_loaded = predictor._model_loaded
            else:
                # 兼容旧版本，检查模型是否存在
                model_loaded = predictor.model is not None
            
            # 确定状态
            if model_loaded:
                status = "ready"
            elif hasattr(predictor, '_load_task') and predictor._load_task:
                status = "loading"
            else:
                status = "error"
            
            return HealthStatus(
                status=status,
                model_loaded=model_loaded,
                device=predictor.device,
                gpu_available=gpu_available
            )
            
        except Exception as e:
            return HealthStatus(
                status="error",
                model_loaded=False,
                device="unknown",
                gpu_available=False,
                error=str(e)
            )
    
    @app.get("/health/detailed", tags=["Monitoring"])
    async def health_check_detailed():
        """
        详细健康检查 - 包含更多系统信息
        """
        try:
            gpu_available = torch.cuda.is_available()
            gpu_count = torch.cuda.device_count() if gpu_available else 0
            
            # 检查内存
            if gpu_available:
                gpu_memory = []
                for i in range(gpu_count):
                    memory = torch.cuda.memory_allocated(i) / 1024**3  # GB
                    gpu_memory.append(f"GPU {i}: {memory:.2f}GB")
            else:
                gpu_memory = []
            
            # 模型状态
            model_loaded = False
            if hasattr(predictor, '_model_loaded'):
                model_loaded = predictor._model_loaded
            
            return {
                "status": "ready" if model_loaded else "loading",
                "model_loaded": model_loaded,
                "device": predictor.device,
                "gpu": {
                    "available": gpu_available,
                    "count": gpu_count,
                    "memory_used": gpu_memory
                },
                "system": {
                    "python_version": "3.9+",
                    "pytorch_version": torch.__version__
                }
            }
        except Exception as e:
            return {"status": "error", "error": str(e)}
    
    @app.get("/health/ready", tags=["Monitoring"])
    async def readiness_probe():
        """
        就绪探针 - 用于K8s readinessProbe
        """
        if hasattr(predictor, '_model_loaded') and predictor._model_loaded:
            return {"status": "ready"}
        else:
            raise HTTPException(status_code=503, detail="Service not ready")
    
    @app.get("/health/alive", tags=["Monitoring"])
    async def liveness_probe():
        """
        存活探针 - 用于K8s livenessProbe
        """
        return {"status": "alive"}
