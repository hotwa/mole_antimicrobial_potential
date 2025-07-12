FROM nvidia/cuda:11.8.0-cudnn8-devel-ubuntu22.04

# 1. 安装依赖
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git wget bzip2 ca-certificates curl \
        libglib2.0-0 libxext6 libsm6 libxrender1 && \
    rm -rf /var/lib/apt/lists/*

# 2. 安装micromamba
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH=$MAMBA_ROOT_PREFIX/bin:$PATH

RUN wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba

# 3. 设置工作目录
WORKDIR /workspace

# 4. 克隆项目
RUN git clone https://github.com/hotwa/mole_antimicrobial_potential.git

WORKDIR /workspace/mole_antimicrobial_potential

# 5. 复制环境文件到容器（如果你本地有更新的environment.yaml，build时会覆盖仓库里的同名文件）
COPY environment.yaml environment.yaml

# 6. 用micromamba创建环境
RUN micromamba env create -f environment.yaml -n mole -y && \
    micromamba clean -a -y

ENV MODE=api
ENV HOST=0.0.0.0
ENV PORT=8000
ENV RELOAD=--reload
ENV WORKERS=1

# 7. 让micromamba shell自动激活环境（bash、ipython、python等都继承环境）
SHELL ["/usr/local/bin/micromamba", "run", "-n", "mole", "/bin/bash", "-c"]

# 8. 默认进入bash
CMD ["/start.sh"]
