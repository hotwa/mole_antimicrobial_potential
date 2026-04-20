# MolE checkpoint 来源说明

这个目录里的 `model.pth` 是本仓库使用的 MolE 预训练权重。

公开来源来自上游 MolE 项目：

- https://github.com/rolayoalarcon/MolE

上游仓库的说明写明了该预训练模型的放置位置：

- 将 `model.pth` 放到 `ckpt/gin_concat_R1000_E8000_lambda0.0001/checkpoints/`

本仓库里当前保留了一个便捷副本：

- `ckpt/model.pth`

以及兼容路径：

- `pretrained_model/model_ginconcat_btwin_100k_d8000_l0.0001/model.pth`

如果本地缺失这个权重，优先从上游 MolE 仓库重新下载，再按上面的路径放置即可。
