<!-- 此文件由机器翻译自 machine_learning.md -->

# 机器学习

## 概述

PathML 为计算病理学提供全面的机器学习功能，包括用于细胞核检测和分割的预构建模型、PyTorch 集成训练工作流程、公共数据集访问和基于 ONNX 的推理部署。该框架将图像预处理与深度学习无缝连接起来，以实现端到端病理学 ML 管道。

## 预建模型

PathML 包括用于细胞核分析的最先进的预训练模型：

### HoVer-Net

**HoVer-Net**（水平和垂直网络）同时执行核实例分割和分类。

**架构：**
- 具有三个预测分支的编码器-解码器结构：
  - **核像素 (NP)** - 核区域的二进制分割
  - **水平-垂直 (HV)** - 到核质心的距离映射
  - **分类 (NC)** - 细胞核类型分类

**细胞核类型：**
1. 上皮细胞
2. 炎症
3. 结缔组织/软组织
4.死亡/坏死
5. 背景

**用途：**
```python
from pathml.ml import HoVerNet
import torch

# Load pre-trained model
model = HoVerNet(
    num_types=5,  # Number of nucleus types
    mode='fast',  # 'fast' or 'original'
    pretrained=True  # Load pre-trained weights
)

# Move to GPU if available
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = model.to(device)

# Inference on tile
tile_image = torch.from_numpy(tile.image).permute(2, 0, 1).unsqueeze(0).float()
tile_image = tile_image.to(device)

with torch.no_grad():
    output = model(tile_image)

# Output contains:
# - output['np']: Nuclear pixel predictions
# - output['hv']: Horizontal-vertical maps
# - output['nc']: Classification predictions
```

**后处理：**
<<<代码块_1>>>

### HACT 网络

**HACTNet**（分层细胞类型网络）通过不确定性量化执行分层细胞核分类。

**特点：**
- 层次分类（粗粒度到细粒度类型）
- 预测的不确定性估计
- 改进了不平衡数据集的性能

<<<代码块_2>>>

## 培训工作流程

### 数据集准备

PathML 提供与 PyTorch 兼容的数据集类：

**平铺数据集：**
<<<代码块_3>>>

**数据模块集成：**
<<<代码块_4>>>

### 训练 HoVer-Net

使用自定义数据训练 HoVer-Net 的完整工作流程：

<<<代码块_5>>>

### PyTorch 闪电集成

PathML 模型与 PyTorch Lightning 集成以简化训练：

<<<代码块_6>>>

## 公共数据集

PathML 提供对公共病理数据集的便捷访问：

### PanNuke 数据集

**PanNuke** 包含来自 19 种组织类型的 7,901 个组织学图像块，以及 5 种细胞类型的细胞核注释。

```python
from pathml.ml.datasets import PanNukeDataModule

# Load PanNuke dataset
pannuke = PanNukeDataModule(
    data_dir='path/to/pannuke',
    batch_size=16,
    num_workers=4,
    tissue_types=None,  # Use all tissue types, or specify list
    fold='all'  # 'fold1', 'fold2', 'fold3', or 'all'
)

# Access dataloaders
train_loader = pannuke.train_dataloader()
val_loader = pannuke.val_dataloader()
test_loader = pannuke.test_dataloader()

# Batch structure
for batch in train_loader:
    images = batch['image']  # Shape: (B, 3, 256, 256)
    inst_map = batch['inst_map']  # Instance segmentation map
    type_map = batch['type_map']  # Cell type map
    np_map = batch['np_map']  # Nuclear pixel map
    hv_map = batch['hv_map']  # Horizontal-vertical distance maps
    tissue_type = batch['tissue_type']  # Tissue category
```

**可用的组织类型：**
乳房、结肠、前列腺、肺、肾、胃、膀胱、食道、子宫颈、肝脏、甲状腺、头颈、睾丸、肾上腺、胰腺、胆管、卵巢、皮肤、子宫

### TCGA 数据集

访问癌症基因组图谱数据集：

```python
from pathml.ml.datasets import TCGADataModule

# Load TCGA dataset
tcga = TCGADataModule(
    data_dir='path/to/tcga',
    cancer_type='BRCA',  # Breast cancer
    batch_size=32,
    tile_size=224
)
```

### 自定义数据集集成

为 PathML 工作流程创建自定义数据集：

```python
from torch.utils.data import Dataset
import numpy as np
from pathlib import Path

class CustomPathologyDataset(Dataset):
    def __init__(self, data_dir, transform=None):
        self.data_dir = Path(data_dir)
        self.image_paths = list(self.data_dir.glob('images/*.png'))
        self.transform = transform

    def __len__(self):
        return len(self.image_paths)

    def __getitem__(self, idx):
        # Load image
        image_path = self.image_paths[idx]
        image = np.array(Image.open(image_path))

        # Load corresponding annotation
        annot_path = self.data_dir / 'annotations' / f'{image_path.stem}.npy'
        annotation = np.load(annot_path)

        # Apply transforms
        if self.transform:
            image = self.transform(image)

        return {
            'image': torch.from_numpy(image).permute(2, 0, 1).float(),
            'annotation': torch.from_numpy(annotation).long(),
            'path': str(image_path)
        }

# Use in PathML workflow
dataset = CustomPathologyDataset('path/to/data')
dataloader = DataLoader(dataset, batch_size=16, shuffle=True, num_workers=4)
```

## 数据增强

应用增强来提高模型泛化能力：

```python
import albumentations as A
from albumentations.pytorch import ToTensorV2

# Define augmentation pipeline
train_transform = A.Compose([
    A.RandomRotate90(p=0.5),
    A.Flip(p=0.5),
    A.ColorJitter(brightness=0.2, contrast=0.2, saturation=0.2, hue=0.1, p=0.5),
    A.GaussianBlur(blur_limit=(3, 7), p=0.3),
    A.ElasticTransform(alpha=1, sigma=50, alpha_affine=50, p=0.3),
    A.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]),
    ToTensorV2()
])

val_transform = A.Compose([
    A.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]),
    ToTensorV2()
])

# Apply to dataset
train_dataset = TileDataset(slide_dataset, transform=train_transform)
val_dataset = TileDataset(val_slide_dataset, transform=val_transform)
```

## 模型评估

### 指标

使用特定于病理学的指标评估模型性能：

```python
from pathml.ml.metrics import (
    dice_coefficient,
    aggregated_jaccard_index,
    panoptic_quality
)

# Dice coefficient for segmentation
dice = dice_coefficient(pred_mask, true_mask)

# Aggregated Jaccard Index (AJI) for instance segmentation
aji = aggregated_jaccard_index(pred_inst, true_inst)

# Panoptic Quality (PQ) for joint segmentation and classification
pq, sq, rq = panoptic_quality(pred_inst, true_inst, pred_types, true_types)

print(f"Dice: {dice:.4f}")
print(f"AJI: {aji:.4f}")
print(f"PQ: {pq:.4f}, SQ: {sq:.4f}, RQ: {rq:.4f}")
```

### 评估循环

```python
from pathml.ml.metrics import evaluate_hovernet

# Comprehensive HoVer-Net evaluation
model.eval()
all_preds = []
all_targets = []

with torch.no_grad():
    for batch in test_loader:
        images = batch['image'].to(device)
        outputs = model(images)

        # Post-process predictions
        for i in range(len(images)):
            inst_pred, type_pred = hovernet_postprocess(
                outputs['np'][i],
                outputs['hv'][i],
                outputs['nc'][i]
            )
            all_preds.append({'inst': inst_pred, 'type': type_pred})
            all_targets.append({
                'inst': batch['inst_map'][i],
                'type': batch['type_map'][i]
            })

# Compute metrics
results = evaluate_hovernet(all_preds, all_targets)

print(f"Detection F1: {results['detection_f1']:.4f}")
print(f"Classification Accuracy: {results['classification_acc']:.4f}")
print(f"Panoptic Quality: {results['pq']:.4f}")
```

## ONNX 推理

使用 ONNX 部署模型进行生产推理：

### 导出到 ONNX

```python
import torch
from pathml.ml import HoVerNet

# Load trained model
model = HoVerNet(num_types=5, pretrained=True)
model.eval()

# Create dummy input
dummy_input = torch.randn(1, 3, 256, 256)

# Export to ONNX
torch.onnx.export(
    model,
    dummy_input,
    'hovernet_model.onnx',
    export_params=True,
    opset_version=11,
    input_names=['input'],
    output_names=['np_output', 'hv_output', 'nc_output'],
    dynamic_axes={
        'input': {0: 'batch_size'},
        'np_output': {0: 'batch_size'},
        'hv_output': {0: 'batch_size'},
        'nc_output': {0: 'batch_size'}
    }
)
```

### ONNX 运行时推理

```python
import onnxruntime as ort
import numpy as np

# Load ONNX model
session = ort.InferenceSession('hovernet_model.onnx')

# Prepare input
input_name = session.get_inputs()[0].name
tile_image = preprocess_tile(tile)  # Normalize, transpose to (1, 3, H, W)

# Run inference
outputs = session.run(None, {input_name: tile_image})
np_output, hv_output, nc_output = outputs

# Post-process
inst_map, type_map = hovernet_postprocess(np_output, hv_output, nc_output)
```

### 批量推理管道

```python
from pathml.core import SlideData
from pathml.preprocessing import Pipeline
import onnxruntime as ort

def run_onnx_inference_pipeline(slide_path, onnx_model_path):
    # Load slide
    wsi = SlideData.from_slide(slide_path)
    wsi.generate_tiles(level=1, tile_size=256, stride=256)

    # Load ONNX model
    session = ort.InferenceSession(onnx_model_path)
    input_name = session.get_inputs()[0].name

    # Inference on all tiles
    results = []
    for tile in wsi.tiles:
        # Preprocess
        tile_array = preprocess_tile(tile.image)

        # Inference
        outputs = session.run(None, {input_name: tile_array})

        # Post-process
        inst_map, type_map = hovernet_postprocess(*outputs)

        results.append({
            'coords': tile.coords,
            'instance_map': inst_map,
            'type_map': type_map
        })

    return results

# Run on slide
results = run_onnx_inference_pipeline('slide.svs', 'hovernet_model.onnx')
```

## 迁移学习

在自定义数据集上微调预训练模型：

```python
from pathml.ml import HoVerNet

# Load pre-trained model
model = HoVerNet(num_types=5, pretrained=True)

# Freeze encoder layers for initial training
for name, param in model.named_parameters():
    if 'encoder' in name:
        param.requires_grad = False

# Fine-tune only decoder and classification heads
optimizer = torch.optim.Adam(
    filter(lambda p: p.requires_grad, model.parameters()),
    lr=1e-4
)

# Train for a few epochs
train_for_n_epochs(model, train_loader, optimizer, num_epochs=10)

# Unfreeze all layers for full fine-tuning
for param in model.parameters():
    param.requires_grad = True

# Continue training with lower learning rate
optimizer = torch.optim.Adam(model.parameters(), lr=1e-5)
train_for_n_epochs(model, train_loader, optimizer, num_epochs=50)
```

## 最佳实践

1. **使用预先训练的模型（如果可用）：**
   - 从 pretrained=True 开始以获得更好的初始化
   - 对特定领域的数据进行微调

2. **应用适当的数据增强：**
   - 旋转、翻转以保持方向不变
   - 颜色抖动以处理染色变化
   - 生物变异性的弹性变形

3. **监控多个指标：**
   - 分别进行轨迹检测、分割和分类
   - 使用超出标准精度的特定领域指标（AJI、PQ）

4. **处理类别不平衡：**
   - 稀有细胞类型的加权损失函数
   - 对少数群体进行过采样
   - 困难示例的焦点损失

5. **在不同的组织类型上进行验证：**
   - 确保跨不同组织的泛化
   - 对保留的解剖部位进行测试

6. **优化推理：**
   - 导出到 ONNX 以加快部署速度
   - 批量切片可有效利用 GPU
   - 尽可能使用混合精度 (FP16)

7. **定期保存检查点：**
   - 根据验证指标保留最佳模型
   - 保存优化器状态以恢复训练

## 常见问题及解决方案

**问题：细胞核边界分割不良**
- 使用 HV 图（水平-垂直）来分离接触的原子核
- 增加 HV 损失项的权重
- 应用形态学后处理

**问题：相似细胞类型的错误分类**
- 增加分类损失权重
- 添加层次分类（HACTNet）
- 增加困惑类别的训练数据

**问题：训练不稳定或不收敛**
- 降低学习率
- 使用渐变裁剪：`torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)`
- 检查数据预处理问题

**问题：训练期间内存不足**
- 减少批量
- 使用梯度累积
- 启用混合精度训练：`torch.cuda.amp`

**问题：模型与训练数据过度拟合**
- 增加数据增强
- 添加滤除层
- 减少模型容量
- 根据验证损失使用提前停止

## 其他资源

- **PathML 机器学习 API：** https://pathml.readthedocs.io/en/latest/api_ml_reference.html
- **HoVer-Net 论文：** Graham 等人，“HoVer-Net：多组织组织学图像中细胞核的同步分割和分类”，医学图像分析，2019 年
- **PanNuke 数据集：** https://warwick.ac.uk/fac/cross_fac/tia/data/pannuke
- **PyTorch 闪电：** https://www.pytorchlightning.ai/
- **ONNX 运行时：** https://onnxruntime.ai/