<!-- 此文件由机器翻译自 slide_management.md -->

# 幻灯片管理

## 概述

`Slide` 类是在 histolab 中处理整个幻灯片图像 (WSI) 的主要接口。它提供了加载、检查和处理以各种格式存储的大型组织病理学图像的方法。

## 初始化

```python
from histolab.slide import Slide

# Initialize a slide with a WSI file and output directory
slide = Slide(processed_path="path/to/processed/output",
              slide_path="path/to/slide.svs")
```

**参数：**
- `slide_path`：整个幻灯片图像文件的路径（支持多种格式：SVS、TIFF、NDPI 等）
- `processed_path`：保存已处理输出（图块、缩略图等）的目录

## 加载样本数据

Histolab 提供来自 TCGA 的内置样本数据集用于测试和演示：

<<<代码块_1>>>

可用的示例数据集：
- `prostate_tissue()`：前列腺组织样本
- `ovarian_tissue()`：卵巢组织样本
- `breast_tissue()`：乳房组织样本
- `heart_tissue()`：心脏组织样本
- `kidney_tissue()`：肾脏组织样本

## 关键属性

### 幻灯片尺寸
<<<代码块_2>>>

### 放大倍率信息
<<<代码块_3>>>

### 幻灯片属性
<<<代码块_4>>>

## 缩略图生成

<<<代码块_5>>>

## 幻灯片可视化

<<<代码块_6>>>

## 提取区域

```python
# Extract region at specific coordinates and level
region = slide.extract_region(
    location=(x, y),  # Top-left coordinates at level 0
    size=(width, height),  # Region size
    level=0  # Pyramid level
)
```

## 使用金字塔级别

WSI 文件使用具有多个分辨率级别的金字塔结构：
- 级别 0：最高分辨率（原始扫描分辨率）
- 1+ 级：逐渐降低分辨率以加快访问速度

```python
# Check available levels
for level in range(slide.levels):
    dims = slide.level_dimensions[level]
    downsample = slide.level_downsamples[level]
    print(f"Level {level}: {dims}, downsample: {downsample}x")
```

## 幻灯片名称和路径

```python
# Get slide filename without extension
slide_name = slide.name

# Get full path to slide file
slide_path = slide.scaled_image
```

## 最佳实践

1. **始终指定processed_path**：在专用目录中组织输出
2. **处理前检查尺寸**：大型幻灯片可能会超出内存限制
3. **使用适当的金字塔级别**：在与您的分析分辨率匹配的级别提取图块
4. **使用缩略图预览**：在进行繁重处理之前使用缩略图进行快速可视化
5. **监控内存使用情况**：大型幻灯片上的 0 级操作需要大量 RAM

## 常见工作流程

### 载玻片检查工作流程
```python
from histolab.slide import Slide

# Load slide
slide = Slide("slide.svs", processed_path="output/")

# Inspect properties
print(f"Dimensions: {slide.dimensions}")
print(f"Levels: {slide.levels}")
print(f"Magnification: {slide.properties.get('openslide.objective-power', 'N/A')}")

# Save thumbnail for review
slide.save_thumbnail()
```

### 多幻灯片处理
```python
import os
from pathlib import Path

slide_dir = Path("slides/")
output_dir = Path("processed/")

for slide_path in slide_dir.glob("*.svs"):
    slide = Slide(slide_path, processed_path=output_dir / slide_path.stem)
    # Process each slide
    print(f"Processing: {slide.name}")
```