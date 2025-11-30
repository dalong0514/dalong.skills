<!-- 此文件由机器翻译自 preprocessing.md -->

# 预处理管道和转换

## 概述

PathML 提供了基于组织成管道的可组合转换的模块化预处理架构。变换是修改图像、创建蒙版或提取特征的单独操作。管道链一起转换，为计算病理学创建可重复、可扩展的预处理工作流程。

## 管道架构

### 管道类

`Pipeline` 类组成了一系列连续应用的变换：

```python
from pathml.preprocessing import Pipeline, Transform1, Transform2

# Create pipeline
pipeline = Pipeline([
    Transform1(param1=value1),
    Transform2(param2=value2),
    # ... more transforms
])

# Run on a single slide
pipeline.run(slide_data)

# Run on a dataset
pipeline.run(dataset, distributed=True, n_workers=8)
```

**主要特点：**
- 转换的顺序执行
- 自动处理瓷砖和掩模
- Dask 的分布式处理支持
- 具有可序列化配置的可重复工作流程

### 转换基类

所有转换均继承自 `Transform` 基类并实现：
- `apply()` - 核心转换逻辑
- `input_type` - 预期输入（图块、掩码等）
- `output_type` - 生成的输出

## 变换类别

PathML 提供六大类的转换：

1. **图像修改** - 模糊、重新缩放、直方图均衡
2. **掩模创建** - 组织检测、细胞核检测、阈值处理
3. **Mask Modification** - 对掩模进行形态学操作
4. **染色处理** - H&E 染色标准化和分离
5. **质量控制** - 伪影检测、空白标记
6. **专业** - 多参数成像、细胞分割

## 图像修改变换

### 模糊操作

应用各种模糊内核来降噪：

**中值模糊：**
<<<代码块_1>>>
- 对椒盐噪音有效
- 比高斯模糊更好地保留边缘

**高斯模糊：**
<<<代码块_2>>>
- 平滑降噪
- 可调节西格玛控制模糊强度

**框模糊：**
<<<代码块_3>>>
- 最快的模糊操作
- 内核内统一平均

### 强度调整

**重新调整强度：**
<<<代码块_4>>>

**直方图均衡：**
<<<代码块_5>>>
- 增强全局对比度
- 分散强度分布

**自适应直方图均衡化 (CLAHE)：**
<<<代码块_6>>>
- 增强局部对比度
- 使用clip_limit防止过度放大
- 更适合具有不同局部对比度的图像

### 超像素处理

**超像素插值：**
```python
from pathml.preprocessing import SuperpixelInterpolation

# Divide into superpixels using SLIC
transform = SuperpixelInterpolation(
    n_segments=100,
    compactness=10.0
)
```
- 将图像分割成感知上有意义的区域
- 对于特征提取和分割有用

## 面具创作转变

### H&E 组织和细胞核检测

**组织检测HE：**
```python
from pathml.preprocessing import TissueDetectionHE

# Detect tissue regions in H&E slides
transform = TissueDetectionHE(
    use_saturation=True,  # Use HSV saturation channel
    threshold=10,  # Intensity threshold
    min_region_size=500  # Minimum tissue region size in pixels
)
```
- 创建二元组织面膜
- 过滤小区域和伪影
- 将掩码存储在 `tile.masks['tissue']` 中

**细胞核检测HE：**
```python
from pathml.preprocessing import NucleusDetectionHE

# Detect nuclei in H&E images
transform = NucleusDetectionHE(
    stain='hematoxylin',  # Use hematoxylin channel
    threshold=0.3,
    min_nucleus_size=10
)
```
- 分离苏木精染色
- 创建核掩模的阈值
- 将掩码存储在 `tile.masks['nucleus']` 中

### 二进制阈值

**二进制阈值：**
```python
from pathml.preprocessing import BinaryThreshold

# Threshold using Otsu's method
transform = BinaryThreshold(
    method='otsu',  # 'otsu' or manual threshold value
    invert=False
)

# Or specify manual threshold
transform = BinaryThreshold(threshold=128)
```

### 前景检测

**前景检测：**
```python
from pathml.preprocessing import ForegroundDetection

# Detect foreground regions
transform = ForegroundDetection(
    threshold=0.5,
    min_region_size=1000,  # Minimum size in pixels
    use_saturation=True
)
```

## 蒙版修改变换

应用形态学操作来清理掩模：

**变形开放：**
```python
from pathml.preprocessing import MorphOpen

# Remove small objects and noise
transform = MorphOpen(
    kernel_size=5,
    mask_name='tissue'  # Which mask to modify
)
```
- 腐蚀后膨胀
- 去除小物体和噪音

**变形关闭：**
```python
from pathml.preprocessing import MorphClose

# Fill small holes
transform = MorphClose(
    kernel_size=5,
    mask_name='tissue'
)
```
- 膨胀后腐蚀
- 填充面罩上的小孔

## 染色标准化

### 染色标准化HE

标准化载玻片上的 H&E 染色，以考虑染色程序和扫描仪的变化：

```python
from pathml.preprocessing import StainNormalizationHE

# Normalize to reference slide
transform = StainNormalizationHE(
    target='normalize',  # 'normalize', 'hematoxylin', or 'eosin'
    stain_estimation_method='macenko',  # 'macenko' or 'vahadane'
    tissue_mask_name=None  # Optional tissue mask for better estimation
)
```

**目标模式：**
- `'normalize'` - 将两个污渍标准化为参考
- `'hematoxylin'` - 仅提取苏木精通道
- `'eosin'` - 仅提取曙红通道

**污点估计方法：**
- `'macenko'` - Macenko 等人。 2009方法（更快、更稳定）
- `'vahadane'` - Vahadane 等人。 2016年方法（更准确，更慢）

**高级参数：**
```python
transform = StainNormalizationHE(
    target='normalize',
    stain_estimation_method='macenko',
    target_od=None,  # Optical density matrix for reference (optional)
    target_concentrations=None,  # Target stain concentrations (optional)
    regularizer=0.1,  # Regularization for vahadane method
    background_intensity=240  # Background intensity level
)
```

**工作流程：**
1. 将RGB转换为光密度(OD)
2. 估计染色矩阵（H&E 向量）
3.分解成污渍浓度
4. 标准化参考染色分布
5. 重建归一化RGB图像

**使用纸巾面膜的示例：**
```python
from pathml.preprocessing import Pipeline, TissueDetectionHE, StainNormalizationHE

pipeline = Pipeline([
    TissueDetectionHE(),  # Create tissue mask first
    StainNormalizationHE(
        target='normalize',
        stain_estimation_method='macenko',
        tissue_mask_name='tissue'  # Use tissue mask for better estimation
    )
])
```

## 质量控制转变

### 伪影检测

**标签ArtifactTileHE：**
```python
from pathml.preprocessing import LabelArtifactTileHE

# Label tiles containing artifacts
transform = LabelArtifactTileHE(
    pen_threshold=0.5,  # Threshold for pen marking detection
    bubble_threshold=0.5  # Threshold for bubble detection
)
```
- 检测钢笔标记、气泡和其他伪影
- 标记受影响的图块以进行过滤

**标签空白HE：**
```python
from pathml.preprocessing import LabelWhiteSpaceHE

# Label tiles with excessive white space
transform = LabelWhiteSpaceHE(
    threshold=0.9,  # Fraction of white pixels
    mask_name='white_space'
)
```
- 识别主要是背景的瓷砖
- 可用于过滤无信息的图块

## 多参数成像变换

### 细胞分割

**分段MIF：**
```python
from pathml.preprocessing import SegmentMIF

# Segment cells using Mesmer deep learning model
transform = SegmentMIF(
    nuclear_channel='DAPI',  # Nuclear marker channel name
    cytoplasm_channel='CD45',  # Cytoplasm marker channel name
    model='mesmer',  # Deep learning segmentation model
    image_resolution=0.5,  # Microns per pixel
    compartment='whole-cell'  # 'nuclear', 'cytoplasm', or 'whole-cell'
)
```
- 使用 DeepCell Mesmer 模型进行细胞分割
- 需要核和细胞质通道规范
- 生成实例分割掩码

**分段MIF远程：**
```python
from pathml.preprocessing import SegmentMIFRemote

# Remote inference using DeepCell API
transform = SegmentMIFRemote(
    nuclear_channel='DAPI',
    cytoplasm_channel='CD45',
    model='mesmer',
    api_url='https://deepcell.org/api'
)
```
- 与 SegmentMIF 功能相同，但使用远程 API
- 无需本地 GPU
- 适合批量处理

### 标记定量

**量化MIF：**
```python
from pathml.preprocessing import QuantifyMIF

# Quantify marker expression per cell
transform = QuantifyMIF(
    segmentation_mask_name='cell_segmentation',
    markers=['CD3', 'CD4', 'CD8', 'CD20', 'CD45'],
    output_format='anndata'  # or 'dataframe'
)
```
- 提取每个分段细胞的平均标记强度
- 计算形态特征（面积、周长等）
- 输出 AnnData 对象用于下游单细胞分析

### CODEX/Vectra 特定

**崩溃运行CODEX：**
```python
from pathml.preprocessing import CollapseRunsCODEX

# Consolidate multi-run CODEX data
transform = CollapseRunsCODEX(
    z_slice=2,  # Select specific z-slice
    run_order=[0, 1, 2]  # Order of acquisition runs
)
```
- 合并来自多个 CODEX 采集运行的通道
- 从 z 堆栈中选择焦平面

**CollapseRunsVectra：**
```python
from pathml.preprocessing import CollapseRunsVectra

# Process Vectra multiplex IF data
transform = CollapseRunsVectra(
    wavelengths=[520, 570, 620, 670, 780]  # Emission wavelengths
)
```

## 建设综合管道

### 基本 H&E 预处理管道

```python
from pathml.preprocessing import (
    Pipeline,
    TissueDetectionHE,
    StainNormalizationHE,
    NucleusDetectionHE,
    MedianBlur,
    LabelWhiteSpaceHE
)

pipeline = Pipeline([
    # 1. Quality control
    LabelWhiteSpaceHE(threshold=0.9),

    # 2. Noise reduction
    MedianBlur(kernel_size=3),

    # 3. Tissue detection
    TissueDetectionHE(min_region_size=500),

    # 4. Stain normalization
    StainNormalizationHE(
        target='normalize',
        stain_estimation_method='macenko',
        tissue_mask_name='tissue'
    ),

    # 5. Nucleus detection
    NucleusDetectionHE(threshold=0.3)
])
```

### CODEX 多参数管道

```python
from pathml.preprocessing import (
    Pipeline,
    CollapseRunsCODEX,
    SegmentMIF,
    QuantifyMIF
)

codex_pipeline = Pipeline([
    # 1. Consolidate multi-run data
    CollapseRunsCODEX(z_slice=2),

    # 2. Cell segmentation
    SegmentMIF(
        nuclear_channel='DAPI',
        cytoplasm_channel='CD45',
        model='mesmer',
        image_resolution=0.377
    ),

    # 3. Quantify markers
    QuantifyMIF(
        segmentation_mask_name='cell_segmentation',
        markers=['CD3', 'CD4', 'CD8', 'CD20', 'PD1', 'PDL1'],
        output_format='anndata'
    )
])
```

### 具有质量控制的先进管道

```python
from pathml.preprocessing import (
    Pipeline,
    LabelWhiteSpaceHE,
    LabelArtifactTileHE,
    TissueDetectionHE,
    MorphOpen,
    MorphClose,
    StainNormalizationHE,
    AdaptiveHistogramEqualization
)

advanced_pipeline = Pipeline([
    # Stage 1: Quality control
    LabelWhiteSpaceHE(threshold=0.85),
    LabelArtifactTileHE(pen_threshold=0.5, bubble_threshold=0.5),

    # Stage 2: Tissue detection
    TissueDetectionHE(threshold=10, min_region_size=1000),
    MorphOpen(kernel_size=5, mask_name='tissue'),
    MorphClose(kernel_size=7, mask_name='tissue'),

    # Stage 3: Stain normalization
    StainNormalizationHE(
        target='normalize',
        stain_estimation_method='vahadane',
        tissue_mask_name='tissue'
    ),

    # Stage 4: Contrast enhancement
    AdaptiveHistogramEqualization(clip_limit=0.03, tile_grid_size=(8, 8))
])
```

## 运行管道

### 单张幻灯片处理

```python
from pathml.core import SlideData

# Load slide
wsi = SlideData.from_slide("slide.svs")

# Generate tiles
wsi.generate_tiles(level=1, tile_size=256, stride=256)

# Run pipeline
pipeline.run(wsi)

# Access processed data
for tile in wsi.tiles:
    normalized_image = tile.image
    tissue_mask = tile.masks.get('tissue')
    nucleus_mask = tile.masks.get('nucleus')
```

### 分布式执行的批处理

```python
from pathml.core import SlideDataset
from dask.distributed import Client
import glob

# Start Dask client
client = Client(n_workers=8, threads_per_worker=2, memory_limit='4GB')

# Create dataset
slide_paths = glob.glob("data/*.svs")
dataset = SlideDataset(
    slide_paths,
    tile_size=512,
    stride=512,
    level=1
)

# Run pipeline in parallel
dataset.run(
    pipeline,
    distributed=True,
    client=client
)

# Save results
dataset.to_hdf5("processed_dataset.h5")

client.close()
```

### 条件管道执行

仅对满足特定条件的图块执行变换：

```python
# Filter tiles before processing
wsi.generate_tiles(level=1, tile_size=256)

# Run pipeline only on tissue tiles
for tile in wsi.tiles:
    if tile.masks.get('tissue') is not None:
        pipeline.run(tile)
```

## 性能优化

### 内存管理

```python
# Process large datasets in batches
batch_size = 100
for i in range(0, len(slide_paths), batch_size):
    batch_paths = slide_paths[i:i+batch_size]
    batch_dataset = SlideDataset(batch_paths)
    batch_dataset.run(pipeline, distributed=True)
    batch_dataset.to_hdf5(f"batch_{i}.h5")
```

### GPU 加速

某些转换会利用 GPU 加速（如果可用）：

```python
import torch

# Check GPU availability
print(f"CUDA available: {torch.cuda.is_available()}")

# Transforms that benefit from GPU:
# - SegmentMIF (Mesmer deep learning model)
# - StainNormalizationHE (matrix operations)
```

### 并行工作线程配置

```python
from dask.distributed import Client

# CPU-bound tasks (image processing)
client = Client(
    n_workers=8,
    threads_per_worker=1,  # Use processes, not threads
    memory_limit='8GB'
)

# GPU tasks (deep learning inference)
client = Client(
    n_workers=2,  # Fewer workers for GPU
    threads_per_worker=4,
    processes=True
)
```

## 自定义转换

通过子类化 `Transform` 创建自定义预处理操作：

```python
from pathml.preprocessing.transforms import Transform
import numpy as np

class CustomTransform(Transform):
    def __init__(self, param1, param2):
        self.param1 = param1
        self.param2 = param2

    def apply(self, tile):
        # Access tile image
        image = tile.image

        # Apply custom operation
        processed = self.custom_operation(image, self.param1, self.param2)

        # Update tile
        tile.image = processed

        return tile

    def custom_operation(self, image, param1, param2):
        # Implement custom logic
        return processed_image

# Use in pipeline
pipeline = Pipeline([
    CustomTransform(param1=10, param2=0.5),
    # ... other transforms
])
```

## 最佳实践

1. **顺序适当变换：**
   - 质量控制第一（LabelWhiteSpace、LabelArtifact）
   - 早期降噪（模糊）
   - 染色标准化前的组织检测
   - 在颜色相关操作之前进行染色归一化

2. **使用组织掩模进行染色标准化：**
   - 通过排除背景提高准确性
   - `TissueDetectionHE()` 然后`StainNormalizationHE(tissue_mask_name='tissue')`

3. **应用形态学操作来清洁掩模：**
   - `MorphOpen` 删除小的误报
   - `MorphClose` 填补小空白

4. **利用大型数据集的分布式处理：**
   - 使用Dask进行并行执行
   - 根据可用资源配置工作人员

5. **保存中间结果：**
   - 将处理后的数据存储到HDF5以供重复使用
   - 避免重新处理计算量大的转换

6. **验证样本图像的预处理：**
   - 可视化中间步骤
   - 在批量处理之前调整代表性样品的参数

7. **处理边缘情况：**
   - 在下游操作之前检查空掩模
   - 在进行昂贵的计算之前验证图块质量

## 常见问题及解决方案

**问题：染色归一化会产生伪影**
- 使用组织遮罩排除背景
- 尝试不同的污渍估计方法（macenko 与 vahadane）
- 验证光密度参数与您的图像匹配

**问题：管道执行期间内存不足**
- 减少 Dask 工作人员数量
- 减少瓷砖尺寸
- 在较低金字塔级别处理图像
- 在 Dask 客户端中启用内存限制参数

**问题：组织检测遗漏组织区域**
- 调整阈值参数
- 使用饱和通道：`use_saturation=True`
- 减少 min_region_size 以捕获较小的组织碎片

**问题：细胞核检测不准确**
- 验证染色分离质量（可视化苏木精通道）
- 调整阈值参数
- 在细胞核检测之前应用染色归一化

## 其他资源

- **PathML 预处理 API：** https://pathml.readthedocs.io/en/latest/api_preprocessing_reference.html
- **染色归一化方法：**
  - 马琴科等人。 2009：“一种用于定量分析的标准化组织学载玻片的方法”
  - 瓦哈丹等人。 2016：“结构保持颜色标准化和稀疏污点分离”
- **DeepCell Mesmer：** https://www.deepcell.org/（细胞分割模型）