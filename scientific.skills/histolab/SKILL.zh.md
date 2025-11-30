<!-- 此文件由机器翻译自 SKILL.md -->

---
名称： 组织实验室
描述：用于整个幻灯片图像（WSI）的数字病理图像处理工具包。在使用组织病理学载玻片、处理 H&E 或 IHC 染色的组织图像、从十亿像素病理图像中提取图块、检测组织区域、分割组织掩模或为计算病理学深度学习管道准备数据集时，请使用此技能。适用于 WSI 格式（SVS、TIFF、NDPI）、基于切片的分析和组织学图像预处理工作流程。
---

# 组织实验室

## 概述

Histolab 是一个用于处理数字病理学中的整个幻灯片图像 (WSI) 的 Python 库。它可以自动执行组织检测，从十亿像素图像中提取信息图块，并为深度学习管道准备数据集。该库可处理多种 WSI 格式，实现复杂的组织分割，并提供灵活的切片提取策略。

## 安装

```bash
uv pip install histolab
```

## 快速入门

从整个幻灯片图像中提取图块的基本工作流程：

<<<代码块_1>>>

## 核心能力

### 1.幻灯片管理

加载、检查和处理各种格式的整个幻灯片图像。

**常用操作：**
- 加载WSI文件（SVS、TIFF、NDPI等）
- 访问幻灯片元数据（尺寸、放大倍率、属性）
- 生成缩略图以进行可视化
- 使用金字塔图像结构
- 提取特定坐标处的区域

**关键类：** `Slide`

**参考：** `references/slide_management.md` 包含以下方面的综合文档：
- 幻灯片初始化和配置
- 内置样本数据集（前列腺、卵巢、乳房、心脏、肾脏组织）
- 访问幻灯片属性和元数据
- 缩略图生成和可视化
- 使用金字塔级别
- 多幻灯片处理工作流程

**工作流程示例：**
<<<代码块_2>>>

### 2. 组织检测和面罩

自动识别组织区域并过滤背景/伪影。

**常用操作：**
- 创建二元组织掩模
- 检测最大的组织区域
- 排除背景和伪影
- 定制组织分割
- 删除笔注释

**关键类：** `TissueMask`、`BiggestTissueBoxMask`、`BinaryMask`

**参考：** `references/tissue_masks.md` 包含以下内容的综合文档：
- TissueMask：使用自动过滤器分割所有组织区域
- BiggestTissueBoxMask：返回最大组织区域的边界框（默认）
- BinaryMask：自定义掩码实现的基类
- 使用 `locate_mask()` 可视化蒙版
- 创建自定义矩形和注释排除蒙版
- 掩模与图块提取集成
- 最佳实践和故障排除

**工作流程示例：**
<<<代码块_3>>>

**何时使用每个面膜：**
- `TissueMask`：多个组织切片，综合分析
- `BiggestTissueBoxMask`：单个主要组织切片，排除伪影（默认）
- 自定义 `BinaryMask`：特定 ROI、排除注释、自定义分段

### 3. 瓦片提取

使用不同的策略从大型 WSI 中提取较小的区域。

**三种提取策略：**

**RandomTiler：** 提取固定数量的随机定位的图块
- 最适合：不同区域采样、探索性分析、训练数据
- 关键参数：`n_tiles`、`seed` 用于再现性

**GridTiler：** 以网格图案系统地跨组织提取图块
- 最适合：完整覆盖、空间分析、重建
- 关键参数：`pixel_overlap` 用于滑动窗口

**ScoreTiler：** 根据评分函数提取排名靠前的图块
- 最适合：信息最丰富的地区，质量驱动的选择
- 关键参数：`scorer`（NucleiScorer、CellularityScorer、自定义）

**常用参数：**
- `tile_size`：平铺尺寸（例如，(512, 512)）
- `level`：提取的金字塔级别（0 = 最高分辨率）
- `check_tissue`：按组织内容过滤图块
- `tissue_percent`：最小组织覆盖率（默认 80%）
- `extraction_mask`：定义提取区域的掩码

**参考：** `references/tile_extraction.md` 包含以下内容的综合文档：
- 每个铺砖策略的详细解释
- 可用的记分器（NucleiScorer、CellularityScorer、自定义）
- 使用 `locate_tiles()` 进行平铺预览
- 提取工作流程和报告
- 高级模式（多级、分层提取）
- 性能优化和故障排除

**工作流程示例：**

<<<代码块_4>>>
**提取前始终预览：**
<<<代码块_5>>>

### 4. 过滤器和预处理

应用图像处理滤波器进行组织检测、质量控制和预处理。

**过滤类别：**

**图像过滤器：** 色彩空间转换、阈值处理、对比度增强
- `RgbToGrayscale`、`RgbToHsv`、`RgbToHed`
- `OtsuThreshold`、`AdaptiveThreshold`
- `StretchContrast`、`HistogramEqualization`

**形态过滤器：** 二值图像上的结构操作
- `BinaryDilation`、`BinaryErosion`
- `BinaryOpening`、`BinaryClosing`
- `RemoveSmallObjects`、`RemoveSmallHoles`

**组成：** 将多个过滤器链接在一起
- `Compose`：创建过滤器管道

**参考：** `references/filters_preprocessing.md` 包含以下内容的综合文档：
- 每种过滤器类型的详细解释
- 过滤器组成和链接
- 常见的预处理流程（组织检测、笔去除、细胞核增强）
- 对瓷砖应用过滤器
- 定制面罩过滤器
- 质量控制滤镜（模糊检测、组织覆盖）
- 最佳实践和故障排除

**工作流程示例：**

<<<代码块_6>>>

### 5.可视化

可视化幻灯片、蒙版、图块位置和提取质量。

**常见可视化任务：**
- 显示幻灯片缩略图
- 可视化组织面罩
- 预览图块位置
- 评估瓷砖质量
- 创建报告和数据

**参考：** `references/visualization.md` 包含以下内容的综合文档：
- 幻灯片缩略图显示和保存
- 使用 `locate_mask()` 屏蔽可视化
- 使用 `locate_tiles()` 进行平铺位置预览
- 显示提取的瓷砖和马赛克
- 质量评估（分数分布、顶部与底部瓷砖）
- 多幻灯片可视化
- 滤镜效果可视化
- 导出高分辨率数据和 PDF 报告
- Jupyter 笔记本中的交互式可视化

**工作流程示例：**

```python
import matplotlib.pyplot as plt
from histolab.masks import TissueMask

# Display slide thumbnail
plt.figure(figsize=(10, 10))
plt.imshow(slide.thumbnail)
plt.title(f"Slide: {slide.name}")
plt.axis('off')
plt.show()

# Visualize tissue mask
tissue_mask = TissueMask()
slide.locate_mask(tissue_mask)

# Preview tile locations
tiler = RandomTiler(tile_size=(512, 512), n_tiles=50)
tiler.locate_tiles(slide, n_tiles=20)

# Display extracted tiles in grid
from pathlib import Path
from PIL import Image

tile_paths = list(Path("output/tiles/").glob("*.png"))[:16]
fig, axes = plt.subplots(4, 4, figsize=(12, 12))
axes = axes.ravel()

for idx, tile_path in enumerate(tile_paths):
    tile_img = Image.open(tile_path)
    axes[idx].imshow(tile_img)
    axes[idx].set_title(tile_path.stem, fontsize=8)
    axes[idx].axis('off')

plt.tight_layout()
plt.show()
```

## 典型工作流程

### 工作流程 1：探索性图块提取

对不同组织区域进行快速采样以进行初步分析。

```python
from histolab.slide import Slide
from histolab.tiler import RandomTiler
import logging

# Enable logging for progress tracking
logging.basicConfig(level=logging.INFO)

# Load slide
slide = Slide("slide.svs", processed_path="output/random_tiles/")

# Inspect slide
print(f"Dimensions: {slide.dimensions}")
print(f"Levels: {slide.levels}")
slide.save_thumbnail()

# Configure random tiler
random_tiler = RandomTiler(
    tile_size=(512, 512),
    n_tiles=100,
    level=0,
    seed=42,
    check_tissue=True,
    tissue_percent=80.0
)

# Preview locations
random_tiler.locate_tiles(slide, n_tiles=20)

# Extract tiles
random_tiler.extract(slide)
```

### 工作流程2：综合网格提取

完整的组织覆盖，用于整个载玻片分析。

```python
from histolab.slide import Slide
from histolab.tiler import GridTiler
from histolab.masks import TissueMask

# Load slide
slide = Slide("slide.svs", processed_path="output/grid_tiles/")

# Use TissueMask for all tissue sections
tissue_mask = TissueMask()
slide.locate_mask(tissue_mask)

# Configure grid tiler
grid_tiler = GridTiler(
    tile_size=(512, 512),
    level=1,  # Use level 1 for faster extraction
    pixel_overlap=0,
    check_tissue=True,
    tissue_percent=70.0
)

# Preview grid
grid_tiler.locate_tiles(slide)

# Extract all tiles
grid_tiler.extract(slide, extraction_mask=tissue_mask)
```

### 工作流程 3：以质量为导向的瓷砖选择

根据原子核密度提取最有信息的图块。

```python
from histolab.slide import Slide
from histolab.tiler import ScoreTiler
from histolab.scorer import NucleiScorer
import pandas as pd
import matplotlib.pyplot as plt

# Load slide
slide = Slide("slide.svs", processed_path="output/scored_tiles/")

# Configure score tiler
score_tiler = ScoreTiler(
    tile_size=(512, 512),
    n_tiles=50,
    level=0,
    scorer=NucleiScorer(),
    check_tissue=True
)

# Preview top tiles
score_tiler.locate_tiles(slide, n_tiles=15)

# Extract with report
score_tiler.extract(slide, report_path="tiles_report.csv")

# Analyze scores
report_df = pd.read_csv("tiles_report.csv")
plt.hist(report_df['score'], bins=20, edgecolor='black')
plt.xlabel('Tile Score')
plt.ylabel('Frequency')
plt.title('Distribution of Tile Scores')
plt.show()
```

### 工作流程 4：多载玻片处理流程

使用一致的参数处理整个幻灯片集合。

```python
from pathlib import Path
from histolab.slide import Slide
from histolab.tiler import RandomTiler
import logging

logging.basicConfig(level=logging.INFO)

# Configure tiler once
tiler = RandomTiler(
    tile_size=(512, 512),
    n_tiles=50,
    level=0,
    seed=42,
    check_tissue=True
)

# Process all slides
slide_dir = Path("slides/")
output_base = Path("output/")

for slide_path in slide_dir.glob("*.svs"):
    print(f"\nProcessing: {slide_path.name}")

    # Create slide-specific output directory
    output_dir = output_base / slide_path.stem
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load and process slide
    slide = Slide(slide_path, processed_path=output_dir)

    # Save thumbnail for review
    slide.save_thumbnail()

    # Extract tiles
    tiler.extract(slide)

    print(f"Completed: {slide_path.name}")
```

### 工作流程 5：自定义组织检测和过滤

处理带有伪影、注释或异常染色的载玻片。

```python
from histolab.slide import Slide
from histolab.masks import TissueMask
from histolab.tiler import RandomTiler
from histolab.filters.compositions import Compose
from histolab.filters.image_filters import RgbToGrayscale, OtsuThreshold
from histolab.filters.morphological_filters import (
    BinaryDilation, RemoveSmallObjects, RemoveSmallHoles
)

# Define custom filter pipeline for aggressive artifact removal
aggressive_filters = Compose([
    RgbToGrayscale(),
    OtsuThreshold(),
    BinaryDilation(disk_size=10),
    RemoveSmallHoles(area_threshold=5000),
    RemoveSmallObjects(area_threshold=3000)  # Remove larger artifacts
])

# Create custom mask
custom_mask = TissueMask(filters=aggressive_filters)

# Load slide and visualize mask
slide = Slide("slide.svs", processed_path="output/")
slide.locate_mask(custom_mask)

# Extract with custom mask
tiler = RandomTiler(tile_size=(512, 512), n_tiles=100)
tiler.extract(slide, extraction_mask=custom_mask)
```

## 最佳实践

### 载玻片装载和检查
1. 加工前务必检查载玻片性能
2. 保存缩略图以便快速查看
3. 检查金字塔级别和尺寸
4. 使用缩略图验证组织是否存在

### 组织检测
1. 在提取之前使用 `locate_mask()` 预览蒙版
2. 对多个部分使用 `TissueMask`，对单个部分使用 `BiggestTissueBoxMask`
3. 针对特定染色定制过滤器（H&E 与 IHC）
4. 使用自定义蒙版处理笔注释
5. 在不同的载玻片上测试掩模

### 瓦片提取
1. **在提取之前始终使用 `locate_tiles()` 进行预览**
2.选择合适的贴砖机：
   - RandomTiler：采样和探索
   - GridTiler：完整覆盖
   - ScoreTiler：质量驱动的选择
3. 设置适当的 `tissue_percent` 阈值（典型值为 70-90%）
4. 在 RandomTiler 中使用种子来实现再现性
5. 在适当的金字塔级别提取以进行分析分辨率
6. 为大型数据集启用日志记录

### 性能
1. 在较低级别（1、2）提取以加快处理速度
2. 在适当的时候使用 `BiggestTissueBoxMask` 而不是 `TissueMask`
3.调整`tissue_percent`以减少无效平铺尝试
4. 限制 `n_tiles` 进行初始探索
5. 对于不重叠的网格使用 `pixel_overlap=0`

### 质量控制
1. 验证图块质量（检查模糊、伪像、焦点）
2. 查看 ScoreTiler 的分数分布
3. 检查顶部和底部的计分块
4. 监测组织覆盖率统计数据
5. 如果需要，通过附加质量指标过滤提取的图块

## 常见用例

### 训练深度学习模型
- 使用 RandomTiler 跨多张幻灯片提取平衡数据集
- 使用 ScoreTiler 和 NucleiScorer 来关注细胞丰富的区域
- 以一致的分辨率提取（0 级或 1 级）
- 生成 CSV 报告以跟踪图块元数据
### 整个幻灯片分析
- 使用 GridTiler 实现完整的组织覆盖
- 在多个金字塔级别提取以进行层次分析
- 维护与网格位置的空间关系
- 使用 `pixel_overlap` 进行滑动窗口方法

### 组织表征
- 使用 RandomTiler 对不同区域进行采样
- 量化面罩的组织覆盖范围
- 通过 HED 分解提取特定于染色的信息
- 比较载玻片上的组织模式

### 质量评估
- 使用 ScoreTiler 确定最佳焦点区域
- 使用自定义掩模和过滤器检测伪影
- 评估整个载玻片集合的染色质量
- 标记有问题的幻灯片以供手动审核

### 数据集管理
- 使用 ScoreTiler 确定信息图块的优先级
- 按组织百分比过滤瓷砖
- 生成包含图块分数和元数据的报告
- 创建跨载玻片和组织类型的分层数据集

## 故障排除

### 未提取任何图块
- 降低`tissue_percent`阈值
- 验证载玻片是否包含组织（检查缩略图）
- 确保 extract_mask 捕获组织区域
- 检查tile_size是否适合幻灯片分辨率

### 许多背景瓷砖
- 启用`check_tissue=True`
- 增加`tissue_percent`阈值
- 使用适当的面罩（TissueMask 与 BiggestTissueBoxMask）
- 定制面罩过滤器以更好地检测组织

### 提取非常慢
- 在较低金字塔级别提取（级别=1或2）
- 减少 RandomTiler/ScoreTiler 的 `n_tiles`
- 使用RandomTiler代替GridTiler进行采样
- 使用 BiggestTissueBoxMask 代替 TissueMask

### 瓷砖有伪影
- 实现自定义注释排除掩码
- 调整滤波器参数以消除伪影
- 增加小物体移除阈值
- 应用提取后质量过滤

### 幻灯片之间的结果不一致
- 对 RandomTiler 使用相同的种子
- 使用预处理过滤器标准化染色
- 根据染色质量调整`tissue_percent`
- 实现幻灯片特定蒙版定制

## 资源

此技能包括 `references/` 目录中的详细参考文档：

### 参考文献/slide_management.md
加载、检查和处理整个幻灯片图像的综合指南：
- 幻灯片初始化和配置
- 内置样本数据集
- 幻灯片属性和元数据
- 缩略图生成和可视化
- 使用金字塔级别
- 多幻灯片处理工作流程
- 最佳实践和常见模式

### 参考文献/tissue_masks.md
有关组织检测和掩蔽的完整文档：
- TissueMask、BiggestTissueBoxMask、BinaryMask 类
- 组织检测过滤器的工作原理
- 使用过滤器链定制面罩
- 可视化面具
- 创建自定义矩形和注释排除蒙版
- 与瓦片提取集成
- 最佳实践和故障排除

### 参考文献/tile_extraction.md
瓦片提取策略详解：
- RandomTiler、GridTiler、ScoreTiler 比较
- 可用的记分器（NucleiScorer、CellularityScorer、自定义）
- 通用和特定于策略的参数
- 使用locate_tiles()进行图块预览
- 提取工作流程和 CSV 报告
- 高级模式（多级、分层）
- 性能优化
- 解决常见问题

### 参考文献/filters_preprocessing.md
完整的过滤器参考和预处理指南：
- 图像过滤器（颜色转换、阈值处理、对比度）
- 形态过滤器（膨胀、腐蚀、打开、关闭）
- 过滤器组成和链接
- 常见的预处理管道
- 对瓷砖应用过滤器
- 定制面罩过滤器
- 质量控制过滤器
- 最佳实践和故障排除

### 参考文献/可视化.md
综合可视化指南：
- 幻灯片缩略图显示和保存
- 掩模可视化技术
- 瓷砖位置预览
- 显示提取的图块并创建马赛克
- 质量评估可视化
- 多幻灯片比较
- 滤镜效果可视化
- 导出高分辨率图形和 PDF
- Jupyter 笔记本中的交互式可视化

**使用模式：** 参考文件包含支持本主要技能文档中描述的工作流程的深入信息。根据需要加载特定的参考文件，以获取详细的实施指南、故障排除或高级功能。