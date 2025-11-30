<!-- 此文件由机器翻译自 tissue_masks.md -->

# 面膜纸

## 概述

组织掩模是识别整个幻灯片图像内的组织区域的二进制表示。它们对于在图块提取过程中过滤掉背景、伪影和非组织区域至关重要。 Histolab 提供多种掩模类别来满足不同的组织分割需求。

## 面具课程

### 二进制掩码

**用途：** 用于创建自定义二进制掩码的通用基类。

```python
from histolab.masks import BinaryMask

class CustomMask(BinaryMask):
    def _mask(self, obj):
        # Implement custom masking logic
        # Return binary numpy array
        pass
```

**使用案例：**
- 自定义组织分割算法
- 特定区域分析（例如，排除注释）
- 与外部细分模型集成

### 组织面膜

**目的：** 使用自动过滤器分割载玻片中的所有组织区域。

<<<代码块_1>>>

**它是如何工作的：**
1. 将图像转换为灰度
2. 应用 Otsu 阈值将组织与背景分离
3. 执行二元扩张以连接附近的组织区域
4. 去除组织区域内的小孔
5. 过滤掉小物体（伪影）

**返回：** 二进制 NumPy 数组，其中：
- `True`（或 1）：组织像素
- `False`（或 0）：背景像素

**最适合：**
- 具有多个独立组织切片的载玻片
- 全面的组织分析
- 当所有组织区域都很重要时

### BiggestTissueBoxMask（默认）

**用途：** 识别并返回最大连接组织区域的边界框。

<<<代码块_2>>>

**它是如何工作的：**
1.应用与TissueMask相同的过滤管道
2. 识别所有连接的组织成分
3. 选择最大的连通分量
4. 返回包围该区域的边界框

**最适合：**
- 具有单个主要组织切片的载玻片
- 排除小伪影或组织碎片
- 专注于主要组织区域（大多数瓷砖工的默认设置）

## 自定义带有滤镜的蒙版

面罩接受定制过滤器链以进行专门的组织检测：

<<<代码块_3>>>

## 可视化掩模

### 使用locate_mask()

<<<代码块_4>>>

这将显示幻灯片缩略图，其中蒙版边界以对比色覆盖。

### 手动可视化

<<<代码块_5>>>

## 创建自定义矩形遮罩

定义特定的感兴趣区域：

<<<代码块_6>>>

## 排除注释

病理学幻灯片通常包含笔标记或数字注释。使用自定义掩码排除它们：

```python
from histolab.masks import TissueMask
from histolab.filters.image_filters import RgbToGrayscale, OtsuThreshold
from histolab.filters.morphological_filters import BinaryDilation

class AnnotationExclusionMask(BinaryMask):
    def _mask(self, obj):
        thumb = obj.thumbnail

        # Convert to HSV to detect pen marks (often blue/green)
        hsv = cv2.cvtColor(np.array(thumb), cv2.COLOR_RGB2HSV)

        # Define color ranges for pen marks
        lower_blue = np.array([100, 50, 50])
        upper_blue = np.array([130, 255, 255])

        # Create mask excluding pen marks
        pen_mask = cv2.inRange(hsv, lower_blue, upper_blue)

        # Apply standard tissue detection
        tissue_mask = TissueMask()(obj)

        # Combine: keep tissue, exclude pen marks
        final_mask = tissue_mask & ~pen_mask.astype(bool)

        return final_mask
```

## 与 Tile Extraction 集成

遮罩通过 `extraction_mask` 参数与平铺器无缝集成：

```python
from histolab.tiler import RandomTiler
from histolab.masks import TissueMask, BiggestTissueBoxMask

# Use TissueMask to extract from all tissue
random_tiler = RandomTiler(
    tile_size=(512, 512),
    n_tiles=100,
    level=0,
    extraction_mask=TissueMask()  # Extract from all tissue regions
)

# Or use default BiggestTissueBoxMask
random_tiler = RandomTiler(
    tile_size=(512, 512),
    n_tiles=100,
    level=0,
    extraction_mask=BiggestTissueBoxMask()  # Default behavior
)
```

## 最佳实践

1. **提取前预览蒙版**：使用`locate_mask()`或手动可视化来验证蒙版质量
2. **选择适当的掩模类型**：对于多个组织切片使用 `TissueMask`，对于单个主要切片使用 `BiggestTissueBoxMask`
3. **针对特定染色进行定制**：不同染色（H&E、IHC）可能需要调整阈值参数
4. **处理伪影**：使用自定义滤镜或遮罩排除笔痕、气泡或折叠
5. **在不同的载玻片上进行测试**：验证具有不同质量和伪影的载玻片的掩模性能
6. **考虑计算成本**：`TissueMask` 比 `BiggestTissueBoxMask` 更全面，但计算量更大

## 常见问题及解决方案

### 问题：蒙版包含太多背景
**解决方案：** 调整Otsu阈值或增加小物体去除阈值

### 问题：面罩不包括有效组织
**解决方案：** 降低小物体去除阈值或修改膨胀参数

### 问题：多个组织切片，但仅捕获最大的组织切片
**解决方案：** 从 `BiggestTissueBoxMask` 切换到 `TissueMask`

### 问题：遮罩中包含笔注释
**解决方案：** 实现自定义注释排除掩码（参见上面的示例）