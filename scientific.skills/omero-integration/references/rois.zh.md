<!-- 此文件由机器翻译自 rois.md -->

# 感兴趣区域 (ROI)

本参考内容涵盖了在 OMERO 中创建、检索和分析 ROI。

## 投资回报率概述

OMERO 中的 ROI（感兴趣区域）是用于标记图像上特定区域的几何形状的容器。每个 ROI 可以包含多个形状，并且形状可以特定于 Z 截面和时间点。

### 支持的形状类型

- **矩形**：矩形区域
- **椭圆**：圆形和椭圆形区域
- **线**：线段
- **点**：单点
- **多边形**：多点多边形
- **掩模**：基于像素的掩模
- **多段线**：多段线

## 创造投资回报率

### 辅助函数

```python
from omero.rtypes import rdouble, rint, rstring
import omero.model

def create_roi(conn, image, shapes):
    """
    Create an ROI and link it to shapes.

    Args:
        conn: BlitzGateway connection
        image: Image object
        shapes: List of shape objects

    Returns:
        Saved ROI object
    """
    roi = omero.model.RoiI()
    roi.setImage(image._obj)

    for shape in shapes:
        roi.addShape(shape)

    updateService = conn.getUpdateService()
    return updateService.saveAndReturnObject(roi)

def rgba_to_int(red, green, blue, alpha=255):
    """
    Convert RGBA values (0-255) to integer encoding for OMERO.

    Args:
        red, green, blue, alpha: Color values (0-255)

    Returns:
        Integer color value
    """
    return int.from_bytes([red, green, blue, alpha],
                          byteorder='big', signed=True)
```

### 矩形 ROI

<<<代码块_1>>>

### 椭圆 ROI

<<<代码块_2>>>

### 线路投资回报率

<<<代码块_3>>>

### 投资回报率点

<<<代码块_4>>>

### 多边形投资回报率

<<<代码块_5>>>

### 掩模投资回报率

<<<代码块_6>>>

## 一个 ROI 中的多个形状

```python
# Create multiple shapes for the same ROI
shapes = []

# Rectangle
rect = omero.model.RectangleI()
rect.x = rdouble(100)
rect.y = rdouble(100)
rect.width = rdouble(50)
rect.height = rdouble(50)
rect.theZ = rint(0)
rect.theT = rint(0)
shapes.append(rect)

# Ellipse
ellipse = omero.model.EllipseI()
ellipse.x = rdouble(125)
ellipse.y = rdouble(125)
ellipse.radiusX = rdouble(20)
ellipse.radiusY = rdouble(20)
ellipse.theZ = rint(0)
ellipse.theT = rint(0)
shapes.append(ellipse)

# Create single ROI with both shapes
roi = create_roi(conn, image, shapes)
```

## 检索 ROI

### 获取图像的所有 ROI

```python
# Get ROI service
roi_service = conn.getRoiService()

# Find all ROIs for image
result = roi_service.findByImage(image_id, None)

print(f"Found {len(result.rois)} ROIs")

for roi in result.rois:
    print(f"ROI ID: {roi.getId().getValue()}")
    print(f"  Number of shapes: {len(roi.copyShapes())}")
```

### 解析 ROI 形状

```python
import omero.model

result = roi_service.findByImage(image_id, None)

for roi in result.rois:
    roi_id = roi.getId().getValue()
    print(f"ROI ID: {roi_id}")

    for shape in roi.copyShapes():
        shape_id = shape.getId().getValue()
        z = shape.getTheZ().getValue() if shape.getTheZ() else None
        t = shape.getTheT().getValue() if shape.getTheT() else None

        # Get label
        label = ""
        if shape.getTextValue():
            label = shape.getTextValue().getValue()

        print(f"  Shape ID: {shape_id}, Z: {z}, T: {t}, Label: {label}")

        # Type-specific parsing
        if isinstance(shape, omero.model.RectangleI):
            x = shape.getX().getValue()
            y = shape.getY().getValue()
            width = shape.getWidth().getValue()
            height = shape.getHeight().getValue()
            print(f"    Rectangle: ({x}, {y}) {width}x{height}")

        elif isinstance(shape, omero.model.EllipseI):
            x = shape.getX().getValue()
            y = shape.getY().getValue()
            rx = shape.getRadiusX().getValue()
            ry = shape.getRadiusY().getValue()
            print(f"    Ellipse: center ({x}, {y}), radii ({rx}, {ry})")

        elif isinstance(shape, omero.model.PointI):
            x = shape.getX().getValue()
            y = shape.getY().getValue()
            print(f"    Point: ({x}, {y})")

        elif isinstance(shape, omero.model.LineI):
            x1 = shape.getX1().getValue()
            y1 = shape.getY1().getValue()
            x2 = shape.getX2().getValue()
            y2 = shape.getY2().getValue()
            print(f"    Line: ({x1}, {y1}) to ({x2}, {y2})")

        elif isinstance(shape, omero.model.PolygonI):
            points = shape.getPoints().getValue()
            print(f"    Polygon: {points}")

        elif isinstance(shape, omero.model.MaskI):
            x = shape.getX().getValue()
            y = shape.getY().getValue()
            width = shape.getWidth().getValue()
            height = shape.getHeight().getValue()
            print(f"    Mask: ({x}, {y}) {width}x{height}")
```

## 分析 ROI 强度

### 获取 ROI 形状的统计数据

```python
# Get all shapes from ROIs
roi_service = conn.getRoiService()
result = roi_service.findByImage(image_id, None)

shape_ids = []
for roi in result.rois:
    for shape in roi.copyShapes():
        shape_ids.append(shape.id.val)

# Define position
z, t = 0, 0
channel_index = 0

# Get statistics
stats = roi_service.getShapeStatsRestricted(
    shape_ids, z, t, [channel_index]
)

# Display statistics
for i, stat in enumerate(stats):
    shape_id = shape_ids[i]
    print(f"Shape {shape_id} statistics:")
    print(f"  Points Count: {stat.pointsCount[channel_index]}")
    print(f"  Min: {stat.min[channel_index]}")
    print(f"  Mean: {stat.mean[channel_index]}")
    print(f"  Max: {stat.max[channel_index]}")
    print(f"  Sum: {stat.sum[channel_index]}")
    print(f"  Std Dev: {stat.stdDev[channel_index]}")
```

### 提取 ROI 内的像素值

```python
import numpy as np

# Get image and ROI
image = conn.getObject("Image", image_id)
result = roi_service.findByImage(image_id, None)

# Get first rectangle shape
roi = result.rois[0]
rect = roi.copyShapes()[0]

# Get rectangle bounds
x = int(rect.getX().getValue())
y = int(rect.getY().getValue())
width = int(rect.getWidth().getValue())
height = int(rect.getHeight().getValue())
z = rect.getTheZ().getValue()
t = rect.getTheT().getValue()

# Get pixel data
pixels = image.getPrimaryPixels()

# Extract region for each channel
for c in range(image.getSizeC()):
    # Get plane
    plane = pixels.getPlane(z, c, t)

    # Extract ROI region
    roi_region = plane[y:y+height, x:x+width]

    print(f"Channel {c}:")
    print(f"  Mean intensity: {np.mean(roi_region)}")
    print(f"  Max intensity: {np.max(roi_region)}")
```

## 修改 ROI

### 更新形状属性

```python
# Get ROI and shape
result = roi_service.findByImage(image_id, None)
roi = result.rois[0]
shape = roi.copyShapes()[0]

# Modify shape (example: change rectangle size)
if isinstance(shape, omero.model.RectangleI):
    shape.setWidth(rdouble(150))
    shape.setHeight(rdouble(100))
    shape.setTextValue(rstring("Updated Rectangle"))

# Save changes
updateService = conn.getUpdateService()
updated_roi = updateService.saveAndReturnObject(roi._obj)
```

### 从 ROI 中删除形状

```python
result = roi_service.findByImage(image_id, None)

for roi in result.rois:
    for shape in roi.copyShapes():
        # Check condition (e.g., remove by label)
        if (shape.getTextValue() and
            shape.getTextValue().getValue() == "test-Ellipse"):

            print(f"Removing shape {shape.getId().getValue()}")
            roi.removeShape(shape)

            # Save modified ROI
            updateService = conn.getUpdateService()
            roi = updateService.saveAndReturnObject(roi)
```

## 删除 ROI

### 删除单个 ROI

```python
# Delete ROI by ID
roi_id = 123
conn.deleteObjects("Roi", [roi_id], wait=True)
print(f"Deleted ROI {roi_id}")
```

### 删除图像的所有 ROI

```python
# Get all ROI IDs for image
result = roi_service.findByImage(image_id, None)
roi_ids = [roi.getId().getValue() for roi in result.rois]

# Delete all
if roi_ids:
    conn.deleteObjects("Roi", roi_ids, wait=True)
    print(f"Deleted {len(roi_ids)} ROIs")
```

## 批量创建 ROI

### 为多个图像创建 ROI

```python
# Get images
dataset = conn.getObject("Dataset", dataset_id)

for image in dataset.listChildren():
    # Create rectangle at center of each image
    x = image.getSizeX() // 2 - 50
    y = image.getSizeY() // 2 - 50

    rect = omero.model.RectangleI()
    rect.x = rdouble(x)
    rect.y = rdouble(y)
    rect.width = rdouble(100)
    rect.height = rdouble(100)
    rect.theZ = rint(0)
    rect.theT = rint(0)
    rect.textValue = rstring("Auto ROI")

    roi = create_roi(conn, image, [rect])
    print(f"Created ROI for image {image.getName()}")
```

### 在 Z-Stack 中创建 ROI

```python
image = conn.getObject("Image", image_id)
size_z = image.getSizeZ()

# Create rectangle on each Z-section
shapes = []
for z in range(size_z):
    rect = omero.model.RectangleI()
    rect.x = rdouble(100)
    rect.y = rdouble(100)
    rect.width = rdouble(50)
    rect.height = rdouble(50)
    rect.theZ = rint(z)
    rect.theT = rint(0)
    shapes.append(rect)

# Single ROI with shapes across Z
roi = create_roi(conn, image, shapes)
```

## 完整示例

```python
from omero.gateway import BlitzGateway
from omero.rtypes import rdouble, rint, rstring
import omero.model

HOST = 'omero.example.com'
PORT = 4064
USERNAME = 'user'
PASSWORD = 'pass'

def rgba_to_int(r, g, b, a=255):
    return int.from_bytes([r, g, b, a], byteorder='big', signed=True)

with BlitzGateway(USERNAME, PASSWORD, host=HOST, port=PORT) as conn:
    # Get image
    image = conn.getObject("Image", image_id)
    print(f"Processing: {image.getName()}")

    # Create multiple ROIs
    updateService = conn.getUpdateService()

    # ROI 1: Rectangle
    roi1 = omero.model.RoiI()
    roi1.setImage(image._obj)

    rect = omero.model.RectangleI()
    rect.x = rdouble(50)
    rect.y = rdouble(50)
    rect.width = rdouble(100)
    rect.height = rdouble(100)
    rect.theZ = rint(0)
    rect.theT = rint(0)
    rect.textValue = rstring("Cell 1")
    rect.strokeColor = rint(rgba_to_int(255, 0, 0, 255))

    roi1.addShape(rect)
    roi1 = updateService.saveAndReturnObject(roi1)
    print(f"Created ROI 1: {roi1.getId().getValue()}")

    # ROI 2: Ellipse
    roi2 = omero.model.RoiI()
    roi2.setImage(image._obj)

    ellipse = omero.model.EllipseI()
    ellipse.x = rdouble(200)
    ellipse.y = rdouble(150)
    ellipse.radiusX = rdouble(40)
    ellipse.radiusY = rdouble(30)
    ellipse.theZ = rint(0)
    ellipse.theT = rint(0)
    ellipse.textValue = rstring("Cell 2")
    ellipse.strokeColor = rint(rgba_to_int(0, 255, 0, 255))

    roi2.addShape(ellipse)
    roi2 = updateService.saveAndReturnObject(roi2)
    print(f"Created ROI 2: {roi2.getId().getValue()}")

    # Retrieve and analyze
    roi_service = conn.getRoiService()
    result = roi_service.findByImage(image_id, None)

    shape_ids = []
    for roi in result.rois:
        for shape in roi.copyShapes():
            shape_ids.append(shape.id.val)

    # Get statistics
    stats = roi_service.getShapeStatsRestricted(shape_ids, 0, 0, [0])

    for i, stat in enumerate(stats):
        print(f"Shape {shape_ids[i]}:")
        print(f"  Mean intensity: {stat.mean[0]:.2f}")
```

## 最佳实践

1. **组织形状**：将相关形状分组到单个 ROI 中
2. **标签形状**：使用textValue进行标识
3. **设置 Z 和 T**：始终指定 Z 部分和时间点
4. **颜色编码**：对形状类型使用一致的颜色
5. **验证坐标**：确保形状位于图像边界内
6. **批量创建**：尽可能在单个事务中创建多个 ROI
7. **删除未使用**：删除临时或测试 ROI
8. **导出数据**：将 ROI 统计数据存储在表格中以供以后分析
9. **版本控制**：在注释中记录 ROI 创建方法
10. **性能**：使用形状统计服务代替手动像素提取