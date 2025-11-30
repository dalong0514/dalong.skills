<!-- 此文件由机器翻译自 transfer_syntaxes.md -->

# DICOM 传输语法参考

本文档提供了 DICOM 传输语法和压缩格式的综合参考。传输语法定义 DICOM 数据的编码方式，包括字节排序、压缩方法和其他编码规则。

## 概述

传输语法 UID 指定：
1. **字节顺序**：Little Endian 或 Big Endian
2. **值表示（VR）**：隐式或显式
3. **压缩**：无，或特定压缩算法

## 未压缩传输语法

### 隐式 VR Little Endian (1.2.840.10008.1.2)
- **默认**传输语法
- 值表示是隐式的（未显式编码）
- 小端字节顺序
- **Pydicom 常数**：`pydicom.uid.ImplicitVRLittleEndian`

**用途：**
```python
import pydicom
ds.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
```

### 显式 VR Little Endian (1.2.840.10008.1.2.1)
- **最常见**传输语法
- 值表示是明确的
- 小端字节顺序
- **Pydicom 常数**：`pydicom.uid.ExplicitVRLittleEndian`

**用途：**
<<<代码块_1>>>

### 显式 VR Big Endian (1.2.840.10008.1.2.2) - 已退休
- 值表示是明确的
- 大端字节顺序
- **已弃用** - 不建议用于新实现
- **Pydicom 常数**：`pydicom.uid.ExplicitVRBigEndian`

## JPEG 压缩

### JPEG 基线（过程 1）(1.2.840.10008.1.2.4.50)
- **有损**压缩
- 仅 8 位样本
- 最广泛支持的 JPEG 格式
- **Pydicom 常数**：`pydicom.uid.JPEGBaseline8Bit`

**依赖项：** 需要 `pylibjpeg` 或 `pillow`

**用途：**
<<<代码块_2>>>

### JPEG 扩展（过程 2 和 4）(1.2.840.10008.1.2.4.51)
- **有损**压缩
- 8 位和 12 位样本
- **Pydicom 常数**：`pydicom.uid.JPEGExtended12Bit`

### JPEG 无损、非分层（过程 14）(1.2.840.10008.1.2.4.57)
- **无损**压缩
- 一阶预测
- **Pydicom 常数**：`pydicom.uid.JPEGLossless`

**依赖项：** 需要 `pylibjpeg-libjpeg` 或 `gdcm`

### JPEG 无损、非分层、一阶预测 (1.2.840.10008.1.2.4.70)
- **无损**压缩
- 使用过程 14 选择值 1
- **Pydicom 常数**：`pydicom.uid.JPEGLosslessSV1`

**用途：**
<<<代码块_3>>>

### JPEG-LS 无损 (1.2.840.10008.1.2.4.80)
- **无损**压缩
- 复杂度低，压缩性好
- **Pydicom 常数**：`pydicom.uid.JPEGLSLossless`

**依赖项：** 需要 `pylibjpeg-libjpeg` 或 `gdcm`

### JPEG-LS 有损（接近无损）(1.2.840.10008.1.2.4.81)
- **近乎无损**压缩
- 允许控制精度损失
- **Pydicom 常数**：`pydicom.uid.JPEGLSNearLossless`

## JPEG 2000 压缩

### 仅 JPEG 2000 无损 (1.2.840.10008.1.2.4.90)
- **无损**压缩
- 基于小波的压缩
- 比 JPEG 无损压缩更好
- **Pydicom 常数**：`pydicom.uid.JPEG2000Lossless`

**依赖项：** 需要 `pylibjpeg-openjpeg`、`gdcm` 或 `pillow`

**用途：**
<<<代码块_4>>>

### JPEG 2000 (1.2.840.10008.1.2.4.91)
- **有损或无损**压缩
- 基于小波的压缩
- 低比特率下的高质量
- **Pydicom 常数**：`pydicom.uid.JPEG2000`

**依赖项：** 需要 `pylibjpeg-openjpeg`、`gdcm` 或 `pillow`

### JPEG 2000 第 2 部分多分量无损 (1.2.840.10008.1.2.4.92)
- **无损**压缩
- 支持多组件图像
- **Pydicom 常数**：`pydicom.uid.JPEG2000MCLossless`

### JPEG 2000 第 2 部分多分量 (1.2.840.10008.1.2.4.93)
- **有损或无损**压缩
- 支持多组件图像
- **Pydicom 常数**：`pydicom.uid.JPEG2000MC`

## RLE 压缩

### RLE 无损 (1.2.840.10008.1.2.5)
- **无损**压缩
- 行程编码
- 简单、快速的算法
- 适用于具有重复值的图像
- **Pydicom 常数**：`pydicom.uid.RLELossless`

**依赖项：** 内置于 pydicom（不需要额外的包）

**用途：**
<<<代码块_5>>>

## 压缩传输语法

### 压缩显式 VR Little Endian (1.2.840.10008.1.2.1.99)
- 对整个数据集使用 ZLIB 压缩
- 不常用
- **Pydicom 常数**：`pydicom.uid.DeflatedExplicitVRLittleEndian`

## MPEG 压缩

### MPEG2 主要配置文件 @ 主级别 (1.2.840.10008.1.2.4.100)
- **有损**视频压缩
- 对于多帧图像/视频
- **Pydicom 常数**：`pydicom.uid.MPEG2MPML`

### MPEG2 主要配置文件 @ 高级别 (1.2.840.10008.1.2.4.101)
- **有损**视频压缩
- 比 MPML 更高的分辨率
- **Pydicom 常数**：`pydicom.uid.MPEG2MPHL`

### MPEG-4 AVC/H.264 高配置 (1.2.840.10008.1.2.4.102-106)
- **有损**视频压缩
- 各种级别（BD、2D、3D、立体）
- 现代视频编解码器

## 检查传输语法

### 识别当前传输语法
<<<代码块_6>>>

### 常见检查
```python
# Check if little endian
if ts_uid.is_little_endian:
    print("Little Endian")

# Check if implicit VR
if ts_uid.is_implicit_VR:
    print("Implicit VR")

# Check compression type
if 'JPEG' in ts_uid.name:
    print("JPEG compressed")
elif 'JPEG2000' in ts_uid.name:
    print("JPEG 2000 compressed")
elif 'RLE' in ts_uid.name:
    print("RLE compressed")
```

## 解压

### 自动减压
Pydicom在访问`pixel_array`时可以自动解压缩像素数据：

```python
import pydicom

# Read compressed DICOM
ds = pydicom.dcmread('compressed.dcm')

# Pixel data is automatically decompressed
pixel_array = ds.pixel_array  # Decompresses if needed
```

### 手动解压
```python
import pydicom

ds = pydicom.dcmread('compressed.dcm')

# Decompress in-place
ds.decompress()

# Now save as uncompressed
ds.save_as('uncompressed.dcm', write_like_original=False)
```

## 压缩

### 压缩 DICOM 文件
```python
import pydicom

ds = pydicom.dcmread('uncompressed.dcm')

# Compress using JPEG 2000 Lossless
ds.compress(pydicom.uid.JPEG2000Lossless)
ds.save_as('compressed_j2k.dcm')

# Compress using RLE Lossless (no additional dependencies)
ds.compress(pydicom.uid.RLELossless)
ds.save_as('compressed_rle.dcm')

# Compress using JPEG Baseline (lossy)
ds.compress(pydicom.uid.JPEGBaseline8Bit)
ds.save_as('compressed_jpeg.dcm')
```

### 使用自定义编码参数进行压缩
```python
import pydicom
from pydicom.encoders import JPEGLSLosslessEncoder

ds = pydicom.dcmread('uncompressed.dcm')

# Compress with custom parameters
ds.compress(pydicom.uid.JPEGLSLossless, encoding_plugin='pylibjpeg')
```

## 安装压缩处理程序

不同的传输语法需要不同的Python包：

### JPEG 基线/扩展
```bash
pip install pylibjpeg pylibjpeg-libjpeg
# Or
pip install pillow
```

### JPEG 无损/JPEG-LS
```bash
pip install pylibjpeg pylibjpeg-libjpeg
# Or
pip install python-gdcm
```

### JPEG 2000
```bash
pip install pylibjpeg pylibjpeg-openjpeg
# Or
pip install python-gdcm
# Or
pip install pillow
```

### RLE
无需额外的软件包 - 内置于 pydicom 中

### 综合安装
```bash
# Install all common handlers
pip install pylibjpeg pylibjpeg-libjpeg pylibjpeg-openjpeg python-gdcm
```

## 检查可用的处理程序

```python
import pydicom

# List available pixel data handlers
from pydicom.pixel_data_handlers.util import get_pixel_data_handlers
handlers = get_pixel_data_handlers()

print("Available handlers:")
for handler in handlers:
    print(f"  - {handler.__name__}")
```

## 最佳实践

1. **在创建新文件时使用显式 VR Little Endian** 以获得最大兼容性
2. **使用 JPEG 2000 Lossless** 获得良好的压缩效果，且不会造成质量损失
3. **如果无法安装额外的依赖项，请使用 RLE Lossless**
4. **在处理之前检查传输语法**，以确保您拥有正确的处理程序
5. **部署前测试解压**，确保安装所有必需的包
6. **尽可能使用 `write_like_original=True` 保留原始**传输语法
7. **选择有损压缩时考虑文件大小**与质量权衡
8. **对诊断图像使用无损压缩**以保持临床质量

## 常见问题

### 问题：“无法解码像素数据”
**原因：** 缺少压缩处理程序
**解决方案：** 安装适当的包（请参阅上面的安装压缩处理程序）

### 问题：“不支持的传输语法”
**原因：** 罕见或不受支持的压缩格式
**解决方案：** 尝试安装支持更多格式的`python-gdcm`

### 问题：“像素数据已解压但看起来错误”
**原因：** 可能需要应用 VOI LUT 或重新缩放
**解决方案：** 使用 `apply_voi_lut()` 或应用 `RescaleSlope`/`RescaleIntercept`

## 参考文献

- DICOM 标准第 5 部分（数据结构和编码）：https://dicom.nema.org/medical/dicom/current/output/chtml/part05/PS3.5.html
- Pydicom 传输语法文档：https://pydicom.github.io/pydicom/stable/guides/user/transfer_syntaxes.html
- Pydicom 压缩指南：https://pydicom.github.io/pydicom/stable/old/image_data_compression.html