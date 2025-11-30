<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：pydicom
描述：用于处理 DICOM（医学数字成像和通信）文件的 Python 库。在读取、写入或修改 DICOM 格式的医学成像数据、从医学图像（CT、MRI、X 射线、超声）中提取像素数据、匿名 DICOM 文件、使用 DICOM 元数据和标签、将 DICOM 图像转换为其他格式、处理压缩的 DICOM 数据或处理医学成像数据集时，请使用此技能。适用于涉及医学图像分析、PACS 系统、放射学工作流程和医疗保健成像应用程序的任务。
---

# 皮迪康姆

## 概述

Pydicom 是一个纯 Python 包，用于处理 DICOM 文件，这是医学成像数据的标准格式。此技能提供有关读取、写入和操作 DICOM 文件的指导，包括使用像素数据、元数据和各种压缩格式。

## 何时使用此技能

在处理以下情况时使用此技能：
- 医学影像文件（CT、MRI、X 射线、超声波、PET 等）
- 需要元数据提取或修改的 DICOM 数据集
- 医学扫描的像素数据提取和图像处理
- 用于研究或数据共享的 DICOM 匿名化
- 将 DICOM 文件转换为标准图像格式
- 需要解压的压缩 DICOM 数据
- DICOM 序列和结构化报告
- 多切片体积重建
- PACS（图片存档和通信系统）集成

## 安装

安装pydicom和常用依赖项：

```bash
uv pip install pydicom
uv pip install pillow  # For image format conversion
uv pip install numpy   # For pixel array manipulation
uv pip install matplotlib  # For visualization
```

为了处理压缩的 DICOM 文件，可能需要额外的包：

<<<代码块_1>>>

## 核心工作流程

### 读取 DICOM 文件

使用 `pydicom.dcmread()` 读取 DICOM 文件：

<<<代码块_2>>>

**要点：**
- `dcmread()` 返回一个 `Dataset` 对象
- 使用属性表示法（例如，`ds.PatientName`）或标签表示法（例如，`ds[0x0010, 0x0010]`）访问数据元素
- 使用 `ds.file_meta` 访问文件元数据，例如传输语法 UID
- 使用 `getattr(ds, 'AttributeName', default_value)` 或 `hasattr(ds, 'AttributeName')` 处理缺失的属性

### 使用像素数据

从 DICOM 文件中提取和操作图像数据：

<<<代码块_3>>>

**处理彩色图像：**

<<<代码块_4>>>

**多帧图像（视频/系列）：**

<<<代码块_5>>>

### 将 DICOM 转换为图像格式

使用提供的 `dicom_to_image.py` 脚本或手动转换：

<<<代码块_6>>>

使用脚本：`python scripts/dicom_to_image.py input.dcm output.png`

### 修改元数据

修改 DICOM 数据元素：

```python
import pydicom
from datetime import datetime

ds = pydicom.dcmread('input.dcm')

# Modify existing elements
ds.PatientName = "Doe^John"
ds.StudyDate = datetime.now().strftime('%Y%m%d')
ds.StudyDescription = "Modified Study"

# Add new elements
ds.SeriesNumber = 1
ds.SeriesDescription = "New Series"

# Remove elements
if hasattr(ds, 'PatientComments'):
    delattr(ds, 'PatientComments')
# Or using del
if 'PatientComments' in ds:
    del ds.PatientComments

# Save modified file
ds.save_as('modified.dcm')
```

### 匿名 DICOM 文件

删除或替换患者身份信息：

```python
import pydicom
from datetime import datetime

ds = pydicom.dcmread('input.dcm')

# Tags commonly containing PHI (Protected Health Information)
tags_to_anonymize = [
    'PatientName', 'PatientID', 'PatientBirthDate',
    'PatientSex', 'PatientAge', 'PatientAddress',
    'InstitutionName', 'InstitutionAddress',
    'ReferringPhysicianName', 'PerformingPhysicianName',
    'OperatorsName', 'StudyDescription', 'SeriesDescription',
]

# Remove or replace sensitive data
for tag in tags_to_anonymize:
    if hasattr(ds, tag):
        if tag in ['PatientName', 'PatientID']:
            setattr(ds, tag, 'ANONYMOUS')
        elif tag == 'PatientBirthDate':
            setattr(ds, tag, '19000101')
        else:
            delattr(ds, tag)

# Update dates to maintain temporal relationships
if hasattr(ds, 'StudyDate'):
    # Shift dates by a random offset
    ds.StudyDate = '20000101'

# Keep pixel data intact
ds.save_as('anonymized.dcm')
```

使用提供的脚本：`python scripts/anonymize_dicom.py input.dcm output.dcm`

### 写入 DICOM 文件

从头开始创建 DICOM 文件：

```python
import pydicom
from pydicom.dataset import Dataset, FileDataset
from datetime import datetime
import numpy as np

# Create file meta information
file_meta = Dataset()
file_meta.MediaStorageSOPClassUID = pydicom.uid.generate_uid()
file_meta.MediaStorageSOPInstanceUID = pydicom.uid.generate_uid()
file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian

# Create the FileDataset instance
ds = FileDataset('new_dicom.dcm', {}, file_meta=file_meta, preamble=b"\0" * 128)

# Add required DICOM elements
ds.PatientName = "Test^Patient"
ds.PatientID = "123456"
ds.Modality = "CT"
ds.StudyDate = datetime.now().strftime('%Y%m%d')
ds.StudyTime = datetime.now().strftime('%H%M%S')
ds.ContentDate = ds.StudyDate
ds.ContentTime = ds.StudyTime

# Add image-specific elements
ds.SamplesPerPixel = 1
ds.PhotometricInterpretation = "MONOCHROME2"
ds.Rows = 512
ds.Columns = 512
ds.BitsAllocated = 16
ds.BitsStored = 16
ds.HighBit = 15
ds.PixelRepresentation = 0

# Create pixel data
pixel_array = np.random.randint(0, 4096, (512, 512), dtype=np.uint16)
ds.PixelData = pixel_array.tobytes()

# Add required UIDs
ds.SOPClassUID = pydicom.uid.CTImageStorage
ds.SOPInstanceUID = file_meta.MediaStorageSOPInstanceUID
ds.SeriesInstanceUID = pydicom.uid.generate_uid()
ds.StudyInstanceUID = pydicom.uid.generate_uid()

# Save the file
ds.save_as('new_dicom.dcm')
```

### 压缩和解压

处理压缩的 DICOM 文件：

```python
import pydicom

# Read compressed DICOM file
ds = pydicom.dcmread('compressed.dcm')

# Check transfer syntax
print(f"Transfer Syntax: {ds.file_meta.TransferSyntaxUID}")
print(f"Transfer Syntax Name: {ds.file_meta.TransferSyntaxUID.name}")

# Decompress and save as uncompressed
ds.decompress()
ds.save_as('uncompressed.dcm', write_like_original=False)

# Or compress when saving (requires appropriate encoder)
ds_uncompressed = pydicom.dcmread('uncompressed.dcm')
ds_uncompressed.compress(pydicom.uid.JPEGBaseline8Bit)
ds_uncompressed.save_as('compressed_jpeg.dcm')
```

**常用传输语法：**
- `ExplicitVRLittleEndian` - 未压缩，最常见
- `JPEGBaseline8Bit` - JPEG 有损压缩
- `JPEGLossless` - JPEG 无损压缩
- `JPEG2000Lossless` - JPEG 2000 无损
- `RLELossless` - 无损游程编码

请参阅 `references/transfer_syntaxes.md` 了解完整列表。

### 使用 DICOM 序列

处理嵌套数据结构：

```python
import pydicom

ds = pydicom.dcmread('file.dcm')

# Access sequences
if 'ReferencedStudySequence' in ds:
    for item in ds.ReferencedStudySequence:
        print(f"Referenced SOP Instance UID: {item.ReferencedSOPInstanceUID}")

# Create a sequence
from pydicom.sequence import Sequence

sequence_item = Dataset()
sequence_item.ReferencedSOPClassUID = pydicom.uid.CTImageStorage
sequence_item.ReferencedSOPInstanceUID = pydicom.uid.generate_uid()

ds.ReferencedImageSequence = Sequence([sequence_item])
```

### 处理 DICOM 系列

处理多个相关的 DICOM 文件：

```python
import pydicom
import numpy as np
from pathlib import Path

# Read all DICOM files in a directory
dicom_dir = Path('dicom_series/')
slices = []

for file_path in dicom_dir.glob('*.dcm'):
    ds = pydicom.dcmread(file_path)
    slices.append(ds)

# Sort by slice location or instance number
slices.sort(key=lambda x: float(x.ImagePositionPatient[2]))
# Or: slices.sort(key=lambda x: int(x.InstanceNumber))

# Create 3D volume
volume = np.stack([s.pixel_array for s in slices])
print(f"Volume shape: {volume.shape}")  # (num_slices, rows, columns)

# Get spacing information for proper scaling
pixel_spacing = slices[0].PixelSpacing  # [row_spacing, col_spacing]
slice_thickness = slices[0].SliceThickness
print(f"Voxel size: {pixel_spacing[0]}x{pixel_spacing[1]}x{slice_thickness} mm")
```

## 帮助脚本

此技能包括 `scripts/` 目录中的实用程序脚本：

### anonymize_dicom.py
通过删除或替换受保护的健康信息 (PHI) 对 DICOM 文件进行匿名化。

```bash
python scripts/anonymize_dicom.py input.dcm output.dcm
```

### dicom_to_image.py
将 DICOM 文件转换为常见图像格式（PNG、JPEG、TIFF）。

```bash
python scripts/dicom_to_image.py input.dcm output.png
python scripts/dicom_to_image.py input.dcm output.jpg --format JPEG
```

### extract_metadata.py
以可读格式提取并显示 DICOM 元数据。

```bash
python scripts/extract_metadata.py file.dcm
python scripts/extract_metadata.py file.dcm --output metadata.txt
```

## 参考资料

详细参考信息可在 `references/` 目录中找到：

- **common_tags.md**：按类别（患者、研究、系列、图像等）组织的常用 DICOM 标签的综合列表
- **transfer_syntaxes.md**：DICOM 传输语法和压缩格式的完整参考

## 常见问题及解决方案

**问题：“无法解码像素数据”**
- 解决方案：安装额外的压缩处理程序：`uv pip install pylibjpeg pylibjpeg-libjpeg python-gdcm`

**问题：访问标签时出现“AttributeError”**
- 解决方案：检查属性是否存在 `hasattr(ds, 'AttributeName')` 或使用 `ds.get('AttributeName', default)`

**问题：图像显示不正确（太暗/太亮）**
- 解决方案：应用 VOI LUT 窗口化：`apply_voi_lut(pixel_array, ds)` 或使用 `WindowCenter` 和 `WindowWidth` 手动调整

**问题：大型系列的内存问题**
- 解决方案：迭代处理文件、使用内存映射数组或对图像进行下采样

## 最佳实践

1. **在使用 `hasattr()` 或 `get()` 访问所需属性之前，始终检查所需属性**
2. **通过使用 `save_as()` 和 `write_like_original=True` 修改文件时保留文件元数据**
3. **在处理像素数据之前使用传输语法 UID** 了解压缩格式
4. **从不受信任的来源读取文件时处理异常**
5. **应用适当的加窗** (VOI LUT) 进行医学图像可视化
6. **处理 3D 体积时维护空间信息**（像素间距、切片厚度）
7. **在共享医疗数据之前彻底验证匿名性**
8. **正确使用UID** - 创建新实例时生成新的UID，修改时保留它们

## 文档

pydicom 官方文档：https://pydicom.github.io/pydicom/dev/
- 用户指南：https://pydicom.github.io/pydicom/dev/guides/user/index.html
- 教程：https://pydicom.github.io/pydicom/dev/tutorials/index.html
- API 参考：https://pydicom.github.io/pydicom/dev/reference/index.html
- 示例：https://pydicom.github.io/pydicom/dev/auto_examples/index.html