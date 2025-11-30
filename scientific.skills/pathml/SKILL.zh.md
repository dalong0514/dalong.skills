<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：pathml
描述：用于分析全切片图像 (WSI) 和多参数成像数据的计算病理学工具包。在处理组织病理学载玻片、H&E 染色图像、多重免疫荧光（CODEX、Vectra）、空间蛋白质组学、细胞核检测/分割、组织图构建或根据病理数据训练 ML 模型时，请使用此技能。支持 160 多种幻灯片格式，包括用于数字病理工作流程的 Aperio SVS、NDPI、DICOM、OME-TIFF。
---

# 路径ML

## 概述

PathML 是一个用于计算病理学工作流程的综合 Python 工具包，旨在促进全幻灯片病理图像的机器学习和图像分析。该框架提供了模块化、可组合的工具，用于加载不同的幻灯片格式、预处理图像、构建空间图、训练深度学习模型以及分析来自 CODEX 和多重免疫荧光等技术的多参数成像数据。

## 何时使用此技能

将此技能应用于：
- 加载和处理各种专有格式的全幻灯片图像（WSI）
- 使用染色归一化预处理 H&E 染色组织图像
- 细胞核检测、分割和分类工作流程
- 构建细胞和组织图以进行空间分析
- 在病理数据上训练或部署机器学习模型（HoVer-Net、HACTNet）
- 分析空间蛋白质组学的多参数成像（CODEX、Vectra、MERFISH）
- 通过多重免疫荧光定量标记物表达
- 使用 HDF5 存储管理大规模病理数据集
- 基于图块的分析和拼接操作

## 核心能力

PathML 提供了参考文件中详细记录的六个主要功能领域：

### 1. 图像加载和格式

从 160 多种专有格式加载整个幻灯片图像，包括 Aperio SVS、Hamamatsu NDPI、Leica SCN、Zeiss ZVI、DICOM 和 OME-TIFF。 PathML 自动处理供应商特定的格式，并提供用于访问图像金字塔、元数据和感兴趣区域的统一接口。

**请参阅：** `references/image_loading.md` 了解支持的格式、加载策略以及使用不同的幻灯片类型。

### 2. 预处理管道

通过组合图像处理、质量控制、染色归一化、组织检测和掩模操作的转换来构建模块化预处理管道。 PathML 的 Pipeline 架构支持跨大型数据集进行可重复、可扩展的预处理。

**关键转变：**
- `StainNormalizationHE` - Macenko/Vahadane 染色归一化
- `TissueDetectionHE`、`NucleusDetectionHE` - 组织/细胞核分割
- `MedianBlur`、`GaussianBlur` - 降噪
- `LabelArtifactTileHE` - 工件的质量控制

**请参阅：** `references/preprocessing.md` 了解完整的转换目录、管道构建和预处理工作流程。

### 3.图构建

构建代表细胞和组织水平关系的空间图。从分段对象中提取特征以创建适合图神经网络和空间分析的基于图的表示。

**请参阅：** `references/graphs.md` 了解图形构建方法、特征提取和空间分析工作流程。

### 4.机器学习

训练和部署用于核检测、分割和分类的深度学习模型。 PathML 将 PyTorch 与预构建模型（HoVer-Net、HACTNet）、自定义 DataLoaders 和 ONNX 支持集成以进行推理。

**主要型号：**
- **HoVer-Net** - 同时进行细胞核分割和分类
- **HACTNet** - 分层细胞类型分类

**请参阅：** `references/machine_learning.md` 了解模型训练、评估、推理工作流程以及使用公共数据集。

### 5. 多参数成像

分析来自 CODEX、Vectra、MERFISH 和其他多重成像平台的空间蛋白质组学和基因表达数据。 PathML 提供专门的幻灯片类和转换，用于处理多参数数据、使用 Mesmer 进行细胞分割以及量化工作流程。

**请参阅：** `references/multiparametric.md` 了解 CODEX/Vectra 工作流程、细胞分割、标记定量以及与 AnnData 的集成。

### 6. 数据管理

使用 HDF5 格式高效存储和管理大型病理数据集。 PathML 在针对机器学习工作流程优化的统一存储结构中处理图块、掩码、元数据和提取的特征。

**请参阅：** `references/data_management.md` 了解 HDF5 集成、切片管理、数据集组织和批处理策略。
## 快速入门

### 安装

```bash
# Install PathML
uv pip install pathml

# With optional dependencies for all features
uv pip install pathml[all]
```

### 基本工作流程示例

<<<代码块_1>>>

### 常见工作流程

**H&E 图像分析：**
1. 用适当的滑动类加载 WSI
2. 应用组织检测和染色标准化
3. 执行核检测或训练分割模型
4. 提取特征并构建空间图
5. 进行下游分析

**多参数成像 (CODEX)：**
1. 使用 `CODEXSlide` 加载 CODEX 幻灯片
2. 折叠多运行通道数据
3. 使用 Mesmer 模型分割细胞
4. 量化标记表达
5. 导出至AnnData进行单细胞分析

**训练机器学习模型：**
1. 使用公共病理数据准备数据集
2. 使用 PathML 数据集创建 PyTorch DataLoader
3. 训练 HoVer-Net 或自定义模型
4. 评估保留的测试集
5. 使用 ONNX 进行部署以进行推理

## 详细文档参考

在处理特定任务时，请参阅相应的参考文件以获取全面信息：

- **加载图像：** `references/image_loading.md`
- **预处理工作流程：** `references/preprocessing.md`
- **空间分析：** `references/graphs.md`
- **模型训练：** `references/machine_learning.md`
- **CODEX/多路复用 IF：** `references/multiparametric.md`
- **数据存储：** `references/data_management.md`

## 资源

该技能包括按功能领域组织的综合参考文档。每个参考文件都包含特定 PathML 功能的详细 API 信息、工作流程示例、最佳实践和故障排除指南。

###参考资料/

深入介绍 PathML 功能的文档文件：

- `image_loading.md` - 整张幻灯片图像格式、加载策略、幻灯片类别
- `preprocessing.md` - 完整的转换目录、管道构建、预处理工作流程
- `graphs.md` - 图构建方法、特征提取、空间分析
- `machine_learning.md` - 模型架构、训练工作流程、评估、推理
- `multiparametric.md` - CODEX、Vectra、多重 IF 分析、细胞分割、定量
- `data_management.md` - HDF5 存储、切片管理、批处理、数据集组织

在处理特定的计算病理学任务时，根据需要加载这些参考文献。