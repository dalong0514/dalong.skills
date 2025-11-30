<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：pyopenms
描述：OpenMS 的 Python 接口，用于质谱数据分析。用于 LC-MS/MS 蛋白质组学和代谢组学工作流程，包括文件处理（mzML、mzXML、mzTab、FASTA、pepXML、protXML、mzIdentML）、信号处理、特征检测、肽鉴定和定量分析。适用于处理质谱数据、分析蛋白质组学实验或处理代谢组学数据集。
---

#PyOpenMS

## 概述

PyOpenMS 提供与 OpenMS 库的 Python 绑定以进行计算质谱分析，从而能够分析蛋白质组学和代谢组学数据。用于处理质谱文件格式、处理光谱数据、检测特征、识别肽/蛋白质以及执行定量分析。

## 安装

使用 uv 安装：

```bash
uv uv pip install pyopenms
```

验证安装：

<<<代码块_1>>>

## 核心能力

PyOpenMS 将功能组织到以下领域：

### 1. 文件 I/O 和数据格式

处理质谱文件格式并在表示形式之间进行转换。

**支持的格式**：mzML、mzXML、TraML、mzTab、FASTA、pepXML、protXML、mzIdentML、featureXML、consensusXML、idXML

基本文件读取：

<<<代码块_2>>>

**有关详细的文件处理**：请参阅 `references/file_io.md`

### 2. 信号处理

通过平滑、过滤、质心和归一化处理原始光谱数据。

基本频谱处理：

<<<代码块_3>>>

**有关算法详细信息**：请参阅`references/signal_processing.md`

### 3.特征检测

检测并链接光谱和样品的特征以进行定量分析。

<<<代码块_4>>>

**完整工作流程**：请参阅 `references/feature_detection.md`

### 4. 肽和蛋白质鉴定

与搜索引擎集成并处理识别结果。

**支持的引擎**：Comet、Mascot、MSGFPlus、XTandem、OMSSA、Myrimatch

基本识别工作流程：

<<<代码块_5>>>

**详细工作流程**：请参阅`references/identification.md`

### 5.代谢组学分析

执行非靶向代谢组学预处理和分析。

典型工作流程：
1. 加载并处理原始数据
2. 检测特征
3. 调整样品的保留时间
4. 将特征链接到共识图
5. 使用复合数据库进行注释

**有关完整的代谢组学工作流程**：请参阅 `references/metabolomics.md`

## 数据结构

PyOpenMS 使用这些主要对象：

- **MSExperiment**：光谱和色谱图的收集
- **MSSpectrum**：具有 m/z 和强度对的单一质谱图
- **MS色谱图**：色谱图
- **功能**：检测到的色谱峰和质量指标
- **FeatureMap**：特征集合
- **肽鉴定**：肽的搜索结果
- **蛋白质识别**：蛋白质搜索结果

**有关详细文档**：请参阅 `references/data_structures.md`

## 常见工作流程

### 快速入门：加载和探索数据

<<<代码块_6>>>

### 参数管理

大多数算法使用参数系统：

```python
# Get algorithm parameters
algo = ms.GaussFilter()
params = algo.getParameters()

# View available parameters
for param in params.keys():
    print(f"{param}: {params.getValue(param)}")

# Modify parameters
params.setValue("gaussian_width", 0.2)
algo.setParameters(params)
```

### 导出到 Pandas

将数据转换为 pandas DataFrame 进行分析：

```python
import pyopenms as ms
import pandas as pd

# Load feature map
fm = ms.FeatureMap()
ms.FeatureXMLFile().load("features.featureXML", fm)

# Convert to DataFrame
df = fm.get_df()
print(df.head())
```

## 与其他工具集成

PyOpenMS 集成：
- **Pandas**：将数据导出到 DataFrames
- **NumPy**：使用峰值数组
- **Scikit-learn**：MS 数据上的机器学习
- **Matplotlib/Seaborn**：可视化
- **R**：通过 rpy2 桥

## 资源

- **官方文档**：https://pyopenms.readthedocs.io
- **OpenMS 文档**：https://www.openms.org
- **GitHub**：https://github.com/OpenMS/OpenMS

## 参考文献

- `references/file_io.md` - 全面的文件格式处理
- `references/signal_processing.md` - 信号处理算法
- `references/feature_detection.md` - 特征检测和链接
- `references/identification.md` - 肽和蛋白质鉴定
- `references/metabolomics.md` - 代谢组学特定工作流程
- `references/data_structures.md` - 核心对象和数据结构