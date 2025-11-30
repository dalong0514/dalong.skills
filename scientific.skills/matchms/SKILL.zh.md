<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：火柴人
描述：“质谱分析。处理 mzML/MGF/MSP、光谱相似性（余弦、修正余弦）、元数据协调、化合物 ID，用于代谢组学和 MS 数据处理。”
---

# 比赛

## 概述

Matchms 是一个用于质谱数据处理和分析的开源 Python 库。从各种格式导入光谱、标准化元数据、过滤峰值、计算光谱相似性并构建可重复的分析工作流程。

## 核心能力

### 1.导入和导出质谱数据

从多种文件格式加载光谱并导出处理后的数据：

```python
from matchms.importing import load_from_mgf, load_from_mzml, load_from_msp, load_from_json
from matchms.exporting import save_as_mgf, save_as_msp, save_as_json

# Import spectra
spectra = list(load_from_mgf("spectra.mgf"))
spectra = list(load_from_mzml("data.mzML"))
spectra = list(load_from_msp("library.msp"))

# Export processed spectra
save_as_mgf(spectra, "output.mgf")
save_as_json(spectra, "output.json")
```

**支持的格式：**
- mzML 和 mzXML（原始质谱格式）
- MGF（吉祥物通用格式）
- MSP（光谱库格式）
- JSON（GNPS 兼容）
- 代谢组学-USI 参考文献
- Pickle（Python 序列化）

有关详细的导入/导出文档，请参阅`references/importing_exporting.md`。

### 2. 频谱过滤和处理

应用全面的过滤器来标准化元数据并细化峰值数据：

<<<代码块_1>>>

**过滤类别：**
- **元数据处理**：协调化合物名称、导出化学结构、标准化加合物、纠正电荷
- **峰过滤**：标准化强度，按 m/z 或强度选择，删除前体峰
- **质量控制**：要求最小峰值，验证前体 m/z，确保元数据完整性
- **化学注释**：添加指纹，派生InChI/SMILES，修复结构不匹配

Matchms 提供 40 多个过滤器。有关完整的过滤器参考，请参阅`references/filtering.md`。

### 3.计算光谱相似度

使用各种相似性度量比较光谱：

<<<代码块_2>>>

**可用的相似度函数：**
- **CosineGreedy/CosineHungarian**：具有不同匹配算法的基于峰值的余弦相似度
- **ModifiedCosine**：考虑前体质量差异的余弦相似度
- **NeutralLossesCosine**：基于中性损失模式的相似性
- **FingerprintSimilarity**：使用指纹的分子结构相似性
- **MetadataMatch**：比较用户定义的元数据字段
- **PrecursorMzMatch/ParentMassMatch**：简单的基于质量的过滤

有关详细的相似性函数文档，请参阅`references/similarity.md`。

### 4. 构建处理管道

创建可重复的多步骤分析工作流程：

<<<代码块_3>>>

### 5. 使用 Spectrum 对象

核心 `Spectrum` 类包含质谱数据：

<<<代码块_4>>>

### 6.元数据管理

标准化和协调频谱元数据：

<<<代码块_5>>>

## 常见工作流程

对于典型的质谱分析工作流程，包括：
- 加载和预处理光谱库
- 将未知光谱与参考库进行匹配
- 质量过滤和数据清理
- 大规模相似度比较
- 基于网络的谱聚类

有关详细示例，请参阅`references/workflows.md`。

## 安装

<<<代码块_6>>>

对于分子结构处理（SMILES、InChI）：
```bash
uv pip install matchms[chemistry]
```

## 参考文档

详细的参考文档可在 `references/` 目录中找到：
- `filtering.md` - 带有描述的完整过滤器函数参考
- `similarity.md` - 所有相似性指标以及何时使用它们
- `importing_exporting.md` - 文件格式详细信息和 I/O 操作
- `workflows.md` - 常见分析模式和示例

根据需要加载这些参考资料，以获取有关特定匹配功能的详细信息。