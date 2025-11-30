<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：生物服务
描述：“用于 40 多种生物信息学服务的主要 Python 工具。多数据库工作流程的首选：UniProt、KEGG、ChEMBL、PubChem、Reactome、QuickGO。用于查询、ID 映射、路径分析的统一 API。要直接 REST 控制，请使用单独的数据库技能（uniprot-database、kegg-database）。”
---

# 生物服务

## 概述

BioServices 是一个 Python 包，提供对大约 40 个生物信息学 Web 服务和数据库的编程访问。检索生物数据、执行跨数据库查询、映射标识符、分析序列以及在 Python 工作流程中集成多种生物资源。该包透明地处理 REST 和 SOAP/WSDL 协议。

## 何时使用此技能

该技能应该在以下情况下使用：
- 从 UniProt、PDB、Pfam 检索蛋白质序列、注释或结构
- 通过KEGG或Reactome分析代谢途径和基因功能
- 搜索化合物数据库（ChEBI、ChEMBL、PubChem）以获取化学信息
- 在不同生物数据库之间转换标识符（KEGG↔UniProt、化合物 ID）
- 运行序列相似性搜索（BLAST、MUSCLE 比对）
- 查询基因本体术语（QuickGO、GO注释）
- 访问蛋白质-蛋白质相互作用数据（PSICQUIC、IntactComplex）
- 挖掘基因组数据（BioMart、ArrayExpress、ENA）
- 在单个工作流程中集成来自多个生物信息学资源的数据

## 核心能力

### 1. 蛋白质分析

检索蛋白质信息、序列和功能注释：

```python
from bioservices import UniProt

u = UniProt(verbose=False)

# Search for protein by name
results = u.search("ZAP70_HUMAN", frmt="tab", columns="id,genes,organism")

# Retrieve FASTA sequence
sequence = u.retrieve("P43403", "fasta")

# Map identifiers between databases
kegg_ids = u.mapping(fr="UniProtKB_AC-ID", to="KEGG", query="P43403")
```

**关键方法：**
- `search()`：使用灵活的搜索词查询 UniProt
- `retrieve()`：获取各种格式的蛋白质条目（FASTA、XML、选项卡）
- `mapping()`：在数据库之间转换标识符

参考：`references/services_reference.md` 了解完整的 UniProt API 详细信息。

### 2. 通路发现与分析

访问基因和生物体的 KEGG 通路信息：

<<<代码块_1>>>

**关键方法：**
- `lookfor_organism()`、`lookfor_pathway()`：按名称搜索
- `get_pathway_by_gene()`：查找包含基因的通路
- `parse_kgml_pathway()`：提取结构化路径数据
- `pathway2sif()`：获取蛋白质相互作用网络

参考：`references/workflow_patterns.md` 了解完整的通路分析工作流程。

### 3. 化合物数据库搜索

跨多个数据库搜索和交叉引用化合物：

<<<代码块_2>>>

**通用工作流程：**
1. 在KEGG中按名称搜索化合物
2.提取KEGG化合物ID
3. 使用 UniChem 进行 KEGG → ChEMBL 作图
4. KEGG 条目中通常提供 ChEBI ID

参考：`references/identifier_mapping.md` 完整的跨数据库映射指南。

### 4. 序列分析

运行 BLAST 搜索和序列比对：

<<<代码块_3>>>

**注意：** BLAST 作业是异步的。在检索结果之前检查状态。

### 5. 标识符映射

在不同生物数据库之间转换标识符：

<<<代码块_4>>>

**支持的映射 (UniProt)：**
- UniProtKB ↔ KEGG
- UniProtKB ↔ Ensebl
- UniProtKB ↔ PDB
- UniProtKB ↔ RefSeq
- 还有更多（请参阅`references/identifier_mapping.md`）

### 6. 基因本体查询

访问 GO 术语和注释：

<<<代码块_5>>>

### 7. 蛋白质-蛋白质相互作用

通过 PSICQUIC 查询交互数据库：

<<<代码块_6>>>

**可用数据库：** MINT、IntAct、BioGRID、DIP 和 30 多个其他数据库。

## 多服务集成工作流程

BioServices 擅长结合多种服务进行全面分析。常见的集成模式：

### 完整的蛋白质分析流程

执行完整的蛋白质表征工作流程：

```bash
python scripts/protein_analysis_workflow.py ZAP70_HUMAN your.email@example.com
```

该脚本演示了：
1.UniProt搜索蛋白质条目
2. FASTA序列检索
3. BLAST相似度搜索
4. KEGG通路发现
5. PSICQUIC交互映射

### 通路网络分析

分析生物体的所有途径：

```bash
python scripts/pathway_analysis.py hsa output_directory/
```

摘录与分析：
- 生物体的所有途径 ID
- 每个途径的蛋白质-蛋白质相互作用
- 交互类型分布
- 导出为 CSV/SIF 格式

### 跨数据库化合物搜索

跨数据库映射化合物标识符：

```bash
python scripts/compound_cross_reference.py Geldanamycin
```

检索：
- KEGG 化合物 ID
- ChEBI 标识符
- ChEMBL 标识符
- 基本化合物特性

### 批量标识符转换

一次转换多个标识符：

```bash
python scripts/batch_id_converter.py input_ids.txt --from UniProtKB_AC-ID --to KEGG
```

## 最佳实践

### 输出格式处理

不同的服务以不同的格式返回数据：
- **XML**：使用 BeautifulSoup 进行解析（大多数 SOAP 服务）
- **制表符分隔 (TSV)**：用于表格数据的 Pandas DataFrames
- **字典/JSON**：直接Python操作
- **FASTA**：用于序列分析的 BioPython 集成

### 速率限制和冗长

控制API请求行为：

```python
from bioservices import KEGG

k = KEGG(verbose=False)  # Suppress HTTP request details
k.TIMEOUT = 30  # Adjust timeout for slow connections
```

### 错误处理

将服务调用包装在 try- except 块中：

```python
try:
    results = u.search("ambiguous_query")
    if results:
        # Process results
        pass
except Exception as e:
    print(f"Search failed: {e}")
```

### 有机体代码

使用标准有机体缩写：
- `hsa`：智人（人类）
- `mmu`：小家鼠（小鼠）
- `dme`：黑腹果蝇
- `sce`：酿酒酵母（酵母）

列出所有生物：`k.list("organism")` 或 `k.organismIds`

### 与其他工具集成

BioServices 可以很好地配合：
- **BioPython**：对检索到的 FASTA 数据进行序列分析
- **Pandas**：表格数据操作
- **PyMOL**：3D 结构可视化（检索 PDB ID）
- **NetworkX**：通路相互作用的网络分析
- **Galaxy**：工作流程平台的自定义工具包装器

## 资源

### 脚本/

演示完整工作流程的可执行 Python 脚本：

- `protein_analysis_workflow.py`：端到端蛋白质表征
- `pathway_analysis.py`：KEGG 通路发现和网络提取
- `compound_cross_reference.py`：多数据库复合搜索
- `batch_id_converter.py`：批量标识符映射实用程序

脚本可以直接执行，也可以针对特定用例进行调整。

###参考资料/

根据需要加载详细文档：

- `services_reference.md`：所有 40 多个服务及其方法的综合列表
- `workflow_patterns.md`：详细的多步骤分析工作流程
- `identifier_mapping.md`：跨数据库ID转换完整指南

处理特定服务或复杂集成任务时加载参考。

## 安装

```bash
uv pip install bioservices
```

依赖关系是自动管理的。包在 Python 3.9-3.12 上进行了测试。

## 附加信息

详细的API文档和高级功能请参考：
- 官方文档：https://bioservices.readthedocs.io/
- 源代码：https://github.com/cokelaer/bioservices
- `references/services_reference.md` 中特定于服务的引用