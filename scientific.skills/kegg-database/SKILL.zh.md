<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：kegg数据库
描述：“直接 REST API 访问 KEGG（仅限学术用途）。通路分析、基因通路映射、代谢通路、药物相互作用、ID 转换。对于具有多个数据库的 Python 工作流程，更喜欢生物服务。使用它进行直接 HTTP/REST 工作或 KEGG 特定控制。”
---

# KEGG 数据库

## 概述

KEGG（京都基因和基因组百科全书）是用于生物途径分析和分子相互作用网络的综合生物信息学资源。

**重要**：KEGG API 仅供学术用户学术使用。

## 何时使用此技能

当使用 KEGG 的 REST API 查询多个生物体的通路、基因、化合物、酶、疾病和药物时，应使用此技能。

## 快速入门

该技能提供：
1. 用于所有 KEGG REST API 操作的 Python 辅助函数 (`scripts/kegg_api.py`)
2. 全面的参考文档 (`references/kegg_reference.md`) 以及详细的 API 规范

当用户请求 KEGG 数据时，确定需要哪个操作并使用 `scripts/kegg_api.py` 中的适当函数。

## 核心运营

### 1. 数据库信息 (`kegg_info`)

检索有关 KEGG 数据库的元数据和统计数据。

**何时使用**：了解数据库结构、检查可用数据、获取发布信息。

**用法**：
```python
from scripts.kegg_api import kegg_info

# Get pathway database info
info = kegg_info('pathway')

# Get organism-specific info
hsa_info = kegg_info('hsa')  # Human genome
```

**常用数据库**：`kegg`、`pathway`、`module`、`brite`、`genes`、`genome`、`compound`、`glycan`、`reaction`、 `enzyme`、`disease`、`drug`

### 2. 列出条目 (`kegg_list`)

列出 KEGG 数据库中的条目标识符和名称。

**何时使用**：获取生物体的所有途径，列出基因，检索化合物目录。

**用法**：
<<<代码块_1>>>

**常见生物体代码**：`hsa`（人类）、`mmu`（小鼠）、`dme`（果蝇）、`sce`（酵母）、`eco`（大肠杆菌）

### 3. 搜索 (`kegg_find`)

通过关键字或分子特性搜索 KEGG 数据库。

**何时使用**：按名称/描述查找基因，按分子式或质量搜索化合物，按关键字发现条目。

**用法**：
<<<代码块_2>>>

**搜索选项**：`formula`（完全匹配）、`exact_mass`（范围）、`mol_weight`（范围）

### 4. 检索条目 (`kegg_get`)

获取完整的数据库条目或特定的数据格式。

**何时使用**：检索通路详细信息、获取基因/蛋白质序列、下载通路图、访问化合物结构。

**用法**：
<<<代码块_3>>>

**输出格式**：`aaseq`（蛋白质 FASTA）、`ntseq`（核苷酸 FASTA）、`mol`（MOL 格式）、`kcf`（KCF 格式）、`image` (PNG)、`kgml` (XML)、 `json`（路径 JSON）

**重要**：图像、KGML 和 JSON 格式一次仅允许一项。

### 5. ID 转换（`kegg_conv`）

在 KEGG 和外部数据库之间转换标识符。

**何时使用**：将 KEGG 数据与其他数据库集成、映射基因 ID、转换化合物标识符。

**用法**：
<<<代码块_4>>>

**支持的转换**：`ncbi-geneid`、`ncbi-proteinid`、`uniprot`、`pubchem`、`chebi`

### 6. 交叉引用 (`kegg_link`)

查找 KEGG 数据库内部和之间的相关条目。

**何时使用**：查找包含基因的通路、获取通路中的基因、将基因映射到 KO 组、查找通路中的化合物。

**用法**：
<<<代码块_5>>>

**常见链接**：基因↔途径、途径↔化合物、途径↔酶、基因↔ko（直系同源）

### 7. 药物间相互作用 (`kegg_ddi`)

检查药物间的相互作用。

**何时使用**：分析药物组合、检查禁忌症、药理学研究。

**用法**：
<<<代码块_6>>>

## 常见分析工作流程

### 工作流程 1：基因到通路图谱

**用例**：寻找与感兴趣的基因相关的通路（例如，用于通路富集分析）。

```python
from scripts.kegg_api import kegg_find, kegg_link, kegg_get

# Step 1: Find gene ID by name
gene_results = kegg_find('genes', 'p53')

# Step 2: Link gene to pathways
pathways = kegg_link('pathway', 'hsa:7157')  # TP53 gene

# Step 3: Get detailed pathway information
for pathway_line in pathways.split('\n'):
    if pathway_line:
        pathway_id = pathway_line.split('\t')[1].replace('path:', '')
        pathway_info = kegg_get(pathway_id)
        # Process pathway information
```

### 工作流程 2：途径丰富背景

**用例**：获取生物体途径中的所有基因以进行富集分析。

```python
from scripts.kegg_api import kegg_list, kegg_link

# Step 1: List all human pathways
pathways = kegg_list('pathway', 'hsa')

# Step 2: For each pathway, get associated genes
for pathway_line in pathways.split('\n'):
    if pathway_line:
        pathway_id = pathway_line.split('\t')[0]
        genes = kegg_link('genes', pathway_id)
        # Process genes for enrichment analysis
```

### 工作流程 3：化合物到通路分析

**用例**：寻找含有感兴趣化合物的代谢途径。
```python
from scripts.kegg_api import kegg_find, kegg_link, kegg_get

# Step 1: Search for compound
compound_results = kegg_find('compound', 'glucose')

# Step 2: Link compound to reactions
reactions = kegg_link('reaction', 'cpd:C00031')  # Glucose

# Step 3: Link reactions to pathways
pathways = kegg_link('pathway', 'rn:R00299')  # Specific reaction

# Step 4: Get pathway details
pathway_info = kegg_get('map00010')  # Glycolysis
```

### 工作流程 4：跨数据库集成

**用例**：将 KEGG 数据与 UniProt、NCBI 或 PubChem 数据库集成。

```python
from scripts.kegg_api import kegg_conv, kegg_get

# Step 1: Convert KEGG gene IDs to external database IDs
uniprot_map = kegg_conv('uniprot', 'hsa')
ncbi_map = kegg_conv('ncbi-geneid', 'hsa')

# Step 2: Parse conversion results
for line in uniprot_map.split('\n'):
    if line:
        kegg_id, uniprot_id = line.split('\t')
        # Use external IDs for integration

# Step 3: Get sequences using KEGG
sequence = kegg_get('hsa:10458', 'aaseq')
```

### 工作流程 5：生物体特异性途径分析

**用例**：比较不同生物体的途径。

```python
from scripts.kegg_api import kegg_list, kegg_get

# Step 1: List pathways for multiple organisms
human_pathways = kegg_list('pathway', 'hsa')
mouse_pathways = kegg_list('pathway', 'mmu')
yeast_pathways = kegg_list('pathway', 'sce')

# Step 2: Get reference pathway for comparison
ref_pathway = kegg_get('map00010')  # Reference glycolysis

# Step 3: Get organism-specific versions
hsa_glycolysis = kegg_get('hsa00010')
mmu_glycolysis = kegg_get('mmu00010')
```

## 衔接课程类别

KEGG 将路径分为七大类。在解释路径 ID 或向用户推荐路径时：

1. **代谢**（例如，`map00010` - 糖酵解，`map00190` - 氧化磷酸化）
2. **遗传信息处理**（例如，`map03010` - 核糖体，`map03040` - 剪接体）
3. **环境信息处理**（例如，`map04010` - MAPK 信号传导，`map02010` - ABC 转运蛋白）
4. **细胞过程**（例如，`map04140` - 自噬，`map04210` - 细胞凋亡）
5. **有机体系统**（例如，`map04610` - 补体级联，`map04910` - 胰岛素信号传导）
6. **人类疾病**（例如，`map05200` - 癌症通路，`map05010` - 阿尔茨海默病）
7. **药物开发**（按时间顺序和基于目标的分类）

参考`references/kegg_reference.md`了解详细的通路列表和分类。

## 重要的标识符和格式

### 路径 ID
- `map#####` - 参考途径（通用，非生物体特异性）
- `hsa#####` - 人类路径
- `mmu#####` - 鼠标通路

### 基因 ID
- 格式：`organism:gene_number`（例如，`hsa:10458`）

### 化合物 ID
- 格式：`cpd:C#####`（例如，对于 ATP，`cpd:C00002`）

### 药物 ID
- 格式：`dr:D#####`（例如，`dr:D00001`）

### 酶 ID
- 格式：`ec:EC_number`（例如，`ec:1.1.1.1`）

### KO (KEGG Orthology) ID
- 格式：`ko:K#####`（例如，`ko:K00001`）

## API 限制

使用 KEGG API 时请遵守这些约束：

1. **条目限制**：每次操作最多 10 个条目（image/kgml/json 除外：仅限 1 个条目）
2. **学术用途**：API仅供学术用途；商业用途需要许可
3. **HTTP状态码**：检查200（成功）、400（错误请求）、404（未找到）
4. **速率限制**：没有明确的限制，但避免快速请求

## 详细参考

有关全面的 API 文档、数据库规范、有机体代码和高级用法，请参阅`references/kegg_reference.md`。这包括：

- KEGG数据库的完整列表
- 详细的API操作语法
- 所有生物体代码
- HTTP 状态代码和错误处理
- 与 Biopython 和 R/Bioconductor 集成
- API使用的最佳实践

## 故障排除

**404 Not Found**：条目或数据库不存在；验证 ID 和有机体代码
**400 Bad Request**：API调用中的语法错误；检查参数格式
**空结果**：搜索词可能与条目不匹配；尝试更广泛的关键字
**图像/KGML 错误**：这些格式仅适用于单个条目；删除批处理

## 附加工具

对于交互式路径可视化和注释：
- **KEGG 映射器**：https://www.kegg.jp/kegg/mapper/
- **BlastKOALA**：自动化基因组注释
- **GhostKOALA**：元基因组/元转录组注释