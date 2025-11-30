<!-- 此文件由机器翻译自 kegg_reference.md -->

# KEGG 数据库参考

## 概述

KEGG（京都基因和基因组百科全书）是一个综合性生物信息学资源，维护手动策划的通路图和分子相互作用网络。它提供了“分子相互作用、反应和关系的接线图”来理解生物系统。

**基本网址**：https://rest.kegg.jp
**官方文档**：https://www.kegg.jp/kegg/rest/keggapi.html
**访问限制**：KEGG API 仅供学术用户用于学术用途。

## KEGG 数据库

KEGG 集成了 16 个主要数据库，分为系统信息、基因组信息、化学信息和健康信息类别：

### 系统信息
- **PATHWAY**：手动绘制新陈代谢、遗传信息处理、环境信息处理、细胞过程、有机体系统、人类疾病和药物开发的路径图
- **模块**：功能单元和路径构建块
- **BRITE**：层次分类和本体

### 基因组信息
- **GENOME**：带有注释的完整基因组
- **GENES**：所有生物体的基因目录
- **ORTHOLOGY**：直系同源组（KO：KEGG Orthology）
- **SSDB**：序列相似性数据库

### 化学信息
- **化合物**：代谢物和其他化学物质
- **GLYCAN**：聚糖结构
- **反应**：化学反应
- **RCLASS**：反应类别（化学结构转换模式）
- **酶**：酶命名法
- **网络**：网络变化

### 健康信息
- **疾病**：与遗传和环境因素有关的人类疾病
- **药物**：已批准的药物，具有化学结构和靶点信息
- **DGROUP**：药物组

### 外部数据库链接
KEGG 对外部数据库的交叉引用包括：
- **PubMed**：文献参考
- **NCBI 基因**：基因数据库
- **UniProt**：蛋白质序列
- **PubChem**：化合物
- **ChEBI**：具有生物学意义的化学实体

## REST API 操作

### 1. INFO - 数据库元数据

**语法**：`/info/<database>`

检索数据库的发布信息和统计信息。

**示例**：
- `/info/kegg` - KEGG 系统信息
- `/info/pathway` - 路径数据库信息
- `/info/hsa` - 人体组织信息

### 2. LIST - 条目列表

**语法**：`/list/<database>[/<organism>]`

列出条目标识符和关联名称。

**参数**：
- `database` - 数据库名称（途径、酶、基因等）或条目 (hsa:10458)
- `organism` - 可选生物体代码（例如，hsa 代表人类，eco 代表大肠杆菌）

**示例**：
- `/list/pathway` - 所有参考路径
- `/list/pathway/hsa` - 人类特异性途径
- `/list/hsa:10458+ece:Z5100` - 特定基因条目（最多 10 个）

**有机体代码**：三个或四个字母代码
- `hsa` - 智人（人类）
- `mmu` - 小家鼠（小鼠）
- `dme` - 黑腹果蝇（果蝇）
- `sce` - 酿酒酵母（酵母）
- `eco` - 大肠杆菌 K-12 MG1655

### 3. FIND - 搜索条目

**语法**：`/find/<database>/<query>[/<option>]`

按关键字或分子属性搜索条目。

**参数**：
- `database` - 要搜索的数据库
- `query` - 搜索词或分子属性
- `option` - 可选：`formula`、`exact_mass`、`mol_weight`

**搜索字段**（取决于数据库）：
- 条目、名称、符号、基因名称、描述、定义
- 有机体、分类学、直系学、通路等。

**示例**：
- `/find/genes/shiga toxin` - 基因中的关键字搜索
- `/find/compound/C7H10N4O2/formula` - 精确公式匹配
- `/find/drug/300-310/exact_mass` - 质量范围搜索 (300-310 Da)
- `/find/compound/300-310/mol_weight` - 分子量范围

### 4. GET - 检索条目

**语法**：`/get/<entry>[+<entry>...][/<option>]`

检索完整的数据库条目或特定的数据格式。

**参数**：
- `entry` - 条目 ID（最多 10 个，用 + 连接）
- `option` - 输出格式（可选）

**输出选项**：
- `aaseq` - 氨基酸序列 (FASTA)
- `ntseq` - 核苷酸序列 (FASTA)
- `mol` - MOL 格式（化合物/药物）
- `kcf` - KCF 格式（KEGG 化学函数、化合物/药物）
- `image` - PNG 图像（路径图，仅限单个条目）
- `kgml` - KGML XML（路径结构，仅限单个条目）
- `json` - JSON 格式（仅限路径，仅限单个条目）

**示例**：
- `/get/hsa00010` - 糖酵解途径（人类）
- `/get/hsa:10458+ece:Z5100` - 多个基因（最多 10 个）
- `/get/hsa:10458/aaseq` - 蛋白质序列
- `/get/cpd:C00002` - ATP 化合物条目
- `/get/hsa05130/json` - JSON 形式的癌症通路
- `/get/hsa05130/image` - PNG 路径图

**图像限制**：仅允许使用图像选项输入一项

### 5. CONV - ID转换

**语法**：`/conv/<target_db>/<source_db>`

在 KEGG 和外部数据库之间转换标识符。

**支持的转换**：
- `ncbi-geneid` ↔ KEGG 基因
- `ncbi-proteinid` ↔ KEGG 基因
- `uniprot` ↔ KEGG 基因
- `pubchem` ↔ KEGG 化合物/药物
- `chebi` ↔ KEGG 化合物/药物

**示例**：
- `/conv/ncbi-geneid/hsa` - 所有人类基因到 NCBI 基因 ID
- `/conv/hsa/ncbi-geneid` - 人类基因的 NCBI 基因 ID（反向）
- `/conv/uniprot/hsa:10458` - UniProt 的特定基因
- `/conv/pubchem/compound` - 所有化合物均为 PubChem ID

### 6. 链接 - 交叉引用

**语法**：`/link/<target_db>/<source_db>`

查找 KEGG 数据库内部和之间的相关条目。

**常用链接**：
- 基因 ↔ 途径
- 途径 ↔ 化合物
- 途径 ↔ 酶
- 基因 ↔ 直系同源 (KO)
- 化合物 ↔ 反应

**示例**：
- `/link/pathway/hsa` - 与人类基因相关的所有途径
- `/link/genes/hsa00010` - 糖酵解途径中的基因
- `/link/pathway/hsa:10458` - 包含特定基因的路径
- `/link/compound/hsa00010` - 途径中的化合物

### 7. DDI - 药物间相互作用

**语法**：`/ddi/<drug>[+<drug>...]`

检索从日本药品标签中提取的药物相互作用信息。

**参数**：
- `drug` - 药品条目 ID（最多 10 个，用 + 连接）

**示例**：
- `/ddi/D00001` - 单一药物的相互作用
- `/ddi/D00001+D00002` - 多种药物之间的相互作用

## 通路分类

KEGG 将路径分为七大类：

### 1.新陈代谢
碳水化合物、能量、脂质、核苷酸、氨基酸、聚糖生物合成和代谢、辅因子和维生素代谢、萜类和聚酮化合物代谢、次级代谢物生物合成、外源生物降解

**途径示例**：
- `map00010` - 糖酵解/糖异生
- `map00020` - 柠檬酸盐循环（TCA 循环）
- `map00190` - 氧化磷酸化

### 2.遗传信息处理
转录、翻译、折叠/排序/降解、复制和修复

**途径示例**：
- `map03010` - 核糖体
- `map03020` - RNA 聚合酶
- `map03040` - 剪接体

### 3.环境信息处理
膜运输、信号转导

**途径示例**：
- `map02010` - ABC 传输器
- `map04010` - MAPK 信号通路

### 4. 细胞过程
运输和分解代谢、细胞生长和死亡、细胞群落、细胞运动

**途径示例**：
- `map04140` - 自噬
- `map04210` - 细胞凋亡

### 5. 有机体系统
免疫、内分泌、循环、消化、神经、感觉、发育、环境适应

**途径示例**：
- `map04610` - 补体和凝血级联
- `map04910` - 胰岛素信号通路

### 6.人类疾病
癌症、免疫疾病、神经退行性疾病、心血管疾病、代谢疾病、传染病

**途径示例**：
- `map05200` - 癌症通路
- `map05010` - 阿尔茨海默病

### 7. 药物开发
按时间顺序分类和基于目标分类

## 通用标识符和命名

### 路径 ID
- `map#####` - 参考途径（通用）
- `hsa#####` - 人类特异性途径
- `mmu#####` - 小鼠特异性途径
- 格式：有机体代码+5位数字

### 基因 ID
- `hsa:10458` - 人类基因（有机体：gene_id）
- 格式：生物体代码+冒号+基因号

### 化合物 ID
- `cpd:C00002` - ATP
- 格式：cpd:C#####

### 药物 ID
- `dr:D00001` - 药物输入
- 格式：dr:D#####

### 酶 ID
- `ec:1.1.1.1` - 乙醇脱氢酶
- 格式：ec:EC_number

### KO (KEGG Orthology) ID
- `ko:K00001` - 直系同源组
- 格式：ko:K#####

## API 限制和最佳实践

### 速率限制和限制
- 每次操作最多 10 个条目（图像/kgml 除外：1 个条目）
- 仅限学术用途 - 商业用途需要单独许可
- 没有明确的速率限制记录，但避免快速请求

### HTTP 状态代码
- `200` - 成功
- `400` - 错误请求（查询中的语法错误）
- `404` - 未找到（条目或数据库不存在）

### 最佳实践
1. 始终检查响应中的 HTTP 状态代码
2. 对于批量操作，使用+批量输入（最多10个）
3、本地缓存结果，减少API调用
4. 尽可能使用特定的生物体代码以获得更快的结果
5. 对于路径可视化，使用 Web 界面或 KGML/JSON 格式
6. 仔细解析制表符分隔的输出（跨操作的格式一致）

## 与其他工具集成

### Biopython 集成
Biopython 提供了 `Bio.KEGG.REST` 模块以方便 Python 集成：
```python
from Bio.KEGG import REST
result = REST.kegg_list("pathway").read()
```

### KEGGREST（R/Bioconductor）
R 用户可以使用 KEGGREEST 包：
<<<代码块_1>>>

## 常见分析工作流程

### 工作流程 1：基因到通路图谱
1. 从您的生物体中获取基因 ID
2. 使用`/link/pathway/<gene_id>`查找关联路径
3. 使用`/get/<pathway_id>`检索详细的通路信息

### 工作流程 2：途径丰富背景
1. 使用`/list/pathway/<org>`获取所有生物体路径
2. 使用`/link/genes/<pathway_id>`获取每个pathway中的基因
3. 进行统计富集分析

### 工作流程 3：化合物到反应图谱
1. 使用`/find/compound/<name>`查找化合物ID
2. 使用 `/link/reaction/<compound_id>` 查找反应
3. 使用`/link/pathway/<reaction_id>`查找包含反应的路径

### 工作流程 4：集成的 ID 转换
1. 使用`/conv/uniprot/<org>`将KEGG基因映射到UniProt
2. 使用 `/conv/ncbi-geneid/<org>` 映射到 NCBI 基因 ID
3.使用转换后的ID与其他数据库集成

## 其他资源

- **KEGG 映射器**：https://www.kegg.jp/kegg/mapper/ - 交互式路径映射
- **BlastKOALA**：测序基因组的自动注释
- **GhostKOALA**：元基因组和元转录组的注释
- **KEGG 模块**：https://www.kegg.jp/kegg/module.html
- **KEGG Brite**：https://www.kegg.jp/kegg/brite.html