<!-- 此文件由机器翻译自 services_reference.md -->

# BioServices：完整服务参考

本文档为 BioServices 中提供的所有主要服务提供了全面的参考，包括关键方法、参数和用例。

## 蛋白质和基因资源

### UniProt

蛋白质序列和功能信息数据库。

**初始化：**
```python
from bioservices import UniProt
u = UniProt(verbose=False)
```

**关键方法：**

- `search(query, frmt="tab", columns=None, limit=None, sort=None, compress=False, include=False, **kwargs)`
  - 使用灵活的查询语法搜索 UniProt
  - `frmt`：“tab”、“fasta”、“xml”、“rdf”、“gff”、“txt”
  - `columns`：逗号分隔列表（例如，“id,genes,organism,length”）
  - 返回：请求格式的字符串

- `retrieve(uniprot_id, frmt="txt")`
  - 检索特定的 UniProt 条目
  - `frmt`：“txt”、“fasta”、“xml”、“rdf”、“gff”
  - 返回：按请求格式输入数据

- `mapping(fr="UniProtKB_AC-ID", to="KEGG", query="P43403")`
  - 在数据库之间转换标识符
  - `fr`/`to`：数据库标识符（参见identifier_mapping.md）
  - `query`：单个 ID 或逗号分隔列表
  - 返回：将输入映射到输出 ID 的字典

- `searchUniProtId(pattern, columns="entry name,length,organism", limit=100)`
  - 基于 ID 的搜索的便捷方法
  - 返回：制表符分隔值

**常用列：** id、条目名称、基因、生物体、蛋白质名称、长度、序列、go-id、ec、途径、相互作用子

**使用案例：**
- BLAST 的蛋白质序列检索
- 功能注释查找
- 跨数据库标识符映射
- 批量蛋白质信息检索

---

### KEGG（京都基因和基因组百科全书）

代谢途径、基因和生物体数据库。

**初始化：**
<<<代码块_1>>>

**关键方法：**

- `list(database)`
  - 列出 KEGG 数据库中的条目
  - `database`：“有机体”、“途径”、“模块”、“疾病”、“药物”、“化合物”
  - 返回：带有条目的多行字符串

- `find(database, query)`
  - 通过关键字搜索数据库
  - 返回：具有 ID 的匹配条目列表

- `get(entry_id)`
  - 通过ID检索条目
  - 支持基因、途径、化合物等。
  - 返回：原始输入文本

- `parse(data)`
  - 将 KEGG 条目解析到字典中
  - 返回：带有结构化数据的字典

- `lookfor_organism(name)`
  - 按名称模式搜索生物体
  - 返回：匹配的有机体代码列表

- `lookfor_pathway(name)`
  - 按名称搜索路径
  - 返回：路径 ID 列表

- `get_pathway_by_gene(gene_id, organism)`
  - 寻找含有基因的途径
  - 返回：路径 ID 列表

- `parse_kgml_pathway(pathway_id)`
  - 解析交互路径 KGML
  - 返回：带有“条目”和“关系”的字典

- `pathway2sif(pathway_id)`
  - 提取简单交互格式数据
  - 激活/抑制过滤器
  - 返回：交互元组列表

**有机体代码：**
- hsa：智人
- mmu：小家鼠
- dme：黑腹果蝇
- sce：酿酒酵母
- 生态：大肠杆菌

**使用案例：**
- 路径分析和可视化
- 基因功能注释
- 代谢网络重建
- 蛋白质-蛋白质相互作用提取

---

### HGNC（人类基因命名委员会）

官方人类基因命名机构。

**初始化：**
<<<代码块_2>>>

**关键方法：**
- `search(query)`：搜索基因符号/名称
- `fetch(format, query)`：检索基因信息

**使用案例：**
- 标准化人类基因名称
- 查找官方基因符号

---

### 我的基因信息

基因注释和查询服务。

**初始化：**
<<<代码块_3>>>

**关键方法：**
- `querymany(ids, scopes, fields, species)`：批量基因查询
- `getgene(geneid)`：获取基因注释

**使用案例：**
- 批量基因注释检索
- 基因ID转换

---

## 化合物资源

### ChEBI（具有生物价值的化学实体）

分子实体词典。

**初始化：**
<<<代码块_4>>>

**关键方法：**
- `getCompleteEntity(chebi_id)`：完整的化合物信息
- `getLiteEntity(chebi_id)`：基本信息
- `getCompleteEntityByList(chebi_ids)`：批量检索

**使用案例：**
- 小分子信息
- 化学结构数据
- 复合属性查找

---

### 化学BL

生物活性药物样化合物数据库。

**初始化：**
<<<代码块_5>>>

**关键方法：**
- `get_molecule_form(chembl_id)`：化合物详细信息
- `get_target(chembl_id)`：目标信息
- `get_similarity(chembl_id)`：获取给定的相似化合物 
- `get_assays()`：生物测定数据

**使用案例：**
- 药物发现数据
- 查找相似的化合物  
- 生物活性信息
- 目标化合物关系

---

### 联合化学

化学标识符映射服务。

**初始化：**
<<<代码块_6>>>

**关键方法：**
- `get_compound_id_from_kegg(kegg_id)`：KEGG → ChEMBL
- `get_all_compound_ids(src_compound_id, src_id)`：获取所有 ID
- `get_src_compound_ids(src_compound_id, from_src_id, to_src_id)`：转换 ID

**来源 ID：**
- 1：化学分子生物学
- 2：药物银行
- 3：PDB
- 6: 凯格
- 7：ChEBI
- 22：PubChem

**使用案例：**
- 跨数据库复合ID映射
- 连接化学数据库

---

### 公共化学

NIH 的化合物数据库。

**初始化：**
```python
from bioservices import PubChem
p = PubChem()
```

**关键方法：**
- `get_compounds(identifier, namespace)`：检索化合物
- `get_properties(properties, identifier, namespace)`：获取属性

**使用案例：**
- 化学结构检索
- 复合属性信息

---

## 序列分析工具

### NCBI 爆炸

序列相似性搜索。

**初始化：**
```python
from bioservices import NCBIblast
s = NCBIblast(verbose=False)
```

**关键方法：**
- `run(program, sequence, stype, database, email, **params)`
  - 提交 BLAST 作业
  - `program`：“blastp”、“blastn”、“blastx”、“tblastn”、“tblastx”
  - `stype`：“蛋白质”或“DNA”
  - `database`：“uniprotkb”、“pdb”、“refseq_ Protein”等。
  - `email`：NCBI 要求
  - 返回：工作 ID

- `getStatus(jobid)`
  - 检查工作状态
  - 返回：“正在运行”、“已完成”、“错误”

- `getResult(jobid, result_type)`
  - 检索结果
  - `result_type`：“out”（默认）、“ids”、“xml”

**重要提示：** BLAST 作业是异步的。在检索结果之前始终检查状态。

**使用案例：**
- 蛋白质同源性搜索
- 序列相似性分析
- 通过同源性进行功能注释

---

## 途径和互动资源

### 反应组

路径数据库。

**初始化：**
```python
from bioservices import Reactome
r = Reactome()
```

**关键方法：**
- `get_pathway_by_id(pathway_id)`：路径详细信息
- `search_pathway(query)`：搜索路径

**使用案例：**
- 人类通路分析
- 生物过程注释

---

### 灵能

蛋白质相互作用查询服务（联合 30 多个数据库）。

**初始化：**
```python
from bioservices import PSICQUIC
s = PSICQUIC()
```

**关键方法：**
- `query(database, query_string)`
  - 查询具体交互数据库
  - 返回：PSI-MI TAB 格式

- `activeDBs`
  - 房产列表可用数据库
  - 返回：数据库名称列表

**可用数据库：** MINT、IntAct、BioGRID、DIP、InnateDB、MatrixDB、MPIDB、UniProt 等 30 多个

**查询语法：** 支持AND、OR、物种过滤器
- 示例：“ZAP70 和物种：9606”

**使用案例：**
- 蛋白质-蛋白质相互作用发现
- 网络分析
- 相互作用组图谱

---

### 完整复合体

蛋白质复合物数据库。

**初始化：**
```python
from bioservices import IntactComplex
i = IntactComplex()
```

**关键方法：**
- `search(query)`：搜索复合体
- `details(complex_ac)`：复杂的细节

**使用案例：**
- 蛋白质复合物成分
- 多蛋白组装分析

---

### OmniPath

综合信号通路数据库。

**初始化：**
```python
from bioservices import OmniPath
o = OmniPath()
```

**关键方法：**
- `interactions(datasets, organisms)`：获取交互
- `ptms(datasets, organisms)`：翻译后修饰

**使用案例：**
- 细胞信号分析
- 监管网络映射

---

## 基因本体论

### 快走

基因本体注释服务。

**初始化：**
```python
from bioservices import QuickGO
g = QuickGO()
```

**关键方法：**
- `Term(go_id, frmt="obo")`
  - 检索GO术语信息
  - 返回：术语定义和元数据

- `Annotation(protein=None, goid=None, format="tsv")`
  - 获取GO注释
  - 返回：按请求格式的注释

**GO类别：**
- 生物过程（BP）
- 分子功能（MF）
- 蜂窝组件 (CC)

**使用案例：**
- 功能注释
- 富集分析
- GO术语查找

---

## 基因组资源

### 生物玛特

基因组数据的数据挖掘工具。

**初始化：**
```python
from bioservices import BioMart
b = BioMart()
```

**关键方法：**
- `datasets(dataset)`：列出可用数据集
- `attributes(dataset)`：列出属性
- `query(query_xml)`：执行 BioMart 查询

**使用案例：**
- 批量基因组数据检索
- 自定义基因组注释
- SNP信息

---

### 数组表达

基因表达数据库。

**初始化：**
```python
from bioservices import ArrayExpress
a = ArrayExpress()
```

**关键方法：**
- `queryExperiments(keywords)`：搜索实验
- `retrieveExperiment(accession)`：获取实验数据

**使用案例：**
- 基因表达数据
- 微阵列分析
- RNA-seq数据检索

---

### ENA（欧洲核苷酸档案）

核苷酸序列数据库。

**初始化：**
```python
from bioservices import ENA
e = ENA()
```

**关键方法：**
- `search_data(query)`：搜索序列
- `retrieve_data(accession)`：检索序列

**使用案例：**
- 核苷酸序列检索
- 基因组组装访问

---

## 结构生物学

### PDB（蛋白质数据库）

3D 蛋白质结构数据库。

**初始化：**
```python
from bioservices import PDB
p = PDB()
```

**关键方法：**
- `get_file(pdb_id, file_format)`：下载结构文件
- `search(query)`：搜索结构

**文件格式：** pdb、cif、xml

**使用案例：**
- 3D结构检索
- 基于结构的分析
- PyMOL 可视化

---

### 普法姆

蛋白质家族数据库。

**初始化：**
```python
from bioservices import Pfam
p = Pfam()
```

**关键方法：**
- `searchSequence(sequence)`：按顺序查找域
- `getPfamEntry(pfam_id)`：域信息

**使用案例：**
- 蛋白质结构域鉴定
- 家庭分类
- 功能基序发现

---

## 专业资源

### 生物模型

系统生物学模型库。

**初始化：**
```python
from bioservices import BioModels
b = BioModels()
```

**关键方法：**
- `get_model_by_id(model_id)`：检索 SBML 模型

**使用案例：**
- 系统生物学建模
- SBML模型检索

---

### COG（直系同源基因簇）

直系同源基因分类。

**初始化：**
```python
from bioservices import COG
c = COG()
```

**使用案例：**
- 直系同源分析
- 功能分类

---

### BiGG 模型

代谢网络模型。

**初始化：**
```python
from bioservices import BiGG
b = BiGG()
```

**关键方法：**
- `list_models()`：可用型号
- `get_model(model_id)`：型号详细信息

**使用案例：**
- 代谢网络分析
- 通量平衡分析

---

## 一般模式

### 错误处理

所有服务都可能抛出异常。将调用包装在 try- except 中：

```python
try:
    result = service.method(params)
    if result:
        # Process result
        pass
except Exception as e:
    print(f"Error: {e}")
```

### 详细控制

大多数服务支持 `verbose` 参数：
```python
service = Service(verbose=False)  # Suppress HTTP logs
```

### 速率限制

服务有超时和速率限制：
```python
service.TIMEOUT = 30  # Adjust timeout
service.DELAY = 1     # Delay between requests (if supported)
```

### 输出格式

常用格式参数：
- `frmt`：“xml”、“json”、“tab”、“txt”、“fasta”
- `format`：特定于服务的变体

### 缓存

一些服务缓存结果：
```python
service.CACHE = True  # Enable caching
service.clear_cache()  # Clear cache
```

## 其他资源

详细API文档：
- 官方文档：https://bioservices.readthedocs.io/
- 从主页链接的个人服务文档
- 源代码：https://github.com/cokelaer/bioservices