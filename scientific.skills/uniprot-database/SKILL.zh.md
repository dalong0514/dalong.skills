<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：uniprot 数据库
描述：“直接 REST API 访问 UniProt。蛋白质搜索、FASTA 检索、ID 映射、Swiss-Prot/TrEMBL。对于具有多个数据库的 Python 工作流程，更喜欢生物服务（40 多个服务的统一接口）。使用它进行直接 HTTP/REST 工作或 UniProt 特定控制。”
---

# UniProt 数据库

## 概述

UniProt 是世界领先的综合蛋白质序列和功能信息资源。按名称、基因或登录号搜索蛋白质，以 FASTA 格式检索序列，跨数据库执行 ID 映射，通过 REST API 访问 Swiss-Prot/TrEMBL 注释以进行蛋白质分析。

## 何时使用此技能

该技能应该在以下情况下使用：
- 按名称、基因符号、登录名或生物体搜索蛋白质条目
- 检索 FASTA 或其他格式的蛋白质序列
- UniProt 和外部数据库（Ensembl、RefSeq、PDB 等）之间的映射标识符
- 访问蛋白质注释，包括 GO 术语、结构域和功能描述
- 高效批量检索多个蛋白质条目
- 查询已审查的 (Swiss-Prot) 与未审查的 (TrEMBL) 蛋白质数据
- 流式传输大型蛋白质数据集
- 使用特定于字段的搜索语法构建自定义查询

## 核心能力

### 1. 寻找蛋白质

使用自然语言查询或结构化搜索语法搜索 UniProt。

**常见搜索模式：**
```python
# Search by protein name
query = "insulin AND organism_name:\"Homo sapiens\""

# Search by gene name
query = "gene:BRCA1 AND reviewed:true"

# Search by accession
query = "accession:P12345"

# Search by sequence length
query = "length:[100 TO 500]"

# Search by taxonomy
query = "taxonomy_id:9606"  # Human proteins

# Search by GO term
query = "go:0005515"  # Protein binding
```

使用 API 搜索端点：`https://rest.uniprot.org/uniprotkb/search?query={query}&format={format}`

**支持的格式：** JSON、TSV、Excel、XML、FASTA、RDF、TXT

### 2. 检索单个蛋白质条目

按登录号检索特定蛋白质条目。

**入藏号格式：**
- 经典：P12345、Q1AAA9、O15530（6 个字符：字母 + 5 个字母数字）
- 扩展：A0A022YWF9（新条目为 10 个字符）

**检索端点：** `https://rest.uniprot.org/uniprotkb/{accession}.{format}`

示例：`https://rest.uniprot.org/uniprotkb/P12345.fasta`

### 3.批量检索和ID映射

在不同数据库系统之间映射蛋白质标识符并有效检索多个条目。

**ID 映射工作流程：**
1. 将映射作业提交至：`https://rest.uniprot.org/idmapping/run`
2. 检查作业状态：`https://rest.uniprot.org/idmapping/status/{jobId}`
3. 检索结果：`https://rest.uniprot.org/idmapping/results/{jobId}`

**支持的映射数据库：**
- UniProtKB AC/ID
- 基因名称
- Ensembl、RefSeq、EMBL
- PDB、AlphaFoldDB
- KEGG、GO 条款
- 还有更多（请参阅`/references/id_mapping_databases.md`）

**限制：**
- 每个作业最多 100,000 个 ID
- 结果保存 7 天

### 4. 流式传输大型结果集

对于超出分页限制的大型查询，请使用流端点：

`https://rest.uniprot.org/uniprotkb/stream?query={query}&format={format}`

流端点返回所有结果，不分页，适合下载完整数据集。

### 5. 自定义检索字段

准确指定要检索哪些字段以实现高效的数据传输。

**常用字段：**
- `accession` - UniProt 登录号
- `id` - 条目名称
- `gene_names` - 基因名称
- `organism_name` - 有机体
- `protein_name` - 蛋白质名称
- `sequence` - 氨基酸序列
- `length` - 序列长度
- `go_*` - 基因本体注释
- `cc_*` - 注释字段（功能、交互等）
- `ft_*` - 特征注释（域、站点等）

**示例：** `https://rest.uniprot.org/uniprotkb/search?query=insulin&fields=accession,gene_names,organism_name,length,sequence&format=tsv`

有关完整字段列表，请参阅 `/references/api_fields.md`。

## Python 实现

对于编程访问，请使用提供的帮助程序脚本 `scripts/uniprot_client.py` 来实现：

- `search_proteins(query, format)` - 使用任何查询搜索 UniProt
- `get_protein(accession, format)` - 检索单个蛋白质条目
- `map_ids(ids, from_db, to_db)` - 标识符类型之间的映射
- `batch_retrieve(accessions, format)` - 检索多个条目
- `stream_results(query, format)` - 流式传输大型结果集

**替代Python包：**
- **Unipressed**：UniProt REST API 的现代类型 Python 客户端
- **bioservices**：综合生物信息学 Web 服务客户端

## 查询语法示例

**布尔运算符：**
<<<代码块_1>>>

**特定领域的搜索：**
<<<代码块_2>>>

**范围查询：**
<<<代码块_3>>>

**通配符：**
<<<代码块_4>>>

请参阅 `/references/query_syntax.md` 以获取全面的语法文档。

## 最佳实践

1. **尽可能使用经过审核的条目**：使用 `reviewed:true` 过滤 Swiss-Prot（手动策划）条目
2. **明确指定格式**：选择最合适的格式（FASTA用于序列，TSV用于表格数据，JSON用于编程解析）
3. **使用字段选择**：仅请求您需要的字段，以减少带宽和处理时间
4. **处理分页**：对于大型结果集，实现适当的分页或使用流端点
5. **缓存结果**：将频繁访问的数据存储在本地，以最大程度地减少 API 调用
6. **速率限制**：尊重API资源；对大批量操作实施延迟
7. **检查数据质量**：TrEMBL 条目是计算预测； Swiss-Prot 条目经过人工审核

## 资源

### 脚本/
`uniprot_client.py` - Python 客户端，具有用于常见 UniProt 操作的辅助函数，包括搜索、检索、ID 映射和流式传输。

###参考资料/
- `api_fields.md` - 用于自定义查询的可用字段的完整列表
- `id_mapping_databases.md` - ID 映射操作支持的数据库
- `query_syntax.md` - 带有高级示例的全面查询语法
- `api_examples.md` - 多种语言的代码示例（Python、curl、R）

## 其他资源

- **API 文档**：https://www.uniprot.org/help/api
- **交互式 API 资源管理器**：https://www.uniprot.org/api-documentation
- **REST 教程**：https://www.uniprot.org/help/uniprot_rest_tutorial
- **查询语法帮助**：https://www.uniprot.org/help/query-fields
- **SPARQL 端点**：https://sparql.uniprot.org/（用于高级图形查询）