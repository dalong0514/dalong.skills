<!-- 此文件由机器翻译自 api_reference.md -->

# ClinVar API 和数据访问参考

## 概述

ClinVar 提供了多种编程数据访问方法：
- **电子实用程序** - NCBI 用于搜索和检索数据的 REST API
- **Entrez Direct** - 用于 UNIX 环境的命令行工具
- **FTP 下载** - XML、VCF 和制表符分隔格式的批量数据文件
- **提交 API** - 用于提交变体解释的 REST API

## 电子实用程序 API

### 基本网址
```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
```

### 支持的操作

#### 1. esearch - 搜索记录
使用与 Web 界面相同的查询语法搜索 ClinVar。

**端点：**
<<<代码块_1>>>

**参数：**
- `db=clinvar` - 数据库名称（必填）
- `term=<query>` - 搜索查询（必填）
- `retmax=<N>` - 返回的最大记录数（默认值：20）
- `retmode=json` - 返回格式（json 或 xml）
- `usehistory=y` - 将大型数据集的结果存储在服务器上

**查询示例：**
<<<代码块_2>>>

**常用搜索字段：**
- `[gene]` - 基因符号
- `[CLNSIG]` - 临床意义（致病性、良性等）
- `[disorder]` - 疾病/病症名称
- `[variant name]` - HGVS 表达式或变体标识符
- `[chr]` - 染色体编号
- `[Assembly]` - GRCh37 或 GRCh38

#### 2. esummary - 检索记录摘要
获取特定 ClinVar 记录的摘要信息。

**端点：**
<<<代码块_3>>>

**参数：**
- `db=clinvar` - 数据库名称（必填）
- `id=<UIDs>` - 以逗号分隔的 ClinVar UID 列表
- `retmode=json` - 返回格式（json 或 xml）
- `version=2.0` - API 版本（推荐用于 JSON）

**示例：**
<<<代码块_4>>>

**摘要输出包括：**
- 加入（RCV/VCV）
- 临床意义
- 审核状态
- 基因符号
- 变体类型
- 基因组位置（GRCh37 和 GRCh38）
- 相关条件
- 等位基因起源（种系/体细胞）

#### 3. efetch - 检索完整记录
下载完整的 XML 记录以进行详细分析。

**端点：**
<<<代码块_5>>>

**参数：**
- `db=clinvar` - 数据库名称（必填）
- `id=<UIDs>` - 逗号分隔的 ClinVar UID
- `rettype=vcv` 或 `rettype=rcv` - 记录类型

**示例：**
<<<代码块_6>>>

#### 4. elink - 查找相关记录
将 ClinVar 记录链接到其他 NCBI 数据库。

**端点：**
```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi
```

**可用链接：**
- clinvar_pubmed - PubMed 引文链接
- clinvar_gene - 链接到基因数据库
- clinvar_medgen - 链接到 MedGen（条件）
- clinvar_snp - 链接到 dbSNP

**示例：**
```bash
# Find PubMed articles for a variant
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=clinvar&db=pubmed&id=12345"
```

### 工作流程示例：完整的搜索和检索

```bash
# Step 1: Search for variants
SEARCH_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=CFTR[gene]+AND+pathogenic[CLNSIG]&retmode=json&retmax=10"

# Step 2: Parse IDs from search results
# (Extract id list from JSON response)

# Step 3: Retrieve summaries
SUMMARY_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=<ids>&retmode=json&version=2.0"

# Step 4: Fetch full records if needed
FETCH_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&id=<ids>&rettype=vcv"
```

## Entrez Direct（命令行）

安装 Entrez Direct 以进行命令行访问：
```bash
sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
```

### 常用命令

**搜索：**
```bash
esearch -db clinvar -query "BRCA1[gene] AND pathogenic[CLNSIG]"
```

**管道搜索到摘要：**
```bash
esearch -db clinvar -query "TP53[gene]" | \
  efetch -format docsum | \
  xtract -pattern DocumentSummary -element AccessionVersion Title
```

**计数结果：**
```bash
esearch -db clinvar -query "breast cancer[disorder]" | \
  efilter -status reviewed | \
  efetch -format docsum
```

## 速率限制和最佳实践

### 速率限制
- **没有 API 密钥：** 3 个请求/秒
- **使用 API 密钥：** 10 个请求/秒
- 大型数据集：使用`usehistory=y`以避免重复查询

### API 密钥设置
1. 在 https://www.ncbi.nlm.nih.gov/account/ 注册 NCBI 帐户
2. 在账户设置中生成API密钥
3. 在所有请求中添加`&api_key=<YOUR_KEY>`

### 最佳实践
- 在自动化之前测试 Web 界面上的查询
- 对于大型结果集（>500 条记录）使用 `usehistory`
- 对速率限制错误实施指数退避
- 适当时缓存结果
- 使用批量请求而不是单独查询
- 尊重 NCBI 服务器 - 不要在美国高峰时段提交大型作业

## 使用 Biopython 的 Python 示例

```python
from Bio import Entrez

# Set email (required by NCBI)
Entrez.email = "your.email@example.com"

# Search ClinVar
def search_clinvar(query, retmax=100):
    handle = Entrez.esearch(db="clinvar", term=query, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

# Get summaries
def get_summaries(id_list):
    ids = ",".join(id_list)
    handle = Entrez.esummary(db="clinvar", id=ids, retmode="json")
    record = Entrez.read(handle)
    handle.close()
    return record

# Example usage
variant_ids = search_clinvar("BRCA2[gene] AND pathogenic[CLNSIG]")
summaries = get_summaries(variant_ids)
```

## 错误处理

### 常见 HTTP 状态代码
- `200` - 成功
- `400` - 错误请求（检查查询语法）
- `429` - 请求过多（速率受限）
- `500` - 服务器错误（使用指数退避重试）

### 错误响应示例
```xml
<ERROR>Empty id list - nothing to do</ERROR>
```

## 其他资源

- NCBI 电子实用程序文档：https://www.ncbi.nlm.nih.gov/books/NBK25501/
- ClinVar 网络服务：https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/
- Entrez Direct 食谱：https://www.ncbi.nlm.nih.gov/books/NBK179288/