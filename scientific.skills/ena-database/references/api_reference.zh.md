<!-- 此文件由机器翻译自 api_reference.md -->

# ENA API 参考

欧洲核苷酸档案 REST API 的综合参考。

## ENA 门户 API

**基本网址：** `https://www.ebi.ac.uk/ena/portal/api`

**官方文档：** https://www.ebi.ac.uk/ena/portal/api/doc

### 搜索端点

**端点：** `/search`

**方法：** 获取

**描述：** 使用灵活的过滤和格式化选项跨 ENA 数据类型执行高级搜索。

**参数：**

|参数|必填|描述 |示例|
|------------|----------|-------------|---------|
| `result` |是的 |要搜索的数据类型 | `sample`、`study`、`read_run`、`assembly`、`sequence`、`analysis`、`taxon` |
| `query` |是的 |使用 ENA 查询语法搜索查询 | `tax_eq(9606)`、`study_accession="PRJNA123456"` |
| `format` |没有 |输出格式（默认：tsv）| `json`、`tsv`、`xml` |
| `fields` |没有 |要返回的以逗号分隔的字段列表 | `accession,sample_title,scientific_name` |
| `limit` |没有 |最大结果数（默认：100000）| `10`、`1000` |
| `offset` |没有 |分页结果偏移 | `0`、`100` |
| `sortFields` |没有 |排序依据的字段（逗号分隔）| `accession`、`collection_date` |
| `sortOrder` |没有 |排序方向 | `asc`、`desc` |
| `dataPortal` |没有 |限制特定数据门户 | `ena`、`pathogen`、`metagenome` |
| `download` |没有 |触发文件下载 | `true`、`false` |
| `includeAccessions` |没有 |以逗号分隔的加入项包括 | `SAMN01,SAMN02` |
| `excludeAccessions` |没有 |要排除的以逗号分隔的加入 | `SAMN03,SAMN04` |

**查询语法：**

ENA 使用带有运算符的专用查询语言：

- **平等：** `field_name="value"` 或 `field_name=value`
- **通配符：** `field_name="*partial*"`（使用 * 表示通配符）
- **范围：** `field_name>=value AND field_name<=value`
- **逻辑：** `query1 AND query2`、`query1 OR query2`、`NOT query`
- **分类：** `tax_eq(taxon_id)` - 完全匹配，`tax_tree(taxon_id)` - 包括后代
- **日期范围：** `collection_date>=2020-01-01 AND collection_date<=2023-12-31`
- **在运算符中：** `study_accession IN (PRJNA1,PRJNA2,PRJNA3)`

**常见结果类型：**

- `study` - 研究项目/研究
- `sample` - 生物样本
- `read_run` - 原始测序运行
- `read_experiment` - 测序实验元数据
- `analysis` - 分析结果
- `assembly` - 基因组/转录组组装
- `sequence` - 组装序列
- `taxon` - 分类记录
- `coding` - 蛋白质编码序列
- `noncoding` - 非编码序列

**请求示例：**

```python
import requests

# Search for human samples
url = "https://www.ebi.ac.uk/ena/portal/api/search"
params = {
    "result": "sample",
    "query": "tax_eq(9606)",
    "format": "json",
    "fields": "accession,sample_title,collection_date",
    "limit": 100
}
response = requests.get(url, params=params)

# Search for RNA-seq experiments in a study
params = {
    "result": "read_experiment",
    "query": 'study_accession="PRJNA123456" AND library_strategy="RNA-Seq"',
    "format": "tsv"
}
response = requests.get(url, params=params)

# Find assemblies for E. coli with minimum contig N50
params = {
    "result": "assembly",
    "query": "tax_tree(562) AND contig_n50>=50000",
    "format": "json"
}
response = requests.get(url, params=params)
```

### 字段端点

**端点：** `/returnFields`

**方法：** 获取

**描述：** 列出特定结果类型的可用字段。

**参数：**

|参数|必填|描述 |示例|
|------------|----------|-------------|---------|
| `result` |是的 |数据类型 | `sample`、`study`、`assembly` |
| `dataPortal` |没有 |按数据门户过滤 | `ena`、`pathogen` |

**示例：**

<<<代码块_1>>>

### 结果端点

**端点：** `/results`

**方法：** 获取

**描述：** 列出可用的结果类型。

**示例：**

<<<代码块_2>>>

### 文件报告端点

**端点：** `/filereport`

**方法：** 获取

**描述：** 获取文件信息和下载 URL 以进行读取和分析。

**参数：**

|参数|必填|描述 |示例|
|------------|----------|-------------|---------|
| `accession` |是的 |运行或分析加入| `ERR123456` |
| `result` |是的 |必须是 `read_run` 或 `analysis` | `read_run` |
| `format` |没有 |输出格式 | `json`、`tsv` |
| `fields` |没有 |要包含的字段 | `run_accession,fastq_ftp,fastq_md5` |

**通用文件报告字段：**

- `run_accession` - 运行入藏号
- `fastq_ftp` - FASTQ 文件的 FTP URL（以分号分隔）
- `fastq_aspera` - FASTQ 文件的 Aspera URL
- `fastq_md5` - MD5 校验和（分号分隔）
- `fastq_bytes` - 文件大小（以字节为单位）（以分号分隔）
- `submitted_ftp` - 最初提交的文件的 FTP URL
- `sra_ftp` - SRA 格式文件的 FTP URL

**示例：**

<<<代码块_3>>>

## ENA 浏览器 API

**基本网址：** `https://www.ebi.ac.uk/ena/browser/api`

**官方文档：** https://www.ebi.ac.uk/ena/browser/api/doc

### XML 检索

**端点：** `/xml/{accession}`

**方法：** 获取

**描述：** 检索 XML 格式的记录元数据。

**参数：**

|参数|类型 |描述 |示例|
|------------|------|-------------|---------|
| `accession` |路径|记录入藏号 | `PRJNA123456`、`SAMEA123456`、`ERR123456` |
| `download` |查询 |设置为 `true` 以触发下载 | `true` |
| `includeLinks` |查询 |包括交叉引用链接 | `true`、`false` |

**示例：**

<<<代码块_4>>>

### 文本检索

**端点：** `/text/{accession}`

**方法：** 获取

**描述：** 检索 EMBL 平面文件格式的序列。

**参数：**

|参数|类型 |描述 |示例|
|------------|------|-------------|---------|
| `accession` |路径|序列加入 | `LN847353` |
| `download` |查询 |触发下载| `true` |
| `expandDataclasses` |查询 |包含相关数据类 | `true` |
| `lineLimit` |查询 |限制输出线 | `1000` |

**示例：**

<<<代码块_5>>>

### FASTA 检索

**端点：** `/fasta/{accession}`

**方法：** 获取

**描述：** 检索 FASTA 格式的序列。

**参数：**

|参数|类型 |描述 |示例|
|------------|------|-------------|---------|
| `accession` |路径|序列加入 | `LN847353` |
| `download` |查询 |触发下载| `true` |
| `range` |查询 |子序列范围 | `100-500` |
| `lineLimit` |查询 |限制输出线 | `1000` |

**示例：**

<<<代码块_6>>>

### 链接检索

**端点：** `/links/{source}/{accession}`

**方法：** 获取

**描述：** 获取外部数据库的交叉引用。

**参数：**

|参数|类型 |描述 |示例|
|------------|------|-------------|---------|
| `source` |路径|源数据库类型 | `sample`、`study`、`sequence` |
| `accession` |路径|入藏号| `SAMEA123456` |
| `target` |查询 |目标数据库过滤| `sra`、`biosample` |

**示例：**

```python
# Get all links for a sample
url = "https://www.ebi.ac.uk/ena/browser/api/links/sample/SAMEA123456"
response = requests.get(url)
```

## ENA 分类 REST API

**基本网址：** `https://www.ebi.ac.uk/ena/taxonomy/rest`

**描述：** 查询分类信息，包括谱系和等级。

### 税号查找

**端点：** `/tax-id/{taxon_id}`

**方法：** 获取

**描述：** 通过 NCBI 分类 ID 获取分类信息。

**示例：**

```python
# Get E. coli taxonomy
taxon_id = "562"
url = f"https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/{taxon_id}"
response = requests.get(url)
taxonomy = response.json()
# Returns: taxId, scientificName, commonName, rank, lineage, etc.
```

### 学名查找

**端点：** `/scientific-name/{name}`

**方法：** 获取

**描述：** 按学名搜索（可能返回多个匹配项）。

**示例：**

```python
# Search by scientific name
name = "Escherichia coli"
url = f"https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/{name}"
response = requests.get(url)
```

### 建议名字

**端点：** `/suggest-for-submission/{partial_name}`

**方法：** 获取

**描述：** 获取提交的分类建议（自动完成）。

**示例：**

```python
# Get suggestions
partial = "Escheri"
url = f"https://www.ebi.ac.uk/ena/taxonomy/rest/suggest-for-submission/{partial}"
response = requests.get(url)
```

## 交叉参考服务

**基本网址：** `https://www.ebi.ac.uk/ena/xref/rest`

**描述：** 访问外部数据库中与 ENA 条目相关的记录。

### 获取交叉引用

**端点：** `/json/{source}/{accession}`

**方法：** 获取

**描述：** 检索 JSON 格式的交叉引用。

**参数：**

|参数|类型 |描述 |示例|
|------------|------|-------------|---------|
| `source` |路径|来源数据库| `ena`、`sra` |
| `accession` |路径|入藏号| `SRR000001` |

**示例：**

```python
# Get cross-references for an SRA accession
url = "https://www.ebi.ac.uk/ena/xref/rest/json/sra/SRR000001"
response = requests.get(url)
xrefs = response.json()
```

## CRAM 参考注册表

**基本网址：** `https://www.ebi.ac.uk/ena/cram`

**描述：** 检索 CRAM 文件中使用的参考序列。

### MD5 查找

**端点：** `/md5/{md5_checksum}`

**方法：** 获取

**描述：** 通过MD5校验和检索参考序列。

**示例：**

```python
# Get reference by MD5
md5 = "7c3f69f0c5f0f0de6d7c34e7c2e25f5c"
url = f"https://www.ebi.ac.uk/ena/cram/md5/{md5}"
response = requests.get(url)
reference_fasta = response.text
```

## 速率限制和错误处理

**速率限制：**
- 最大：每秒 50 个请求
- 超出限制返回 HTTP 429（请求过多）
- 收到 429 响应时实施指数退避
**常见 HTTP 状态代码：**

- `200 OK` - 成功
- `204 No Content` - 成功但没有返回数据
- `400 Bad Request` - 无效参数
- `404 Not Found` - 未找到加入
- `429 Too Many Requests` - 超出速率限制
- `500 Internal Server Error` - 服务器错误（退避重试）

**错误处理模式：**

```python
import time
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

def create_session_with_retries():
    """Create requests session with retry logic"""
    session = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"]
    )
    adapter = HTTPAdapter(max_retries=retries)
    session.mount("https://", adapter)
    return session

# Usage
session = create_session_with_retries()
response = session.get(url, params=params)
```

## 批量下载推荐

用于下载大量文件或大型数据集：

1. **直接使用FTP**而不是API进行文件下载
   - 基本 FTP：`ftp://ftp.sra.ebi.ac.uk/vol1/fastq/`
   - Aspera 高速：`era-fasp@fasp.sra.ebi.ac.uk:`

2. **使用 enaBrowserTools** 命令行实用程序
   ```bash
   # Download by accession
   enaDataGet ERR123456

   # Download all runs from a study
   enaGroupGet PRJEB1234
   ```

3. **批量 API 请求**并具有适当的延迟
   ```python
   import time

   accessions = ["ERR001", "ERR002", "ERR003"]
   for acc in accessions:
       response = requests.get(f"{base_url}/xml/{acc}")
       # Process response
       time.sleep(0.02)  # 50 req/sec = 0.02s between requests
   ```

## 查询优化技巧

1. **使用特定的结果类型**而不是广泛的搜索
2. **使用 `fields` 参数将字段**限制为您需要的内容
3. **对大型结果集使用分页**（限制+偏移量）
4. **在本地缓存分类查找**
5. **尽可能选择 JSON/TSV** 而不是 XML（更小、更快）
6. **使用includeAccessions/excludeAccessions**高效过滤大型结果集
7. **尽可能将相似的查询批量处理**在一起