<!-- 此文件由机器翻译自 api_reference.md -->

# RCSB PDB API 参考

本文档提供有关 RCSB 蛋白质数据库 API 的详细信息，包括高级使用模式、数据模式和最佳实践。

## API 概述

RCSB PDB 提供多种编程接口：

1. **数据 API** - 当您有标识符时检索 PDB 数据
2. **搜索 API** - 查找匹配特定搜索条件的标识符
3. **ModelServer API** - 访问大分子模型子集
4. **VolumeServer API** - 检索体积数据子集
5. **序列坐标 API** - 获取结构和序列数据库之间的比对
6. **Alignment API** - 执行结构对齐计算

## 数据API

### 核心数据对象

数据 API 按层次结构组织信息：

- **core_entry**：PDB条目或计算结构模型（CSM ID以AF_或MA_开头）
- **core_polymer_entity**：蛋白质、DNA 和 RNA 实体
- **核心_非聚合物_实体**：配体、辅因子、离子
- **核心支化实体**：低聚糖
- **core_ assembly**：生物组件
- **core_polymer_entity_instance**：单个链
- **core_chem_comp**：化学成分

### REST API 端点

基本网址：`https://data.rcsb.org/rest/v1/`

**输入数据：**
```
GET https://data.rcsb.org/rest/v1/core/entry/{entry_id}
```

**聚合物实体：**
<<<代码块_1>>>

**组装：**
<<<代码块_2>>>

**示例：**
<<<代码块_3>>>

### GraphQL API

端点：`https://data.rcsb.org/graphql`

GraphQL API 支持灵活的数据检索，允许您在单个查询中从层次结构的任何级别获取任何数据。

**查询示例：**
<<<代码块_4>>>

**Python 示例：**
<<<代码块_5>>>

### 通用数据字段

**入门级：**
- `struct.title` - 结构标题/描述
- `exptl[].method` - 实验方法（X射线衍射、核磁共振、电子显微镜等）
- `rcsb_entry_info.resolution_combined` - 分辨率（以 Ångströms 为单位）
- `rcsb_entry_info.deposited_atom_count` - 原子总数
- `rcsb_accession_info.deposit_date` - 沉积日期
- `rcsb_accession_info.initial_release_date` - 发布日期

**聚合物实体级别：**
- `entity_poly.pdbx_seq_one_letter_code` - 主序列
- `rcsb_polymer_entity.formula_weight` - 分子量
- `rcsb_entity_source_organism.scientific_name` - 源生物体
- `rcsb_entity_source_organism.ncbi_taxonomy_id` - NCBI 分类 ID

**装配级别：**
- `rcsb_assembly_info.polymer_entity_count` - 聚合物实体数量
- `rcsb_assembly_info.assembly_id` - 程序集标识符

## 搜索API

### 查询类型

搜索 API 支持七种主要查询类型：

1. **TextQuery** - 全文搜索
2. **AttributeQuery** - 基于属性的搜索
3. **SequenceQuery** - 序列相似度搜索
4. **SequenceMotifQuery** - Motif 模式搜索
5. **StructSimilarityQuery** - 3D结构相似度
6. **StructMotifQuery** - 结构主题搜索
7. **ChemSimilarityQuery** - 化学相似性搜索

### 属性查询运算符

AttributeQuery 的可用运算符：

- `exact_match` - 精确字符串匹配
- `contains_words` - 包含所有单词
- `contains_phrase` - 包含确切的短语
- `equals` - 数值相等
- `greater` - 大于（数字）
- `greater_or_equal` - 大于或等于
- `less` - 小于（数字）
- `less_or_equal` - 小于或等于
- `range` - 数值范围（闭区间）
- `exists` - 字段有一个值
- `in` - 列表中的值

### 常见的可搜索属性

**分辨率和质量：**
<<<代码块_6>>>

**实验方法：**
```python
from rcsbapi.search.attrs import exptl

query = AttributeQuery(
    attribute=exptl.method,
    operator="exact_match",
    value="X-RAY DIFFRACTION"
)
```

**生物体：**
```python
from rcsbapi.search.attrs import rcsb_entity_source_organism

query = AttributeQuery(
    attribute=rcsb_entity_source_organism.scientific_name,
    operator="exact_match",
    value="Homo sapiens"
)
```

**分子量：**
```python
from rcsbapi.search.attrs import rcsb_polymer_entity

query = AttributeQuery(
    attribute=rcsb_polymer_entity.formula_weight,
    operator="range",
    value=(10000, 50000)  # 10-50 kDa
)
```

**发布日期：**
```python
from rcsbapi.search.attrs import rcsb_accession_info

# Structures released in 2024
query = AttributeQuery(
    attribute=rcsb_accession_info.initial_release_date,
    operator="range",
    value=("2024-01-01", "2024-12-31")
)
```

### 序列相似度搜索

使用 MMseqs2 搜索具有相似序列的结构：

```python
from rcsbapi.search import SequenceQuery

# Basic sequence search
query = SequenceQuery(
    value="MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM",
    evalue_cutoff=0.1,
    identity_cutoff=0.9
)

# With sequence type specified
query = SequenceQuery(
    value="ACGTACGTACGT",
    evalue_cutoff=1e-5,
    identity_cutoff=0.8,
    sequence_type="dna"  # or "rna" or "protein"
)
```

### 结构相似性搜索

使用 BioZernike 查找具有相似 3D 几何形状的结构：

```python
from rcsbapi.search import StructSimilarityQuery

# Search by entry
query = StructSimilarityQuery(
    structure_search_type="entry",
    entry_id="4HHB"
)

# Search by chain
query = StructSimilarityQuery(
    structure_search_type="chain",
    entry_id="4HHB",
    chain_id="A"
)

# Search by assembly
query = StructSimilarityQuery(
    structure_search_type="assembly",
    entry_id="4HHB",
    assembly_id="1"
)
```

### 组合查询

使用Python按位运算符组合查询：

```python
from rcsbapi.search import TextQuery, AttributeQuery
from rcsbapi.search.attrs import rcsb_entry_info, rcsb_entity_source_organism

# AND operation (&)
query1 = TextQuery("kinase")
query2 = AttributeQuery(
    attribute=rcsb_entity_source_organism.scientific_name,
    operator="exact_match",
    value="Homo sapiens"
)
combined = query1 & query2

# OR operation (|)
organism1 = AttributeQuery(
    attribute=rcsb_entity_source_organism.scientific_name,
    operator="exact_match",
    value="Homo sapiens"
)
organism2 = AttributeQuery(
    attribute=rcsb_entity_source_organism.scientific_name,
    operator="exact_match",
    value="Mus musculus"
)
combined = organism1 | organism2

# NOT operation (~)
all_structures = TextQuery("protein")
low_res = AttributeQuery(
    attribute=rcsb_entry_info.resolution_combined,
    operator="greater",
    value=3.0
)
high_res_only = all_structures & (~low_res)

# Complex combinations
high_res_human_kinases = (
    TextQuery("kinase") &
    AttributeQuery(
        attribute=rcsb_entity_source_organism.scientific_name,
        operator="exact_match",
        value="Homo sapiens"
    ) &
    AttributeQuery(
        attribute=rcsb_entry_info.resolution_combined,
        operator="less",
        value=2.5
    )
)
```

### 返回类型

控制返回哪些信息：

```python
from rcsbapi.search import TextQuery, ReturnType

query = TextQuery("hemoglobin")

# Return PDB IDs (default)
results = list(query())  # ['4HHB', '1A3N', ...]

# Return entry IDs with scores
results = list(query(return_type=ReturnType.ENTRY, return_scores=True))
# [{'identifier': '4HHB', 'score': 0.95}, ...]

# Return polymer entities
results = list(query(return_type=ReturnType.POLYMER_ENTITY))
# ['4HHB_1', '4HHB_2', ...]
```

## 文件下载地址

### 结构文件

**PDB 格式（旧版）：**
```
https://files.rcsb.org/download/{PDB_ID}.pdb
```

**mmCIF 格式（现代标准）：**
```
https://files.rcsb.org/download/{PDB_ID}.cif
```

**结构因素：**
```
https://files.rcsb.org/download/{PDB_ID}-sf.cif
```

**生物组装：**
```
https://files.rcsb.org/download/{PDB_ID}.pdb1  # Assembly 1
https://files.rcsb.org/download/{PDB_ID}.pdb2  # Assembly 2
```

**FASTA 序列：**
```
https://www.rcsb.org/fasta/entry/{PDB_ID}
```

### Python 下载助手

```python
import requests

def download_pdb_file(pdb_id, format="pdb", output_dir="."):
    """
    Download PDB structure file.

    Args:
        pdb_id: 4-character PDB ID
        format: 'pdb' or 'cif'
        output_dir: Directory to save file
    """
    base_url = "https://files.rcsb.org/download"
    url = f"{base_url}/{pdb_id}.{format}"

    response = requests.get(url)
    if response.status_code == 200:
        output_path = f"{output_dir}/{pdb_id}.{format}"
        with open(output_path, "w") as f:
            f.write(response.text)
        print(f"Downloaded {pdb_id}.{format}")
        return output_path
    else:
        print(f"Error downloading {pdb_id}: {response.status_code}")
        return None

# Usage
download_pdb_file("4HHB", format="pdb")
download_pdb_file("4HHB", format="cif")
```

## 速率限制和最佳实践

### 速率限制

- API实施速率限制以确保公平使用
- 如果超出限制，您将收到 429 HTTP 错误代码
- 建议的起点：每秒几个请求
- 使用指数退避来找到可接受的请求率

### 指数退避实现

```python
import time
import requests

def fetch_with_retry(url, max_retries=5, initial_delay=1):
    """
    Fetch URL with exponential backoff on rate limit errors.

    Args:
        url: URL to fetch
        max_retries: Maximum number of retry attempts
        initial_delay: Initial delay in seconds
    """
    delay = initial_delay

    for attempt in range(max_retries):
        response = requests.get(url)

        if response.status_code == 200:
            return response
        elif response.status_code == 429:
            print(f"Rate limited. Waiting {delay}s before retry...")
            time.sleep(delay)
            delay *= 2  # Exponential backoff
        else:
            response.raise_for_status()

    raise Exception(f"Failed after {max_retries} retries")
```

### 批处理最佳实践

1. **先使用Search API**获取ID列表，然后获取数据
2. **缓存结果**以避免冗余查询
3. **分块处理**而不是一次全部处理
4. **在请求之间添加延迟**以遵守速率限制
5. **使用 GraphQL** 进行复杂查询以最大程度地减少请求

```python
import time
from rcsbapi.search import TextQuery
from rcsbapi.data import fetch, Schema

def batch_fetch_structures(query, delay=0.5):
    """
    Fetch structures matching a query with rate limiting.

    Args:
        query: Search query object
        delay: Delay between requests in seconds
    """
    # Get list of IDs
    pdb_ids = list(query())
    print(f"Found {len(pdb_ids)} structures")

    # Fetch data for each
    results = {}
    for i, pdb_id in enumerate(pdb_ids):
        try:
            data = fetch(pdb_id, schema=Schema.ENTRY)
            results[pdb_id] = data
            print(f"Fetched {i+1}/{len(pdb_ids)}: {pdb_id}")
            time.sleep(delay)  # Rate limiting
        except Exception as e:
            print(f"Error fetching {pdb_id}: {e}")

    return results
```

## 高级用例

### 寻找药物靶点复合物

```python
from rcsbapi.search import AttributeQuery
from rcsbapi.search.attrs import rcsb_polymer_entity, rcsb_nonpolymer_entity_instance_container_identifiers

# Find structures with specific drug molecule
query = AttributeQuery(
    attribute=rcsb_nonpolymer_entity_instance_container_identifiers.comp_id,
    operator="exact_match",
    value="ATP"  # or other ligand code
)

results = list(query())
print(f"Found {len(results)} structures with ATP")
```

### 按分辨率和 R 因子过滤

```python
from rcsbapi.search import AttributeQuery
from rcsbapi.search.attrs import rcsb_entry_info, refine

# High-quality X-ray structures
resolution_query = AttributeQuery(
    attribute=rcsb_entry_info.resolution_combined,
    operator="less",
    value=2.0
)

rfactor_query = AttributeQuery(
    attribute=refine.ls_R_factor_R_free,
    operator="less",
    value=0.25
)

high_quality = resolution_query & rfactor_query
results = list(high_quality())
```

### 查找最近的结构

```python
from rcsbapi.search import AttributeQuery
from rcsbapi.search.attrs import rcsb_accession_info

# Structures released in last month
import datetime

one_month_ago = (datetime.date.today() - datetime.timedelta(days=30)).isoformat()
today = datetime.date.today().isoformat()

query = AttributeQuery(
    attribute=rcsb_accession_info.initial_release_date,
    operator="range",
    value=(one_month_ago, today)
)

recent_structures = list(query())
```

## 故障排除

### 常见错误

**404 未找到：**
- PDB ID不存在或已过时
- 检查ID是否正确（区分大小写）
- 验证条目尚未被取代

**429 请求过多：**
- 超出速率限制
- 实施指数退避
- 减少请求频率

**500 内部服务器错误：**
- 临时服务器问题
- 短暂延迟后重试
- 检查RCSB PDB状态页面

**空结果：**
- 查询限制过多
- 检查属性名称和运算符
- 验证搜索字段的数据是否存在

### 调试技巧

```python
# Enable verbose output for searches
from rcsbapi.search import TextQuery

query = TextQuery("hemoglobin")
print(query.to_dict())  # See query structure

# Check query JSON
import json
print(json.dumps(query.to_dict(), indent=2))

# Test with curl
import subprocess
result = subprocess.run(
    ["curl", "https://data.rcsb.org/rest/v1/core/entry/4HHB"],
    capture_output=True,
    text=True
)
print(result.stdout)
```

## 其他资源

- **API 文档：** https://www.rcsb.org/docs/programmatic-access/web-apis-overview
- **数据 API 重述：** https://data.rcsb.org/redoc/index.html
- **GraphQL 架构：** https://data.rcsb.org/graphql
- **Python 包文档：** https://rcsbapi.readthedocs.io/
- **GitHub 问题：** https://github.com/rcsb/py-rcsb-api/issues
- **社区论坛：** https://www.rcsb.org/help