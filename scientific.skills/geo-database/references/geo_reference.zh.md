<!-- 此文件由机器翻译自 geo_reference.md -->

# GEO 数据库参考文档

## 完整的电子公用事业 API 规范

### 概述

NCBI Entrez 编程实用程序（E-实用程序）通过一组九个服务器端程序提供对 GEO 元数据的编程访问。默认情况下，所有电子实用程序都以 XML 格式返回结果。

### 基本网址

```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
```

### 核心电子实用程序

#### eSearch - 对 ID 列表的文本查询

**用途：** 搜索数据库并返回与查询匹配的 UID 列表。

**网址模式：**
<<<代码块_1>>>

**参数：**
- `db`（必需）：要搜索的数据库（例如“gds”、“geoprofiles”）
- `term`（必需）：搜索查询字符串
- `retmax`：返回 UID 的最大数量（默认值：20，最大：10000）
- `retstart`：结果集中的起始位置（用于分页）
- `usehistory`：设置为“y”以将结果存储在历史服务器上
- `sort`：排序顺序（例如“相关性”、“pub_date”）
- `field`：将搜索限制到特定字段
- `datetype`：要限制的日期类型
- `reldate`：仅限今天 N 天内的商品
- `mindate`、`maxdate`：日期范围限制 (YYYY/MM/DD)

**示例：**
<<<代码块_2>>>

#### eSummary - 文档摘要

**用途：** 检索 UID 列表的文档摘要。

**网址模式：**
<<<代码块_3>>>

**参数：**
- `db`（必需）：数据库
- `id`（必需）：以逗号分隔的 UID 列表或 query_key+WebEnv
- `retmode`：返回格式（“xml”或“json”）
- `version`：摘要版本（推荐“2.0”）

**示例：**
<<<代码块_4>>>

#### eFetch - 全记录

**用途：** 检索 UID 列表的完整记录。

**网址模式：**
<<<代码块_5>>>

**参数：**
- `db`（必需）：数据库
- `id`（必需）：以逗号分隔的 UID 列表
- `retmode`：返回格式（“xml”，“文本”）
- `rettype`：记录类型（特定于数据库）

**示例：**
<<<代码块_6>>>

#### eLink - 跨数据库链接

**目的：**查找相同或不同数据库中的相关记录。

**网址模式：**
```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi
```

**参数：**
- `dbfrom`（必需）：源数据库
- `db`（必需）：目标数据库
- `id`（必需）：来自源数据库的 UID
- `cmd`：链接命令类型
  - “neighbor”：返回链接的 UID（默认）
  - “neighbor_score”：返回评分链接
  - “acheck”：检查链接
  - “ncheck”：链接计数
  - “llinks”：返回 LinkOut 资源的 URL

**示例：**
```python
from Bio import Entrez
Entrez.email = "your@email.com"

# Find PubMed articles linked to a GEO dataset
handle = Entrez.elink(
    dbfrom="gds",
    db="pubmed",
    id="200000001"
)
links = Entrez.read(handle)
handle.close()
```

#### ePost - 上传 UID 列表

**用途：** 将 UID 列表上传到历史服务器以供后续请求使用。

**网址模式：**
```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi
```

**参数：**
- `db`（必需）：数据库
- `id`（必需）：以逗号分隔的 UID 列表

**示例：**
```python
from Bio import Entrez
Entrez.email = "your@email.com"

# Post large list of IDs
large_id_list = [str(i) for i in range(200000001, 200000101)]
handle = Entrez.epost(db="gds", id=",".join(large_id_list))
result = Entrez.read(handle)
handle.close()

# Use returned QueryKey and WebEnv in subsequent calls
query_key = result["QueryKey"]
webenv = result["WebEnv"]
```

#### eInfo - 数据库信息

**目的：** 获取有关可用数据库及其字段的信息。

**网址模式：**
```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi
```

**参数：**
- `db`：数据库名称（省略获取所有数据库的列表）
- `version`：设置为“2.0”以获取详细字段信息

**示例：**
```python
from Bio import Entrez
Entrez.email = "your@email.com"

# Get information about gds database
handle = Entrez.einfo(db="gds", version="2.0")
info = Entrez.read(handle)
handle.close()

# Returns:
# - Database description
# - Last update date
# - Record count
# - Available search fields
# - Link information
```

### GEO 的搜索字段限定符

用于构建目标查询的常见搜索字段：

**一般领域：**
- `[Accession]`：GEO 登录号
- `[Title]`：数据集标题
- `[Author]`：作者姓名
- `[Organism]`：源生物体
- `[Entry Type]`：条目类型（例如，“按数组进行表达式分析”）
- `[Platform]`：平台加入或名称
- `[PubMed ID]`：关联的 PubMed ID

**日期字段：**
- `[Publication Date]`：发布日期（YYYY 或 YYYY/MM/DD）
- `[Submission Date]`：提交日期
- `[Modification Date]`：上次修改日期

**MeSH 术语：**
- `[MeSH Terms]`：医学主题标题
- `[MeSH Major Topic]`：主要 MeSH 主题

**研究类型字段：**
- `[DataSet Type]`：研究类型（例如“RNA-seq”、“ChIP-seq”）
- `[Sample Type]`：样本类型

**复杂查询示例：**
```python
query = """
    (breast cancer[MeSH] OR breast neoplasms[Title]) AND
    Homo sapiens[Organism] AND
    expression profiling by array[Entry Type] AND
    2020:2024[Publication Date] AND
    GPL570[Platform]
"""
```

## SOFT 文件格式规范

### 概述

SOFT（简单综合文本格式）是 GEO 的主要数据交换格式。文件的结构为带有数据表的键值对。

### 文件类型

**家庭软件文件：**
- 文件名：`GSExxxxx_family.soft.gz`
- 包含：包含所有示例和平台的完整系列
- 大小：可以非常大（压缩后 100 MB）
- 用途：完整的数据提取

**系列矩阵文件：**
- 文件名：`GSExxxxx_series_matrix.txt.gz`
- 包含：具有最少元数据的表达矩阵
- 大小：比家庭文件小
- 用途：快速访问表达数据

**平台软件文件：**
- 文件名：`GPLxxxxx.soft`
- 包含：平台注释和探测信息
- 用途：将探针映射到基因

### 软件文件结构

```
^DATABASE = GeoMiame
!Database_name = Gene Expression Omnibus (GEO)
!Database_institute = NCBI NLM NIH
!Database_web_link = http://www.ncbi.nlm.nih.gov/geo
!Database_email = geo@ncbi.nlm.nih.gov

^SERIES = GSExxxxx
!Series_title = Study Title Here
!Series_summary = Study description and background...
!Series_overall_design = Experimental design...
!Series_type = Expression profiling by array
!Series_pubmed_id = 12345678
!Series_submission_date = Jan 01 2024
!Series_last_update_date = Jan 15 2024
!Series_contributor = John,Doe
!Series_contributor = Jane,Smith
!Series_sample_id = GSMxxxxxx
!Series_sample_id = GSMxxxxxx

^PLATFORM = GPLxxxxx
!Platform_title = Platform Name
!Platform_distribution = commercial or custom
!Platform_organism = Homo sapiens
!Platform_manufacturer = Affymetrix
!Platform_technology = in situ oligonucleotide
!Platform_data_row_count = 54675
#ID = Probe ID
#GB_ACC = GenBank accession
#SPOT_ID = Spot identifier
#Gene Symbol = Gene symbol
#Gene Title = Gene title
!platform_table_begin
ID    GB_ACC    SPOT_ID    Gene Symbol    Gene Title
1007_s_at    U48705    -    DDR1    discoidin domain receptor...
1053_at    M87338    -    RFC2    replication factor C...
!platform_table_end

^SAMPLE = GSMxxxxxx
!Sample_title = Sample name
!Sample_source_name_ch1 = cell line XYZ
!Sample_organism_ch1 = Homo sapiens
!Sample_characteristics_ch1 = cell type: epithelial
!Sample_characteristics_ch1 = treatment: control
!Sample_molecule_ch1 = total RNA
!Sample_label_ch1 = biotin
!Sample_platform_id = GPLxxxxx
!Sample_data_processing = normalization method
#ID_REF = Probe identifier
#VALUE = Expression value
!sample_table_begin
ID_REF    VALUE
1007_s_at    8.456
1053_at    7.234
!sample_table_end
```

### 解析 SOFT 文件

**使用 GEOparse：**
```python
import GEOparse

# Parse series
gse = GEOparse.get_GEO(filepath="GSE123456_family.soft.gz")

# Access metadata
metadata = gse.metadata
phenotype_data = gse.phenotype_data

# Access samples
for gsm_name, gsm in gse.gsms.items():
    sample_data = gsm.table
    sample_metadata = gsm.metadata

# Access platforms
for gpl_name, gpl in gse.gpls.items():
    platform_table = gpl.table
    platform_metadata = gpl.metadata
```

**手动解析：**
```python
import gzip

def parse_soft_file(filename):
    """Basic SOFT file parser"""
    sections = {}
    current_section = None
    current_metadata = {}
    current_table = []
    in_table = False

    with gzip.open(filename, 'rt') as f:
        for line in f:
            line = line.strip()

            # New section
            if line.startswith('^'):
                if current_section:
                    sections[current_section] = {
                        'metadata': current_metadata,
                        'table': current_table
                    }
                parts = line[1:].split(' = ')
                current_section = parts[1] if len(parts) > 1 else parts[0]
                current_metadata = {}
                current_table = []
                in_table = False

            # Metadata
            elif line.startswith('!'):
                if in_table:
                    in_table = False
                key_value = line[1:].split(' = ', 1)
                if len(key_value) == 2:
                    key, value = key_value
                    if key in current_metadata:
                        if isinstance(current_metadata[key], list):
                            current_metadata[key].append(value)
                        else:
                            current_metadata[key] = [current_metadata[key], value]
                    else:
                        current_metadata[key] = value

            # Table data
            elif line.startswith('#') or in_table:
                in_table = True
                current_table.append(line)

    return sections
```

## MINIML 文件格式

### 概述

MINiML（标记语言中的MIAME 表示法）是GEO 基于XML 的数据交换格式。

### 文件结构

```xml
<?xml version="1.0" encoding="UTF-8"?>
<MINiML xmlns="http://www.ncbi.nlm.nih.gov/geo/info/MINiML"
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <Series iid="GDS123">
    <Status>
      <Submission-Date>2024-01-01</Submission-Date>
      <Release-Date>2024-01-15</Release-Date>
      <Last-Update-Date>2024-01-15</Last-Update-Date>
    </Status>
    <Title>Study Title</Title>
    <Summary>Study description...</Summary>
    <Overall-Design>Experimental design...</Overall-Design>
    <Type>Expression profiling by array</Type>
    <Contributor>
      <Person>
        <First>John</First>
        <Last>Doe</Last>
      </Person>
    </Contributor>
  </Series>

  <Platform iid="GPL123">
    <Title>Platform Name</Title>
    <Distribution>commercial</Distribution>
    <Technology>in situ oligonucleotide</Technology>
    <Organism taxid="9606">Homo sapiens</Organism>
    <Data-Table>
      <Column position="1">
        <Name>ID</Name>
        <Description>Probe identifier</Description>
      </Column>
      <Data>
        <Row>
          <Cell column="1">1007_s_at</Cell>
          <Cell column="2">U48705</Cell>
        </Row>
      </Data>
    </Data-Table>
  </Platform>

  <Sample iid="GSM123">
    <Title>Sample name</Title>
    <Source>cell line XYZ</Source>
    <Organism taxid="9606">Homo sapiens</Organism>
    <Characteristics tag="cell type">epithelial</Characteristics>
    <Characteristics tag="treatment">control</Characteristics>
    <Platform-Ref ref="GPL123"/>
    <Data-Table>
      <Column position="1">
        <Name>ID_REF</Name>
      </Column>
      <Column position="2">
        <Name>VALUE</Name>
      </Column>
      <Data>
        <Row>
          <Cell column="1">1007_s_at</Cell>
          <Cell column="2">8.456</Cell>
        </Row>
      </Data>
    </Data-Table>
  </Sample>
</MINiML>
```

## FTP 目录结构

### 系列文件

**图案：**
```
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE{nnn}nnn/GSE{xxxxx}/
```

其中 `{nnn}` 表示用“nnn”替换最后 3 位数字，`{xxxxx}` 是完整加入。

**示例：**
- GSE123456 → `/geo/series/GSE123nnn/GSE123456/`
- GSE1234 → `/geo/series/GSE1nnn/GSE1234/`
- GSE100001 → `/geo/series/GSE100nnn/GSE100001/`

**子目录：**
- `/matrix/` - 系列矩阵文件
- `/soft/` - 系列软件文件
- `/miniml/` - MINIML XML 文件
- `/suppl/` - 补充文件

**文件类型：**
```
matrix/
  └── GSE123456_series_matrix.txt.gz

soft/
  └── GSE123456_family.soft.gz

miniml/
  └── GSE123456_family.xml.tgz

suppl/
  ├── GSE123456_RAW.tar
  ├── filelist.txt
  └── [various supplementary files]
```

### 示例文件

**图案：**
```
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM{nnn}nnn/GSM{xxxxx}/
```

**子目录：**
- `/suppl/` - 示例特定的补充文件

### 平台文件

**图案：**
```
ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL{nnn}nnn/GPL{xxxxx}/
```

**文件类型：**
```
soft/
  └── GPL570.soft.gz

miniml/
  └── GPL570.xml

annot/
  └── GPL570.annot.gz  # Enhanced annotation (if available)
```

## 高级 GEOparse 用法

### 自定义解析选项

```python
import GEOparse

# Parse with custom options
gse = GEOparse.get_GEO(
    geo="GSE123456",
    destdir="./data",
    silent=False,  # Show progress
    how="full",  # Parse mode: "full", "quick", "brief"
    annotate_gpl=True,  # Include platform annotation
    geotype="GSE"  # Explicit type
)

# Access specific sample
gsm = gse.gsms['GSM1234567']

# Get expression values for specific probe
probe_id = "1007_s_at"
if hasattr(gsm, 'table'):
    probe_data = gsm.table[gsm.table['ID_REF'] == probe_id]

# Get all characteristics
characteristics = {}
for key, values in gsm.metadata.items():
    if key.startswith('characteristics'):
        for value in (values if isinstance(values, list) else [values]):
            if ':' in value:
                char_key, char_value = value.split(':', 1)
                characteristics[char_key.strip()] = char_value.strip()
```

### 使用平台注释

```python
import GEOparse
import pandas as pd

gse = GEOparse.get_GEO(geo="GSE123456", destdir="./data")

# Get platform
gpl = list(gse.gpls.values())[0]

# Extract annotation table
if hasattr(gpl, 'table'):
    annotation = gpl.table

    # Common annotation columns:
    # - ID: Probe identifier
    # - Gene Symbol: Gene symbol
    # - Gene Title: Gene description
    # - GB_ACC: GenBank accession
    # - Gene ID: Entrez Gene ID
    # - RefSeq: RefSeq accession
    # - UniGene: UniGene cluster

    # Map probes to genes
    probe_to_gene = dict(zip(
        annotation['ID'],
        annotation['Gene Symbol']
    ))

    # Handle multiple probes per gene
    gene_to_probes = {}
    for probe, gene in probe_to_gene.items():
        if gene and gene != '---':
            if gene not in gene_to_probes:
                gene_to_probes[gene] = []
            gene_to_probes[gene].append(probe)
```

### 处理大型数据集

```python
import GEOparse
import pandas as pd
import numpy as np

def process_large_gse(gse_id, chunk_size=1000):
    """Process large GEO series in chunks"""
    gse = GEOparse.get_GEO(geo=gse_id, destdir="./data")

    # Get sample list
    sample_list = list(gse.gsms.keys())

    # Process in chunks
    for i in range(0, len(sample_list), chunk_size):
        chunk_samples = sample_list[i:i+chunk_size]

        # Extract data for chunk
        chunk_data = {}
        for gsm_id in chunk_samples:
            gsm = gse.gsms[gsm_id]
            if hasattr(gsm, 'table'):
                chunk_data[gsm_id] = gsm.table['VALUE']

        # Process chunk
        chunk_df = pd.DataFrame(chunk_data)

        # Save chunk results
        chunk_df.to_csv(f"chunk_{i//chunk_size}.csv")

        print(f"Processed {i+len(chunk_samples)}/{len(sample_list)} samples")
```

## 常见问题故障排除

### 问题：GEOparse 无法下载

**症状：** 超时错误、连接失败

**解决方案：**
1. 检查互联网连接
2.先尝试直接通过FTP下载
3.解析本地文件：
```python
gse = GEOparse.get_GEO(filepath="./local/GSE123456_family.soft.gz")
```
4.增加超时（如果需要，修改GEOparse源）

### 问题：缺少表达数据

**症状：** `pivot_samples()` 失败或返回空

**原因：** 并非所有系列都有系列矩阵文件（较旧的提交）

**解决方案：** 解析各个样本表：
```python
expression_data = {}
for gsm_name, gsm in gse.gsms.items():
    if hasattr(gsm, 'table') and 'VALUE' in gsm.table.columns:
        expression_data[gsm_name] = gsm.table.set_index('ID_REF')['VALUE']

expression_df = pd.DataFrame(expression_data)
```

### 问题：探针 ID 不一致

**症状：** 样品之间的探针 ID 不匹配

**原因：** 平台版本不同或样本处理不同

**解决方案：** 使用平台注释进行标准化：
```python
# Get common probe set
all_probes = set()
for gsm in gse.gsms.values():
    if hasattr(gsm, 'table'):
        all_probes.update(gsm.table['ID_REF'].values)

# Create standardized matrix
standardized_data = {}
for gsm_name, gsm in gse.gsms.items():
    if hasattr(gsm, 'table'):
        sample_data = gsm.table.set_index('ID_REF')['VALUE']
        standardized_data[gsm_name] = sample_data.reindex(all_probes)

expression_df = pd.DataFrame(standardized_data)
```

### 问题：电子公用事业费率限制

**症状：** HTTP 429 错误、响应缓慢

**解决方案：**
1. 从 NCBI 获取 API 密钥
2. 实施速率限制：
```python
import time
from functools import wraps

def rate_limit(calls_per_second=3):
    min_interval = 1.0 / calls_per_second

    def decorator(func):
        last_called = [0.0]

        @wraps(func)
        def wrapper(*args, **kwargs):
            elapsed = time.time() - last_called[0]
            wait_time = min_interval - elapsed
            if wait_time > 0:
                time.sleep(wait_time)
            result = func(*args, **kwargs)
            last_called[0] = time.time()
            return result
        return wrapper
    return decorator

@rate_limit(calls_per_second=3)
def safe_esearch(query):
    handle = Entrez.esearch(db="gds", term=query)
    results = Entrez.read(handle)
    handle.close()
    return results
```

### 问题：大型数据集的内存错误

**症状：** 内存错误、系统速度变慢

**解决方案：**
1. 分块处理数据
2. 对表达数据使用稀疏矩阵
3. 仅加载必要的列
4. 使用内存高效的数据类型：
```python
import pandas as pd

# Read with specific dtypes
expression_df = pd.read_csv(
    "expression_matrix.csv",
    dtype={'ID': str, 'GSM1': np.float32}  # Use float32 instead of float64
)

# Or use sparse format for mostly-zero data
import scipy.sparse as sp
sparse_matrix = sp.csr_matrix(expression_df.values)
```

## 特定于平台的注意事项

### Affymetrix 阵列

- 探针 ID 格式：`1007_s_at`、`1053_at`
- 每个基因共有多个探针组
- 检查 `_at`、`_s_at`、`_x_at` 后缀
- 可能需要 RMA 或 MAS5 标准化

### Illumina 阵列

- 探针 ID 格式：`ILMN_1234567`
- 留意重复的探针
- 可能需要 BeadChip 特定处理

### RNA测序

- 可能没有传统的“探针”
- 检查基因 ID（Ensembl、Entrez）
- 计数与 FPKM/TPM 值
- 可能需要单独的计数文件

### 两通道阵列

- 在元数据中查找 `_ch1` 和 `_ch2` 后缀
- VALUE_ch1、VALUE_ch2 列
- 可能需要比率或强度值
- 检查染料交换实验

## 最佳实践总结

1. **在使用电子实用程序之前始终设置 Entrez.email**
2. **使用API密钥**以获得更好的速率限制
3. **在本地缓存下载的文件**
4. **分析前检查数据质量**
5. **验证平台注释**是最新的
6. **文档数据处理**步骤
7. **使用数据时引用原始研究**
8. **在荟萃分析中检查批次效应**
9. **使用独立数据集验证结果**
10. **遵循 NCBI 使用指南**