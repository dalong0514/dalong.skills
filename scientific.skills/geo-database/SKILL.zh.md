<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：地理数据库
描述：“访问 NCBI GEO 获取基因表达/基因组数据。搜索/下载微阵列和 RNA-seq 数据集（GSE、GSM、GPL），检索 SOFT/Matrix 文件，用于转录组学和表达分析。”
---

# 地理数据库

## 概述

Gene Expression Omnibus (GEO) 是 NCBI 的高通量基因表达和功能基因组学数据的公共存储库。 GEO 包含超过 264,000 项研究，其中来自基于阵列和基于序列的实验的超过 800 万个样本。

## 何时使用此技能

在搜索基因表达数据集、检索实验数据、下载原始文件和处理后的文件、查询表达谱或将 GEO 数据集成到计算分析工作流程中时，应使用此技能。

## 核心能力

### 1. 了解 GEO 数据组织

GEO 使用不同的加入类型分层组织数据：

**系列（GSE）：** 使用一组相关样本进行的完整实验
- 示例：GSE123456
- 包含实验设计、样本和整体研究信息
- GEO 最大的组织单位
- 当前数量：264,928+ 系列

**样品 (GSM)：** 单个实验样品或生物复制品
- 示例：GSM987654
- 包含单独的样本数据、协议和元数据
- 链接到平台和系列
- 当前计数：8,068,632+ 样本

**平台 (GPL)：** 使用的微阵列或测序平台
- 示例：GPL570（Affymetrix 人类基因组 U133 Plus 2.0 阵列）
- 描述技术和探针/功能注释
- 在多个实验中共享
- 当前数量：27,739+ 个平台

**数据集 (GDS)：** 具有一致格式的精选集合
- 示例：GDS5678
- 按研究设计组织的实验可比样本
- 进行差异分析处理
- GEO 数据子集（4,348 个精选数据集）
- 快速比较分析的理想选择

**概况：**与序列特征相关的基因特异性表达数据
- 可通过基因名称或注释查询
- Entrez 基因交叉引用
- 在所有研究中实现以基因为中心的搜索

### 2. 搜索GEO数据

**GEO 数据集搜索：**

按关键词、生物体或实验条件搜索研究：

```python
from Bio import Entrez

# Configure Entrez (required)
Entrez.email = "your.email@example.com"

# Search for datasets
def search_geo_datasets(query, retmax=20):
    """Search GEO DataSets database"""
    handle = Entrez.esearch(
        db="gds",
        term=query,
        retmax=retmax,
        usehistory="y"
    )
    results = Entrez.read(handle)
    handle.close()
    return results

# Example searches
results = search_geo_datasets("breast cancer[MeSH] AND Homo sapiens[Organism]")
print(f"Found {results['Count']} datasets")

# Search by specific platform
results = search_geo_datasets("GPL570[Accession]")

# Search by study type
results = search_geo_datasets("expression profiling by array[DataSet Type]")
```

**GEO 档案搜索：**

查找基因特异性表达模式：

<<<代码块_1>>>

**高级搜索模式：**

<<<代码块_2>>>

### 3. 使用 GEOparse 检索 GEO 数据（推荐）

**GEOparse** 是用于访问 GEO 数据的主要 Python 库：

**安装：**
<<<代码块_3>>>

**基本用法：**

<<<代码块_4>>>

**使用表达数据：**

<<<代码块_5>>>

**访问补充文件：**

<<<代码块_6>>>

**过滤和子集数据：**

```python
import GEOparse

gse = GEOparse.get_GEO(geo="GSE123456", destdir="./data")

# Filter samples by metadata
control_samples = [
    gsm_name for gsm_name, gsm in gse.gsms.items()
    if 'control' in gsm.metadata.get('title', [''])[0].lower()
]

treatment_samples = [
    gsm_name for gsm_name, gsm in gse.gsms.items()
    if 'treatment' in gsm.metadata.get('title', [''])[0].lower()
]

print(f"Control samples: {len(control_samples)}")
print(f"Treatment samples: {len(treatment_samples)}")

# Extract subset expression matrix
expression_df = gse.pivot_samples('VALUE')
control_expr = expression_df[control_samples]
treatment_expr = expression_df[treatment_samples]
```

### 4. 使用 NCBI 电子实用程序进行 GEO 访问

**电子实用程序**提供对 GEO 元数据的较低级别的编程访问：

**基本电子公用事业工作流程：**

```python
from Bio import Entrez
import time

Entrez.email = "your.email@example.com"

# Step 1: Search for GEO entries
def search_geo(query, db="gds", retmax=100):
    """Search GEO using E-utilities"""
    handle = Entrez.esearch(
        db=db,
        term=query,
        retmax=retmax,
        usehistory="y"
    )
    results = Entrez.read(handle)
    handle.close()
    return results

# Step 2: Fetch summaries
def fetch_geo_summaries(id_list, db="gds"):
    """Fetch document summaries for GEO entries"""
    ids = ",".join(id_list)
    handle = Entrez.esummary(db=db, id=ids)
    summaries = Entrez.read(handle)
    handle.close()
    return summaries

# Step 3: Fetch full records
def fetch_geo_records(id_list, db="gds"):
    """Fetch full GEO records"""
    ids = ",".join(id_list)
    handle = Entrez.efetch(db=db, id=ids, retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    return records

# Example workflow
search_results = search_geo("breast cancer AND Homo sapiens")
id_list = search_results['IdList'][:5]

summaries = fetch_geo_summaries(id_list)
for summary in summaries:
    print(f"GDS: {summary.get('Accession', 'N/A')}")
    print(f"Title: {summary.get('title', 'N/A')}")
    print(f"Samples: {summary.get('n_samples', 'N/A')}")
    print()
```

**使用电子实用程序进行批处理：**

```python
from Bio import Entrez
import time

Entrez.email = "your.email@example.com"

def batch_fetch_geo_metadata(accessions, batch_size=100):
    """Fetch metadata for multiple GEO accessions"""
    results = {}

    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i + batch_size]

        # Search for each accession
        for accession in batch:
            try:
                query = f"{accession}[Accession]"
                search_handle = Entrez.esearch(db="gds", term=query)
                search_results = Entrez.read(search_handle)
                search_handle.close()

                if search_results['IdList']:
                    # Fetch summary
                    summary_handle = Entrez.esummary(
                        db="gds",
                        id=search_results['IdList'][0]
                    )
                    summary = Entrez.read(summary_handle)
                    summary_handle.close()
                    results[accession] = summary[0]

                # Be polite to NCBI servers
                time.sleep(0.34)  # Max 3 requests per second

            except Exception as e:
                print(f"Error fetching {accession}: {e}")

    return results

# Fetch metadata for multiple datasets
gse_list = ["GSE100001", "GSE100002", "GSE100003"]
metadata = batch_fetch_geo_metadata(gse_list)
```

### 5. 直接通过 FTP 访问数据文件

**GEO 数据的 FTP URL：**

GEO数据可以直接通过FTP下载：

```python
import ftplib
import os

def download_geo_ftp(accession, file_type="matrix", dest_dir="./data"):
    """Download GEO files via FTP"""
    # Construct FTP path based on accession type
    if accession.startswith("GSE"):
        # Series files
        gse_num = accession[3:]
        base_num = gse_num[:-3] + "nnn"
        ftp_path = f"/geo/series/GSE{base_num}/{accession}/"

        if file_type == "matrix":
            filename = f"{accession}_series_matrix.txt.gz"
        elif file_type == "soft":
            filename = f"{accession}_family.soft.gz"
        elif file_type == "miniml":
            filename = f"{accession}_family.xml.tgz"

    # Connect to FTP server
    ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
    ftp.login()
    ftp.cwd(ftp_path)

    # Download file
    os.makedirs(dest_dir, exist_ok=True)
    local_file = os.path.join(dest_dir, filename)

    with open(local_file, 'wb') as f:
        ftp.retrbinary(f'RETR {filename}', f.write)

    ftp.quit()
    print(f"Downloaded: {local_file}")
    return local_file

# Download series matrix file
download_geo_ftp("GSE123456", file_type="matrix")

# Download SOFT format file
download_geo_ftp("GSE123456", file_type="soft")
```

**使用 wget 或 curl 进行下载：**

```bash
# Download series matrix file
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/matrix/GSE123456_series_matrix.txt.gz

# Download all supplementary files for a series
wget -r -np -nd ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/suppl/

# Download SOFT format family file
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/soft/GSE123456_family.soft.gz
```

### 6. 分析 GEO 数据

**质量控制和预处理：**

```python
import GEOparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load dataset
gse = GEOparse.get_GEO(geo="GSE123456", destdir="./data")
expression_df = gse.pivot_samples('VALUE')

# Check for missing values
print(f"Missing values: {expression_df.isnull().sum().sum()}")

# Log transformation (if needed)
if expression_df.min().min() > 0:  # Check if already log-transformed
    if expression_df.max().max() > 100:
        expression_df = np.log2(expression_df + 1)
        print("Applied log2 transformation")

# Distribution plots
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
expression_df.plot.box(ax=plt.gca())
plt.title("Expression Distribution per Sample")
plt.xticks(rotation=90)

plt.subplot(1, 2, 2)
expression_df.mean(axis=1).hist(bins=50)
plt.title("Gene Expression Distribution")
plt.xlabel("Average Expression")

plt.tight_layout()
plt.savefig("geo_qc.png", dpi=300, bbox_inches='tight')
```

**差异表达分析：**

```python
import GEOparse
import pandas as pd
import numpy as np
from scipy import stats

gse = GEOparse.get_GEO(geo="GSE123456", destdir="./data")
expression_df = gse.pivot_samples('VALUE')

# Define sample groups
control_samples = ["GSM1", "GSM2", "GSM3"]
treatment_samples = ["GSM4", "GSM5", "GSM6"]

# Calculate fold changes and p-values
results = []
for gene in expression_df.index:
    control_expr = expression_df.loc[gene, control_samples]
    treatment_expr = expression_df.loc[gene, treatment_samples]

    # Calculate statistics
    fold_change = treatment_expr.mean() - control_expr.mean()
    t_stat, p_value = stats.ttest_ind(treatment_expr, control_expr)

    results.append({
        'gene': gene,
        'log2_fold_change': fold_change,
        'p_value': p_value,
        'control_mean': control_expr.mean(),
        'treatment_mean': treatment_expr.mean()
    })

# Create results DataFrame
de_results = pd.DataFrame(results)

# Multiple testing correction (Benjamini-Hochberg)
from statsmodels.stats.multitest import multipletests
_, de_results['q_value'], _, _ = multipletests(
    de_results['p_value'],
    method='fdr_bh'
)

# Filter significant genes
significant_genes = de_results[
    (de_results['q_value'] < 0.05) &
    (abs(de_results['log2_fold_change']) > 1)
]

print(f"Significant genes: {len(significant_genes)}")
significant_genes.to_csv("de_results.csv", index=False)
```

**相关性和聚类分析：**

```python
import GEOparse
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist

gse = GEOparse.get_GEO(geo="GSE123456", destdir="./data")
expression_df = gse.pivot_samples('VALUE')

# Sample correlation heatmap
sample_corr = expression_df.corr()

plt.figure(figsize=(10, 8))
sns.heatmap(sample_corr, cmap='coolwarm', center=0,
            square=True, linewidths=0.5)
plt.title("Sample Correlation Matrix")
plt.tight_layout()
plt.savefig("sample_correlation.png", dpi=300, bbox_inches='tight')

# Hierarchical clustering
distances = pdist(expression_df.T, metric='correlation')
linkage = hierarchy.linkage(distances, method='average')

plt.figure(figsize=(12, 6))
hierarchy.dendrogram(linkage, labels=expression_df.columns)
plt.title("Hierarchical Clustering of Samples")
plt.xlabel("Samples")
plt.ylabel("Distance")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("sample_clustering.png", dpi=300, bbox_inches='tight')
```

### 7. 批处理多个数据集

**下载并处理多个系列：**

```python
import GEOparse
import pandas as pd
import os

def batch_download_geo(gse_list, destdir="./geo_data"):
    """Download multiple GEO series"""
    results = {}

    for gse_id in gse_list:
        try:
            print(f"Processing {gse_id}...")
            gse = GEOparse.get_GEO(geo=gse_id, destdir=destdir)

            # Extract key information
            results[gse_id] = {
                'title': gse.metadata.get('title', ['N/A'])[0],
                'organism': gse.metadata.get('organism', ['N/A'])[0],
                'platform': list(gse.gpls.keys())[0] if gse.gpls else 'N/A',
                'num_samples': len(gse.gsms),
                'submission_date': gse.metadata.get('submission_date', ['N/A'])[0]
            }

            # Save expression data
            if hasattr(gse, 'pivot_samples'):
                expr_df = gse.pivot_samples('VALUE')
                expr_df.to_csv(f"{destdir}/{gse_id}_expression.csv")
                results[gse_id]['num_genes'] = len(expr_df)

        except Exception as e:
            print(f"Error processing {gse_id}: {e}")
            results[gse_id] = {'error': str(e)}

    # Save summary
    summary_df = pd.DataFrame(results).T
    summary_df.to_csv(f"{destdir}/batch_summary.csv")

    return results

# Process multiple datasets
gse_list = ["GSE100001", "GSE100002", "GSE100003"]
results = batch_download_geo(gse_list)
```

**跨研究的荟萃分析：**

```python
import GEOparse
import pandas as pd
import numpy as np

def meta_analysis_geo(gse_list, gene_of_interest):
    """Perform meta-analysis of gene expression across studies"""
    results = []

    for gse_id in gse_list:
        try:
            gse = GEOparse.get_GEO(geo=gse_id, destdir="./data")

            # Get platform annotation
            gpl = list(gse.gpls.values())[0]

            # Find gene in platform
            if hasattr(gpl, 'table'):
                gene_probes = gpl.table[
                    gpl.table['Gene Symbol'].str.contains(
                        gene_of_interest,
                        case=False,
                        na=False
                    )
                ]

                if not gene_probes.empty:
                    expr_df = gse.pivot_samples('VALUE')

                    for probe_id in gene_probes['ID']:
                        if probe_id in expr_df.index:
                            expr_values = expr_df.loc[probe_id]

                            results.append({
                                'study': gse_id,
                                'probe': probe_id,
                                'mean_expression': expr_values.mean(),
                                'std_expression': expr_values.std(),
                                'num_samples': len(expr_values)
                            })

        except Exception as e:
            print(f"Error in {gse_id}: {e}")

    return pd.DataFrame(results)

# Meta-analysis for TP53
gse_studies = ["GSE100001", "GSE100002", "GSE100003"]
meta_results = meta_analysis_geo(gse_studies, "TP53")
print(meta_results)
```

## 安装和设置

### Python 库

```bash
# Primary GEO access library (recommended)
uv pip install GEOparse

# For E-utilities and programmatic NCBI access
uv pip install biopython

# For data analysis
uv pip install pandas numpy scipy

# For visualization
uv pip install matplotlib seaborn

# For statistical analysis
uv pip install statsmodels scikit-learn
```

### 配置

设置 NCBI 电子实用程序访问：

```python
from Bio import Entrez

# Always set your email (required by NCBI)
Entrez.email = "your.email@example.com"

# Optional: Set API key for increased rate limits
# Get your API key from: https://www.ncbi.nlm.nih.gov/account/
Entrez.api_key = "your_api_key_here"

# With API key: 10 requests/second
# Without API key: 3 requests/second
```

## 常见用例

### 转录组学研究
- 下载特定条件下的基因表达数据
- 比较不同研究的表达谱
- 识别差异表达基因
- 跨多个数据集执行荟萃分析

### 药物反应研究
- 分析药物治疗后基因表达变化
- 识别药物反应的生物标志物
- 比较不同细胞系或患者的药物效果
- 建立药物敏感性预测模型

### 疾病生物学
- 研究疾病与正常组织中的基因表达
- 识别疾病相关的表达特征
- 比较患者亚组和疾病阶段
- 将表达与临床结果相关联

### 生物标志物发现
- 筛查诊断或预后标志物
- 验证独立队列中的生物标志物
- 比较跨平台的标记性能
- 将表达与临床数据相结合

## 关键概念

**SOFT（文本中的简单综合格式）：** GEO 基于文本的主要格式，包含元数据和数据表。通过 GEOparse 轻松解析。

**MINiML（标记语言中的 MIAME 表示法）：** GEO 数据的 XML 格式，用于编程访问和数据交换。

**系列矩阵：** 制表符分隔的表达矩阵，其中样本为列，基因/探针为行。获取表达数据的最快格式。

**MIAME 合规性：** 有关微阵列实验的最低信息 - GEO 对所有提交强制执行的标准化注释。

**表达值类型：**不同类型的表达测量（原始信号、标准化、对数转换）。经常检查平台和处理方法。

**平台注释：** 将探针/特征 ID 映射到基因。对于表达数据的生物学解释至关重要。

## GEO2R 网络工具

要快速分析而无需编码，请使用 GEO2R：

- 集成到 GEO 中的基于网络的统计分析工具
- 可访问：https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSExxxxx
- 进行差异表达分析
- 生成 R 脚本以实现可重复性
- 对于下载数据之前的探索性分析很有用

## 速率限制和最佳实践

**NCBI 电子公用事业费率限制：**
- 没有 API 密钥：每秒 3 个请求
- 使用 API 密钥：每秒 10 个请求
- 在请求之间实现延迟：`time.sleep(0.34)`（无 API 密钥）或 `time.sleep(0.1)`（有 API 密钥）

**FTP 访问：**
- FTP 下载没有速率限制
- 批量下载的首选方法
- 可以使用 wget -r 下载整个目录

**GEOparse 缓存：**
- GEOparse 自动将下载的文件缓存在 destdir 中
- 后续调用使用缓存数据
- 定期清理缓存以节省磁盘空间

**最佳实践：**
- 使用 GEOparse 进行系列级访问（最简单）
- 使用电子实用程序进行元数据搜索和批量查询
- 使用 FTP 进行直接文件下载和批量操作
- 本地缓存数据，避免重复下载
- 使用 Biopython 时始终设置 Entrez.email

## 资源

### 参考文献/geo_reference.md

综合参考文档涵盖：
- 详细的电子实用程序 API 规范和端点
- 完整的 SOFT 和 MINIML 文件格式文档
- 高级 GEOparse 使用模式和示例
- FTP目录结构和文件命名约定
- 数据处理管道和标准化方法
- 常见问题故障排除和错误处理
- 特定于平台的注意事项和怪癖

有关深入的技术细节、复杂的查询模式或使用不常见的数据格式时，请参阅此参考。

## 重要提示

### 数据质量注意事项

- GEO接受用户提交的具有不同质量标准的数据
- 经常检查平台注释和处理方法
- 验证样本元数据和实验设计
- 谨慎对待研究中的批次效应
- 考虑重新处理原始数据以保持一致性

### 文件大小警告

- 系列矩阵文件可能很大（大型研究>1 GB）
- 补充文件（例如 CEL 文件）可能非常大
- 下载前规划足够的磁盘空间
- 考虑增量下载示例

### 数据使用和引用

- GEO数据可免费供研究使用
- 使用 GEO 数据时始终引用原始研究
- 引用 GEO 数据库：Barrett 等人。 (2013) 核酸研究
- 检查单个数据集使用限制（如果有）
- 遵循 NCBI 指南进行编程访问

### 常见陷阱

- 不同平台使用不同的探针ID（需要注释映射）
- 表达式值可以是原始的、标准化的或对数转换的（检查元数据）
- 样本元数据在不同研究中的格式可能不一致
- 并非所有系列都有系列矩阵文件（较旧的提交）
- 平台注释可能已过时（基因已重命名、ID 已弃用）

## 其他资源

- **GEO 网站：** https://www.ncbi.nlm.nih.gov/geo/
- **GEO 提交指南：** https://www.ncbi.nlm.nih.gov/geo/info/submission.html
- **GEOparse 文档：** https://geoparse.readthedocs.io/
- **电子实用程序文档：** https://www.ncbi.nlm.nih.gov/books/NBK25501/
- **GEO FTP 站点：** ftp://ftp.ncbi.nlm.nih.gov/geo/
- **GEO2R 工具：** https://www.ncbi.nlm.nih.gov/geo/geo2r/
- **NCBI API 密钥：** https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
- **Biopython 教程：** https://biopython.org/DIST/docs/tutorial/Tutorial.html