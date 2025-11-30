<!-- 此文件由机器翻译自 SKILL.md -->

---
名称： alphafold 数据库
描述：“访问 AlphaFold 超过 2 亿个人工智能预测的蛋白质结构。通过 UniProt ID 检索结构，下载 PDB/mmCIF 文件，分析置信度指标（pLDDT、PAE），用于药物发现和结构生物学。”
---

# AlphaFold 数据库

## 概述

AlphaFold DB 是一个公共存储库，包含超过 2 亿个蛋白质的 AI 预测 3D 蛋白质结构，由 DeepMind 和 EMBL-EBI 维护。使用置信度指标访问结构预测、下载坐标文件、检索批量数据集并将预测集成到计算工作流程中。

## 何时使用此技能

在以下场景中处理人工智能预测的蛋白质结构时，应使用此技能：

- 通过 UniProt ID 或蛋白质名称检索蛋白质结构预测
- 下载PDB/mmCIF坐标文件进行结构分析
- 分析预测置信度指标（pLDDT、PAE）以评估可靠性
- 通过 Google Cloud Platform 访问批量蛋白质组数据集
- 将预测结构与实验数据进行比较
- 进行基于结构的药物发现或蛋白质工程
- 为缺乏实验结构的蛋白质建立结构模型
- 将 AlphaFold 预测集成到计算管道中

## 核心能力

### 1. 搜索和检索预测

**使用Biopython（推荐）：**

Biopython 库提供了用于检索 AlphaFold 结构的最简单的接口：

```python
from Bio.PDB import alphafold_db

# Get all predictions for a UniProt accession
predictions = list(alphafold_db.get_predictions("P00520"))

# Download structure file (mmCIF format)
for prediction in predictions:
    cif_file = alphafold_db.download_cif_for(prediction, directory="./structures")
    print(f"Downloaded: {cif_file}")

# Get Structure objects directly
from Bio.PDB import MMCIFParser
structures = list(alphafold_db.get_structural_models_for("P00520"))
```

**直接API访问：**

使用 REST 端点的查询预测：

<<<代码块_1>>>

**使用 UniProt 查找加入物：**

首先搜索 UniProt 以查找蛋白质种质：

<<<代码块_2>>>

### 2.下载结构文件

AlphaFold 为每个预测提供多种文件格式：

**可用的文件类型：**

- **模型坐标** (`model_v4.cif`)：mmCIF/PDBx 格式的原子坐标
- **置信度分数** (`confidence_v4.json`)：每个残基 pLDDT 分数 (0-100)
- **预测对齐误差** (`predicted_aligned_error_v4.json`)：残基对置信度的 PAE 矩阵

**下载地址：**

<<<代码块_3>>>

**PDB 格式（替代）：**

<<<代码块_4>>>

### 3. 使用置信度指标

AlphaFold 预测包括对解释至关重要的置信度估计：

**pLDDT（每个残基置信度）：**

<<<代码块_5>>>

**PAE（预测对齐误差）：**

PAE 表示相对域位置的置信度：

<<<代码块_6>>>

### 4. 通过 Google Cloud 进行批量数据访问

对于大规模分析，请使用 Google Cloud 数据集：

**谷歌云存储：**

```bash
# Install gsutil
uv pip install gsutil

# List available data
gsutil ls gs://public-datasets-deepmind-alphafold-v4/

# Download entire proteomes (by taxonomy ID)
gsutil -m cp gs://public-datasets-deepmind-alphafold-v4/proteomes/proteome-tax_id-9606-*.tar .

# Download specific files
gsutil cp gs://public-datasets-deepmind-alphafold-v4/accession_ids.csv .
```

**BigQuery 元数据访问：**

```python
from google.cloud import bigquery

# Initialize client
client = bigquery.Client()

# Query metadata
query = """
SELECT
  entryId,
  uniprotAccession,
  organismScientificName,
  globalMetricValue,
  fractionPlddtVeryHigh
FROM `bigquery-public-data.deepmind_alphafold.metadata`
WHERE organismScientificName = 'Homo sapiens'
  AND fractionPlddtVeryHigh > 0.8
LIMIT 100
"""

results = client.query(query).to_dataframe()
print(f"Found {len(results)} high-confidence human proteins")
```

**按物种下载：**

```python
import subprocess

def download_proteome(taxonomy_id, output_dir="./proteomes"):
    """Download all AlphaFold predictions for a species"""
    pattern = f"gs://public-datasets-deepmind-alphafold-v4/proteomes/proteome-tax_id-{taxonomy_id}-*_v4.tar"
    cmd = f"gsutil -m cp {pattern} {output_dir}/"
    subprocess.run(cmd, shell=True, check=True)

# Download E. coli proteome (tax ID: 83333)
download_proteome(83333)

# Download human proteome (tax ID: 9606)
download_proteome(9606)
```

### 5.解析和分析结构

使用 BioPython 处理下载的 AlphaFold 结构：

```python
from Bio.PDB import MMCIFParser, PDBIO
import numpy as np

# Parse mmCIF file
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("protein", "AF-P00520-F1-model_v4.cif")

# Extract coordinates
coords = []
for model in structure:
    for chain in model:
        for residue in chain:
            if 'CA' in residue:  # Alpha carbons only
                coords.append(residue['CA'].get_coord())

coords = np.array(coords)
print(f"Structure has {len(coords)} residues")

# Calculate distances
from scipy.spatial.distance import pdist, squareform
distance_matrix = squareform(pdist(coords))

# Identify contacts (< 8 Å)
contacts = np.where((distance_matrix > 0) & (distance_matrix < 8))
print(f"Number of contacts: {len(contacts[0]) // 2}")
```

**提取 B 因子（pLDDT 值）：**

AlphaFold 将 pLDDT 分数存储在 B 因子列中：

```python
from Bio.PDB import MMCIFParser

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("protein", "AF-P00520-F1-model_v4.cif")

# Extract pLDDT from B-factors
plddt_scores = []
for model in structure:
    for chain in model:
        for residue in chain:
            if 'CA' in residue:
                plddt_scores.append(residue['CA'].get_bfactor())

# Identify high-confidence regions
high_conf_regions = [(i, score) for i, score in enumerate(plddt_scores, 1) if score > 90]
print(f"High confidence residues: {len(high_conf_regions)}")
```

### 6. 批量处理多种蛋白质

高效处理多个预测：

```python
from Bio.PDB import alphafold_db
import pandas as pd

uniprot_ids = ["P00520", "P12931", "P04637"]  # Multiple proteins
results = []

for uniprot_id in uniprot_ids:
    try:
        # Get prediction
        predictions = list(alphafold_db.get_predictions(uniprot_id))

        if predictions:
            pred = predictions[0]

            # Download structure
            cif_file = alphafold_db.download_cif_for(pred, directory="./batch_structures")

            # Get confidence data
            alphafold_id = pred['entryId']
            conf_url = f"https://alphafold.ebi.ac.uk/files/{alphafold_id}-confidence_v4.json"
            conf_data = requests.get(conf_url).json()

            # Calculate statistics
            plddt_scores = conf_data['confidenceScore']
            avg_plddt = np.mean(plddt_scores)
            high_conf_fraction = sum(1 for s in plddt_scores if s > 90) / len(plddt_scores)

            results.append({
                'uniprot_id': uniprot_id,
                'alphafold_id': alphafold_id,
                'avg_plddt': avg_plddt,
                'high_conf_fraction': high_conf_fraction,
                'length': len(plddt_scores)
            })
    except Exception as e:
        print(f"Error processing {uniprot_id}: {e}")

# Create summary DataFrame
df = pd.DataFrame(results)
print(df)
```

## 安装和设置

### Python 库

```bash
# Install Biopython for structure access
uv pip install biopython

# Install requests for API access
uv pip install requests

# For visualization and analysis
uv pip install numpy matplotlib pandas scipy

# For Google Cloud access (optional)
uv pip install google-cloud-bigquery gsutil
```

### 3D 信标 API 替代方案

AlphaFold 还可以通过 3D-Beacons 联合 API 访问：

```python
import requests

# Query via 3D-Beacons
uniprot_id = "P00520"
url = f"https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/uniprot/summary/{uniprot_id}.json"
response = requests.get(url)
data = response.json()

# Filter for AlphaFold structures
af_structures = [s for s in data['structures'] if s['provider'] == 'AlphaFold DB']
```

## 常见用例

### 结构蛋白质组学
- 下载完整的蛋白质组预测以进行分析
- 识别蛋白质的高可信度结构区域
- 将预测结构与实验数据进行比较
- 建立蛋白质家族的结构模型

### 药物发现
- 检索目标蛋白质结构以进行对接研究
- 分析结合位点构象
- 识别预测结构中的可成药口袋
- 比较同源物的结构

### 蛋白质工程
- 使用 pLDDT 识别稳定/不稳定区域
- 在高置信度区域设计突变
- 使用 PAE 分析域架构
- 模型蛋白质变异和突变

### 进化研究
- 比较跨物种的直系同源结构
- 分析结构特征的保护
- 研究领域演化模式
- 识别功能上重要的区域

## 关键概念

**UniProt 登录号：** 蛋白质的主要标识符（例如“P00520”）。查询 AlphaFold DB 时需要。

**AlphaFold ID：** 内部标识符格式：`AF-[UniProt accession]-F[fragment number]`（例如“AF-P00520-F1”）。

**pLDDT（预测局部距离差异测试）：** 每个残基置信度度量 (0-100)。值越高表示预测越可信。
**PAE（预测对齐误差）：** 表示残基对之间相对位置置信度的矩阵。低值 (<5 Å) 表明有信心的相对定位。

**数据库版本：** 当前版本是 v4。文件 URL 包含版本后缀（例如 `model_v4.cif`）。

**片段数：** 大蛋白质可能会分裂成片段。片段编号出现在 AlphaFold ID 中（例如 F1、F2）。

## 置信度解释指南

**pLDDT 阈值：**
- **>90**：置信度非常高 - 适合详细分析
- **70-90**：高置信度 - 一般可靠的骨干结构
- **50-70**：低置信度 - 谨慎使用，灵活区域
- **<50**：置信度非常低 - 可能无序或不可靠

**PAE 指南：**
- **<5 Å**：域的可靠相对定位
- **5-10 Å**：对安排有中等信心
- **>15 Å**：相对位置不确定，域可能是移动的

## 资源

### 参考文献/api_reference.md

全面的 API 文档涵盖：
- 完整的REST API端点规范
- 文件格式详细信息和数据模式
- Google Cloud数据集结构和访问模式
- 高级查询示例和批处理策略
- 速率限制、缓存和最佳实践
- 解决常见问题

有关详细 API 信息、批量下载策略或处理大型数据集时，请参阅此参考。

## 重要提示

### 数据使用和归因

- AlphaFold DB 在 CC-BY-4.0 许可下免费提供
- 引用：Jumper 等人。 (2021) 自然和瓦拉迪等人。 (2022) 核酸研究
- 预测是计算模型，而不是实验结构
- 在下游分析之前始终评估置信度指标

### 版本管理

- 当前数据库版本：v4（截至 2024-2025 年）
- 文件 URL 包含版本后缀（例如 `_v4.cif`）
- 定期检查数据库更新
- 随着时间的推移，旧版本可能会被弃用

### 数据质量注意事项

- 高 pLDDT 并不能保证功能准确性
- 低置信度区域可能在体内紊乱
- PAE表示相对领域置信度，而不是绝对定位
- 预测缺乏配体、翻译后修饰和辅因子
- 不预测多链复合物（仅单链）

### 性能提示

- 使用 Biopython 进行简单的单一蛋白质访问
- 使用 Google Cloud 进行批量下载（比单个文件快得多）
- 将下载的文件缓存在本地，避免重复下载
- BigQuery 免费套餐：每月处理 1 TB 数据
- 考虑大规模下载的网络带宽

## 其他资源

- **AlphaFold DB 网站：** https://alphafold.ebi.ac.uk/
- **API 文档：** https://alphafold.ebi.ac.uk/api-docs
- **Google 云数据集：** https://cloud.google.com/blog/products/ai-machine-learning/alphafold-protein-structure-database
- **3D 信标 API：** https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/
- **AlphaFold 纸张：**
  - 自然（2021）：https://doi.org/10.1038/s41586-021-03819-2
  - 核酸研究（2024）：https://doi.org/10.1093/nar/gkad1011
- **Biopython 文档：** https://biopython.org/docs/dev/api/Bio.PDB.alphafold_db.html
- **GitHub 存储库：** https://github.com/google-deepmind/alphafold