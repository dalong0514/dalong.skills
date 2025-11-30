<!-- 此文件由机器翻译自 api_reference.md -->

# AlphaFold 数据库 API 参考

本文档提供了以编程方式访问 AlphaFold 蛋白质结构数据库的全面技术文档。

## 目录

1. [REST API 端点](#rest-api-endpoints)
2. [文件访问模式](#file-access-patterns)
3. [数据模式](#data-schemas)
4. [谷歌云访问](#google-cloud-access)
5.[BigQuery 架构](#bigquery-schema)
6. [最佳实践](#best-practices)
7. [错误处理](#error-handling)
8. [速率限制](#rate-limiting)

---

## REST API 端点

### 基本网址

```
https://alphafold.ebi.ac.uk/api/
```

### 1. 通过 UniProt Accession 获取预测

**端点：** `/prediction/{uniprot_id}`

**方法：** 获取

**描述：** 检索给定 UniProt 加入的 AlphaFold 预测元数据。

**参数：**
- `uniprot_id`（必需）：UniProt 加入（例如“P00520”）

**请求示例：**
<<<代码块_1>>>

**响应示例：**
<<<代码块_2>>>

**响应字段：**
- `entryId`：AlphaFold 内部标识符（格式：AF-{uniprot}-F{fragment}）
- `gene`：基因符号
- `uniprotAccession`：UniProt 加入
- `uniprotId`：UniProt 条目名称
- `uniprotDescription`：蛋白质描述
- `taxId`：NCBI 分类标识符
- `organismScientificName`：物种学名
- `uniprotStart/uniprotEnd`：覆盖的残基范围
- `uniprotSequence`：完整蛋白质序列
- `modelCreatedDate`：初始预测日期
- `latestVersion`：当前模型版本号
- `allVersions`：可用版本列表
- `cifUrl/bcifUrl/pdbUrl`：结构文件下载 URL
- `paeImageUrl`：PAE 可视化图像 URL
- `paeDocUrl`：PAE 数据 JSON URL

### 2. 3D 信标集成

AlphaFold 集成到 3D-Beacons 网络中，用于联合结构访问。

**端点：** `https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/uniprot/summary/{uniprot_id}.json`

**示例：**
<<<代码块_3>>>

---

## 文件访问模式

### 直接文件下载

所有 AlphaFold 文件均可通过直接 URL 访问，无需身份验证。

**网址模式：**
<<<代码块_4>>>

**组件：**
- `{alphafold_id}`：条目标识符（例如“AF-P00520-F1”）
- `{file_type}`：文件类型（见下文）
- `{version}`：数据库版本（例如“v4”）
- `{extension}`：文件格式扩展名

### 可用的文件类型

#### 1.模型坐标

**mmCIF 格式（推荐）：**
<<<代码块_5>>>
- 标准晶体格式
- 包含完整的元数据
- 支持大型结构
- 文件大小：可变（典型值为 100KB - 10MB）

**二进制 CIF 格式：**
<<<代码块_6>>>
- mmCIF 的压缩二进制版本
- 文件大小更小（减少约 70%）
- 更快的解析
- 需要专门的解析器

**PDB 格式（旧版）：**
```
https://alphafold.ebi.ac.uk/files/AF-P00520-F1-model_v4.pdb
```
- 传统的PDB文本格式
- 仅限于 99,999 个原子
- 受到旧工具的广泛支持
- 文件大小：与 mmCIF 类似

#### 2. 置信度指标

**每个残基置信度 (JSON)：**
```
https://alphafold.ebi.ac.uk/files/AF-P00520-F1-confidence_v4.json
```

**结构：**
```json
{
  "confidenceScore": [87.5, 91.2, 93.8, ...],
  "confidenceCategory": ["high", "very_high", "very_high", ...]
}
```

**字段：**
- `confidenceScore`：每个残基的 pLDDT 值 (0-100) 数组
- `confidenceCategory`：分类（very_low、low、high、very_high）

#### 3. 预测对齐误差 (JSON)

```
https://alphafold.ebi.ac.uk/files/AF-P00520-F1-predicted_aligned_error_v4.json
```

**结构：**
```json
{
  "distance": [[0, 2.3, 4.5, ...], [2.3, 0, 3.1, ...], ...],
  "max_predicted_aligned_error": 31.75
}
```

**字段：**
- `distance`：PAE 值的 N×N 矩阵（以 Ångströms 为单位）
- `max_predicted_aligned_error`：矩阵中的最大 PAE 值

#### 4.PAE 可视化（PNG）

```
https://alphafold.ebi.ac.uk/files/AF-P00520-F1-predicted_aligned_error_v4.png
```
- 预渲染的 PAE 热图
- 可用于快速视觉评估
- 分辨率：根据蛋白质大小而变化

### 批量下载策略

为了有效地下载多个文件，请使用并发下载并进行适当的错误处理和速率限制，以尊重服务器资源。

---

## 数据模式

### 坐标文件 (mmCIF) 架构

AlphaFold mmCIF 文件包含：

**关键数据类别：**
- `_entry`：入门级元数据
- `_struct`：结构标题和描述
- `_entity`：分子实体信息
- `_atom_site`：原子坐标和属性
- `_pdbx_struct_assembly`：生物组装信息

**`_atom_site` 中的重要字段：**
- `group_PDB`：所有记录的“ATOM”
- `id`：原子序列号
- `label_atom_id`：原子名称（例如“CA”、“N”、“C”）
- `label_comp_id`：残基名称（例如“ALA”、“GLY”）
- `label_seq_id`：残基序列号
- `Cartn_x/y/z`：笛卡尔坐标 (Ångströms)
- `B_iso_or_equiv`：B 因子（包含 pLDDT 分数）
**B 因子列中的 pLDDT：**
AlphaFold 将每个残基置信度 (pLDDT) 存储在 B 因子字段中。这允许标准结构查看器自动按置信度着色。

### 置信度 JSON 架构

```json
{
  "confidenceScore": [
    87.5,   // Residue 1 pLDDT
    91.2,   // Residue 2 pLDDT
    93.8    // Residue 3 pLDDT
    // ... one value per residue
  ],
  "confidenceCategory": [
    "high",      // Residue 1 category
    "very_high", // Residue 2 category
    "very_high"  // Residue 3 category
    // ... one category per residue
  ]
}
```

**置信类别：**
- `very_high`：pLDDT > 90
- `high`：70 < pLDDT ≤ 90
- `low`：50 < pLDDT ≤ 70
- `very_low`：pLDDT ≤ 50

### PAE JSON 架构

```json
{
  "distance": [
    [0.0, 2.3, 4.5, ...],     // PAE from residue 1 to all residues
    [2.3, 0.0, 3.1, ...],     // PAE from residue 2 to all residues
    [4.5, 3.1, 0.0, ...]      // PAE from residue 3 to all residues
    // ... N×N matrix for N residues
  ],
  "max_predicted_aligned_error": 31.75
}
```

**释义：**
- `distance[i][j]`：如果预测结构和真实结构在残基 i 上对齐，则残基 j 的预期位置误差 (Ångströms)
- 值越低表示相对定位越有信心
- 对角线始终为 0（残基与其自身对齐）
- 矩阵不对称：距离[i][j] ≠ 距离[j][i]

---

## 谷歌云访问

AlphaFold DB 托管在 Google Cloud Platform 上，用于批量访问。

### 云存储桶

**存储桶：** `gs://public-datasets-deepmind-alphafold-v4`

**目录结构：**
```
gs://public-datasets-deepmind-alphafold-v4/
├── accession_ids.csv              # Index of all entries (13.5 GB)
├── sequences.fasta                # All protein sequences (16.5 GB)
└── proteomes/                     # Grouped by species (1M+ archives)
```

### 安装 gsutil

```bash
# Using pip
pip install gsutil

# Or install Google Cloud SDK
curl https://sdk.cloud.google.com | bash
```

### 下载蛋白质组

**按分类 ID：**

```bash
# Download all archives for a species
TAX_ID=9606  # Human
gsutil -m cp gs://public-datasets-deepmind-alphafold-v4/proteomes/proteome-tax_id-${TAX_ID}-*_v4.tar .
```

---

## BigQuery 架构

AlphaFold 元数据可在 BigQuery 中用于基于 SQL 的查询。

**数据集：** `bigquery-public-data.deepmind_alphafold`
**表：** `metadata`

### 关键字段

|领域 |类型 |描述 |
|--------|------|-------------|
| `entryId` |字符串 | AlphaFold 条目 ID |
| `uniprotAccession` |字符串 | UniProt 加入 |
| `gene` |字符串 |基因符号|
| `organismScientificName` |字符串 |物种学名|
| `taxId` |整数 | NCBI 分类 ID |
| `globalMetricValue` |浮动|整体质量指标|
| `fractionPlddtVeryHigh` |浮动| pLDDT ≥ 90 的分数 |
| `isReviewed` |布尔 | Swiss-Prot 审核状态 |
| `sequenceLength` |整数 |蛋白质序列长度|

### 查询示例

```sql
SELECT
  entryId,
  uniprotAccession,
  gene,
  fractionPlddtVeryHigh
FROM `bigquery-public-data.deepmind_alphafold.metadata`
WHERE
  taxId = 9606  -- Homo sapiens
  AND fractionPlddtVeryHigh > 0.8
  AND isReviewed = TRUE
ORDER BY fractionPlddtVeryHigh DESC
LIMIT 100;
```

---

## 最佳实践

### 1. 缓存策略

始终在本地缓存下载的文件，以避免重复下载。

### 2. 错误处理

通过针对暂时性故障的重试逻辑，对 API 请求实施强大的错误处理。

### 3.批量处理

为了处理许多蛋白质，请使用具有适当速率限制的并发下载。

### 4.版本管理

始终在代码中指定和跟踪数据库版本（当前：v4）。

---

## 错误处理

### 常见 HTTP 状态代码

|代码|意义|行动|
|------|---------|--------|
| 200 | 200成功|正常处理响应 |
| 404 | 404未找到 |没有针对此 UniProt ID 的 AlphaFold 预测 |
| 429 | 429太多请求 |实施速率限制并通过退避重试 |
| 500 | 500服务器错误 |使用指数退避重试 |
| 503 | 503服务不可用 |等待并稍后重试 |

---

## 速率限制

### 建议

- 最多限制为 **10 个并发请求**
- 在顺序请求之间添加 **100-200ms 延迟**
- 使用 Google Cloud 代替 REST API 进行批量下载
- 将所有下载的数据缓存在本地

---

## 其他资源

- **AlphaFold GitHub：** https://github.com/google-deepmind/alphafold
- **Google Cloud 文档：** https://cloud.google.com/datasets/alphafold
- **3D 信标文档：** https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/docs
- **Biopython 教程：** https://biopython.org/wiki/AlphaFold

## 版本历史

- **v1** (2021)：具有约 350K 结构的初始版本
- **v2** (2022)：扩展到 2 亿+ 结构
- **v3** (2023)：更新模型并扩大覆盖范围
- **v4** (2024)：具有改进的置信度指标的当前版本

## 引文

在出版物中使用 AlphaFold DB 时，请引用：

1.Jumper, J. 等人。使用 AlphaFold 进行高度准确的蛋白质结构预测。自然 596, 583–589 (2021)。
2.瓦拉迪，M.等人。 2024 年 AlphaFold 蛋白质结构数据库：提供超过 2.14 亿个蛋白质序列的结构覆盖。核酸研究。 52，D368–D375（2024）。