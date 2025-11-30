<!-- 此文件由机器翻译自 data_formats.md -->

# ClinVar 数据格式和 FTP 访问

## 概述

ClinVar 提供多种格式的批量数据下载，以支持不同的研究工作流程。数据通过 FTP 分发并定期更新。

## FTP 访问

### 基本网址
```
ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/
```

### 更新时间表

- **每月发布**：每个月的第一个星期四
  - 完整的数据集和全面的文档
  - 无限期存档以确保可重复性
  - 包括发行说明

- **每周更新**：每周一
  - 每月发布的增量更新
  - 保留到下个月发布
  - 允许与 ClinVar 网站同步

### 目录结构

<<<代码块_1>>>

## 数据格式

### 1. XML 格式（初次分发）

XML 提供了最全面的数据，包括完整的提交详细信息、证据和元数据。

#### VCV（变体）文件
- **目的**：以变体为中心的聚合
- **位置**：`xml/clinvar_variation/`
- **登录格式**：VCV000000001.1
- **最适合**：无论条件如何，查询都集中于特定变体
- **文件命名**：`ClinVarVariationRelease_YYYY-MM-DD.xml.gz`

**VCV 记录结构：**
<<<代码块_2>>>

#### RCV（记录）文件
- **目的**：变体条件对聚合
- **位置**：`xml/RCV/`
- **登录格式**：RCV000000001.1
- **最适合**：专注于变异疾病关系的查询
- **文件命名**：`ClinVarRCVRelease_YYYY-MM-DD.xml.gz`

**与 VCV 的主要区别：**
- 每个变体条件组合一个 RCV
- 单个变体可能有多个RCV记录（不同条件）
- 更注重每种疾病的临床解释

#### SCV（提交）记录
- **格式**：VCV/RCV 记录中的个人提交
- **登录格式**：SCV000000001.1
- **内容**：提交者特定的解释和证据

### 2.VCF格式

用于基因组分析流程的变体调用格式文件。

#### 地点
- **GRCh37/hg19**：`vcf_GRCh37/clinvar.vcf.gz`
- **GRCh38/hg38**：`vcf_GRCh38/clinvar.vcf.gz`

#### 内容限制
- **包含**：具有精确基因组坐标的简单等位基因
- **排除**：
  - 变体 >10 kb
  - 细胞遗传学变异
  - 复杂的结构变体
  - 没有精确断点的变体

#### VCF 信息字段

ClinVar VCF 中的关键信息字段：

|领域 |描述 |
|--------|-------------|
| **等位基因** | ClinVar 等位基因标识符 |
| **CLNSIG** |临床意义 |
| **CLNREVSTAT** |审核状态 |
| **CLNDN** |条件名称 |
| **CLNVC** |变体类型（SNV、删除等）|
| **CLNVCSO** |序列本体术语|
| **基因信息** |基因符号：基因ID |
| **主持人** |分子后果|
| **RS** | dbSNP rsID | dbSNP rsID |
| **AF_ESP** |等位基因频率 (ESP) |
| **AF_EXAC** |等位基因频率 (ExAC) |
| **AF_TGP** |等位基因频率（1000 个基因组）|

#### VCF 线示例
<<<代码块_3>>>

### 3. 制表符分隔格式

用于快速分析和数据库加载的摘要文件。

####variant_summary.txt
主要摘要文件，包含所有基因组映射变体的选定元数据。

**关键栏：**
- `VariationID` - ClinVar 变异标识符
- `Type` - 变体类型（SNV、indel、CNV 等）
- `Name` - 变体名称（通常为 HGVS）
- `GeneID` - NCBI 基因 ID
- `GeneSymbol` - 基因符号
- `ClinicalSignificance` - 分类
- `ReviewStatus` - 星级
- `LastEvaluated` - 上次审核日期
- `RS# (dbSNP)` - dbSNP rsID（如果可用）
- `Chromosome` - 染色体
- `PositionVCF` - 位置 (GRCh38)
- `ReferenceAlleleVCF` - 参考等位基因
- `AlternateAlleleVCF` - 替代等位基因
- `Assembly` - 参考程序集 (GRCh37/GRCh38)
- `PhenotypeIDS` - MedGen/OMIM/Orphanet ID
- `Origin` - 种系、体细胞、从头等。
- `SubmitterCategories` - 提交者类型（临床、研究等）

**用法示例：**
<<<代码块_4>>>

#### var_引用.txt
对 PubMed 文章、dbSNP 和 dbVar 的交叉引用。

**列：**
- `AlleleID` - ClinVar 等位基因 ID
- `VariationID` - ClinVar 变异 ID
- `rs` - dbSNP rsID
- `nsv/esv` - dbVar ID
- `PubMedID` - PubMed 引文

#### cross_references.txt
带有修改日期的数据库交叉引用。

**列：**
- `VariationID`
- `Database`（OMIM、UniProtKB、GTR 等）
- `Identifier`
- `DateLastModified`

## 选择正确的格式

### 在以下情况下使用 XML：
- 需要完整的提交详细信息
- 想要跟踪证据和标准
- 建立全面的变异数据库
- 需要完整的元数据和关系

### 在以下情况下使用 VCF：
- 与基因组分析流程集成
- 注释测序中的变体调用
- 需要基因组坐标进行重叠分析
- 使用标准生物信息学工具

### 在以下情况下使用制表符分隔：
- 快速数据库查询和过滤
- 加载到电子表格或数据库中
- 简单的数据提取和统计
- 不需要完整的证据细节

## 种质类型和标识符

### VCV（变体存档）
- **格式**：VCV000012345.6（ID.版本）
- **范围**：聚合单个变体的所有数据
- **版本控制**：变体数据更改时增加

### RCV（记录）
- **格式**：RCV000056789.4
- **范围**：一种变体条件解释
- **版本控制**：解释更改时增加

### SCV（提交）
- **格式**：SCV000098765.2
- **范围**：个别提交者的解释
- **版本控制**：提交更新时增加

### 其他标识符
- **VariationID**：变体的稳定数字标识符
- **AlleleID**：等位基因的稳定数字标识符
- **dbSNP rsID**：对 dbSNP 的交叉引用（如果可用）

## 文件处理技巧

### XML 处理

**Python 与 xml.etree:**
<<<代码块_5>>>

**带有 xmllint 的命令行：**
<<<代码块_6>>>

### VCF 处理

**使用 bcftools：**
```bash
# Filter by clinical significance
bcftools view -i 'INFO/CLNSIG~"Pathogenic"' clinvar.vcf.gz

# Extract specific genes
bcftools view -i 'INFO/GENEINFO~"BRCA"' clinvar.vcf.gz

# Annotate your VCF
bcftools annotate -a clinvar.vcf.gz -c INFO your_variants.vcf
```

**使用 PyVCF：**
```python
import vcf

vcf_reader = vcf.Reader(filename='clinvar.vcf.gz')
for record in vcf_reader:
    clnsig = record.INFO.get('CLNSIG', [])
    if 'Pathogenic' in clnsig:
        print(f"{record.CHROM}:{record.POS} - {clnsig}")
```

### 制表符分隔处理

**使用熊猫：**
```python
import pandas as pd

# Read variant summary
df = pd.read_csv('variant_summary.txt.gz', sep='\t', compression='gzip')

# Filter pathogenic variants
pathogenic = df[df['ClinicalSignificance'].str.contains('Pathogenic', na=False)]

# Group by gene
gene_counts = pathogenic.groupby('GeneSymbol').size().sort_values(ascending=False)
```

## 数据质量注意事项

### 已知限制

1. **VCF 文件排除大变体** - 不包括 >10 kb 的变体
2. **历史数据可能不太准确** - 较早提交的标准化要求较少
3. **存在冲突的解释** - 多个提交者可能不同意
4. **并非所有变体都有基因组坐标** - 某些 HGVS 表达无法映射

### 验证建议

- 尽可能交叉引用多种数据格式
- 检查评论状态（首选 ★★★ 或 ★★★★ 评级）
- 根据当前基因组构建验证基因组坐标
- 考虑人口频率数据（gnomAD）作为背景
- 审查提交日期 - 更新的数据可能更准确

## 批量下载脚本

### 下载最新的每月版本

```bash
#!/bin/bash
# Download latest ClinVar monthly XML release

BASE_URL="ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/clinvar_variation"

# Get latest file
LATEST=$(curl -s ${BASE_URL}/ | \
         grep -oP 'ClinVarVariationRelease_\d{4}-\d{2}\.xml\.gz' | \
         tail -1)

# Download
wget ${BASE_URL}/${LATEST}
```

### 下载所有格式

```bash
#!/bin/bash
# Download ClinVar in all formats

FTP_BASE="ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar"

# XML
wget ${FTP_BASE}/xml/clinvar_variation/ClinVarVariationRelease_00-latest.xml.gz

# VCF (both assemblies)
wget ${FTP_BASE}/vcf_GRCh37/clinvar.vcf.gz
wget ${FTP_BASE}/vcf_GRCh38/clinvar.vcf.gz

# Tab-delimited
wget ${FTP_BASE}/tab_delimited/variant_summary.txt.gz
wget ${FTP_BASE}/tab_delimited/var_citations.txt.gz
```

## 其他资源

- ClinVar FTP 入门：https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/
- XML 架构文档：https://www.ncbi.nlm.nih.gov/clinvar/docs/xml_schemas/
- VCF 规范：https://samtools.github.io/hts-specs/VCFv4.3.pdf
- 发行说明：https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/README.txt