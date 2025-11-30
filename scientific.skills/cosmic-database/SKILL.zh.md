<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：宇宙数据库
描述：“访问 COSMIC 癌症突变数据库。查询体细胞突变、癌症基因普查、突变特征、基因融合，用于癌症研究和精准肿瘤学。需要身份验证。”
---

# 宇宙数据库

## 概述

COSMIC（癌症体细胞突变目录）是世界上最大、最全面的探索人类癌症体细胞突变的数据库。以编程方式访问 COSMIC 广泛收集的癌症基因组数据，包括数千种癌症类型的数百万个突变、精心策划的基因列表、突变特征和临床注释。

## 何时使用此技能

该技能应该在以下情况下使用：
- 从 COSMIC 下载癌症突变数据
- 访问癌症基因普查以获取精选的癌症基因列表
- 检索突变特征谱
- 查询结构变异、拷贝数改变或基因融合
- 分析耐药突变
- 使用癌细胞系基因组数据
- 将癌症突变数据整合到生物信息学流程中
- 研究癌症背景下的特定基因或突变

## 先决条件

### 账户注册
COSMIC 需要身份验证才能下载数据：
- **学术用户**：通过 https://cancer.sanger.ac.uk/cosmic/register 注册即可免费访问
- **商业用户**：需要许可证（联系 QIAGEN）

### Python 要求
```bash
uv pip install requests pandas
```

## 快速入门

### 1. 基本文件下载

使用 `scripts/download_cosmic.py` 脚本下载 COSMIC 数据文件：

<<<代码块_1>>>

### 2. 命令行使用

<<<代码块_2>>>

### 3. 使用下载的数据

<<<代码块_3>>>

## 可用数据类型

### 核心突变
下载全面的突变数据，包括点突变、插入缺失和基因组注释。

**常见数据类型**：
- `mutations` - 完整编码突变（TSV 格式）
- `mutations_vcf` - 以 VCF 格式编码突变
- `sample_info` - 样本元数据和肿瘤信息

<<<代码块_4>>>

### 癌症基因普查
访问专家整理的约 700 多个癌症基因列表，其中包含癌症参与的大量证据。

<<<代码块_5>>>

**用例**：
- 识别已知的癌症基因
- 按癌症相关性过滤变体
- 了解基因作用（癌基因与肿瘤抑制基因）
- 研究目的基因选择

### 突变特征
下载签名配置文件以进行突变签名分析。

<<<代码块_6>>>

**签名类型**：
- 单碱基取代 (SBS) 签名
- 双碱基取代 (DBS) 签名
- 插入/删除（ID）签名

### 结构变体和融合
访问基因融合数据和结构重排。

**可用数据类型**：
- `structural_variants` - 结构断点
- `fusion_genes` - 基因融合事件

```python
# Download gene fusions
download_cosmic_file(
    email="user@email.com",
    password="password",
    filepath="GRCh38/cosmic/latest/CosmicFusionExport.tsv.gz"
)
```

### 拷贝数和表达式
检索拷贝数改变和基因表达数据。

**可用数据类型**：
- `copy_number` - 复制数字增益/损失
- `gene_expression` - 表达过度/不足的数据

```python
# Download copy number data
download_cosmic_file(
    email="user@email.com",
    password="password",
    filepath="GRCh38/cosmic/latest/CosmicCompleteCNA.tsv.gz"
)
```

### 抗性突变
通过临床注释获取耐药突变数据。

```python
# Download resistance mutations
download_cosmic_file(
    email="user@email.com",
    password="password",
    filepath="GRCh38/cosmic/latest/CosmicResistanceMutations.tsv.gz"
)
```

## 使用 COSMIC 数据

### 基因组组装
COSMIC 提供两个参考基因组的数据：
- **GRCh38**（推荐，当前标准）
- **GRCh37**（旧版，适用于较旧的管道）

在文件路径中指定程序集：
```python
# GRCh38 (recommended)
filepath="GRCh38/cosmic/latest/CosmicMutantExport.tsv.gz"

# GRCh37 (legacy)
filepath="GRCh37/cosmic/latest/CosmicMutantExport.tsv.gz"
```

### 版本控制
- 在文件路径中使用 `latest` 始终获取最新版本
- COSMIC 每季度更新一次（当前版本：v102，2025 年 5 月）
- 可以使用特定版本来实现再现性：`v102`、`v101` 等。

### 文件格式
- **TSV/CSV**：制表符/逗号分隔，gzip 压缩，用 pandas 读取
- **VCF**：标准变体格式，与 pysam、bcftools 或 GATK 一起使用
- 所有文件都包含描述列内容的标题

### 常见分析模式

**按基因过滤突变**：
```python
import pandas as pd

mutations = pd.read_csv('cosmic_mutations.tsv.gz', sep='\t', compression='gzip')
tp53_mutations = mutations[mutations['Gene name'] == 'TP53']
```

**按角色识别癌症基因**：
```python
gene_census = pd.read_csv('cancer_gene_census.csv')
oncogenes = gene_census[gene_census['Role in Cancer'].str.contains('oncogene', na=False)]
tumor_suppressors = gene_census[gene_census['Role in Cancer'].str.contains('TSG', na=False)]
```

**按癌症类型提取突变**：
```python
mutations = pd.read_csv('cosmic_mutations.tsv.gz', sep='\t', compression='gzip')
lung_mutations = mutations[mutations['Primary site'] == 'lung']
```

**使用 VCF 文件**：
```python
import pysam

vcf = pysam.VariantFile('CosmicCodingMuts.vcf.gz')
for record in vcf.fetch('17', 7577000, 7579000):  # TP53 region
    print(record.id, record.ref, record.alts, record.info)
```

## 数据参考

有关 COSMIC 数据结构、可用文件和字段描述的全面信息，请参阅`references/cosmic_data_reference.md`。该参考资料包括：

- 可用数据类型和文件的完整列表
- 每种文件类型的详细字段描述
- 文件格式规范
- 通用文件路径和命名约定
- 数据更新时间表和版本控制
- 引文信息

在以下情况下使用此参考：
- 探索 COSMIC 中可用的数据
- 理解特定字段的含义
- 确定数据类型的正确文件路径
- 使用 COSMIC 数据规划分析工作流程

## 辅助函数

下载脚本包含常见操作的辅助函数：

### 获取常用文件路径
```python
from scripts.download_cosmic import get_common_file_path

# Get path for mutations file
path = get_common_file_path('mutations', genome_assembly='GRCh38')
# Returns: 'GRCh38/cosmic/latest/CosmicMutantExport.tsv.gz'

# Get path for gene census
path = get_common_file_path('gene_census')
# Returns: 'GRCh38/cosmic/latest/cancer_gene_census.csv'
```

**可用的快捷键**：
- `mutations` - 核心编码突变
- `mutations_vcf` - VCF 格式突变
- `gene_census` - 癌症基因普查
- `resistance_mutations` - 耐药性数据
- `structural_variants` - 结构变体
- `gene_expression` - 表达式数据
- `copy_number` - 副本编号更改
- `fusion_genes` - 基因融合
- `signatures` - 突变签名
- `sample_info` - 示例元数据

## 故障排除

### 身份验证错误
- 验证电子邮件和密码是否正确
- 确保帐户已在 cancer.sanger.ac.uk/cosmic 注册
- 检查您的用例是否需要商业许可证

### 未找到文件
- 验证文件路径是否正确
- 检查请求的版本是否存在
- 使用`latest`作为最新版本
- 确认基因组组装（GRCh37 与 GRCh38）正确

### 大文件下载
- COSMIC 文件大小可达数 GB
- 确保有足够的磁盘空间
- 下载可能需要几分钟，具体取决于连接情况
- 脚本显示大文件的下载进度

### 商业用途
- 商业用户必须通过 QIAGEN 获得 COSMIC 许可
- 联系方式：cosmic-translation@sanger.ac.uk
- 学术访问是免费的，但需要注册

## 与其他工具集成

COSMIC 数据与以下各项完美集成：
- **变异注释**：VEP、ANNOVAR、SnpEff
- **签名分析**：SigProfiler、deconstructSigs、MuSiCa
- **癌症基因组学**：cBioPortal、OncoKB、CIViC
- **生物信息学**：Bioconductor、TCGA 分析工具
- **数据科学**：pandas、scikit-learn、PyTorch

## 其他资源

- **COSMIC 网站**：https://cancer.sanger.ac.uk/cosmic
- **文档**：https://cancer.sanger.ac.uk/cosmic/help
- **发行说明**：https://cancer.sanger.ac.uk/cosmic/release_notes
- **联系方式**：cosmic@sanger.ac.uk

## 引文

使用 COSMIC 数据时，引用：
Tate JG、Bamford S、Jubb HC 等。 COSMIC：癌症体细胞突变目录。核酸研究。 2019；47（D1）：D941-D947。