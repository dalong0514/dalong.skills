<!-- 此文件由机器翻译自 cosmic_data_reference.md -->

# COSMIC 数据库参考

## 概述

COSMIC（癌症体细胞突变目录）是世界上最大、最全面的资源，用于探索体细胞突变对人类癌症的影响。它由 Wellcome Sanger Institute 维护，对数千种癌症类型的数百万个突变进行了分类。

**网站**：https://cancer.sanger.ac.uk/cosmic
**发布时间表**：季度更新
**当前版本**：v102（2025 年 5 月），在 API 调用中使用“最新”表示最新版本

## 数据访问

### 身份验证
- **学术用户**：免费访问（需要注册）
- **商业用户**：需要许可证（联系 QIAGEN）
- **注册**：https://cancer.sanger.ac.uk/cosmic/register

### 下载方法
1. **网络浏览器**：在 https://cancer.sanger.ac.uk/cosmic 进行交互式搜索
2. **文件下载**：通过下载 API 进行编程访问
3. **数据文件**：TSV、CSV 和 VCF 格式

## 可用数据类型

### 1. 核心突变数据
**主要文件**：
- `CosmicMutantExport.tsv.gz` - 完整的编码突变
- `CosmicCodingMuts.vcf.gz` - VCF 格式的突变
- `CosmicNonCodingVariants.vcf.gz` - 非编码变体
- `CosmicMutantExportCensus.tsv.gz` - 仅癌症基因普查基因突变

**内容**：
- 点突变（SNV）
- 小插入和缺失（indels）
- 基因组坐标
- 变体注释
- 样本信息
- 肿瘤类型关联

### 2.癌症基因普查
**文件**：`cancer_gene_census.csv`

**内容**：
- 专家整理的癌症基因列表
- 约 700 多个基因有大量证据表明与癌症有关
- 基因作用（癌基因、抑癌基因、融合基因）
- 突变类型
- 组织协会
- 分子遗传学信息

### 3. 突变特征
**文件**：位于 `signatures/` 目录中
- `signatures.tsv` - 签名定义
- 单碱基取代 (SBS) 签名
- 双碱基取代 (DBS) 签名
- 插入/删除（ID）签名

**当前版本**：v3.4（在 COSMIC v98 中发布）

**内容**：
- 签名配置文件（96通道、78通道、83通道）
- 病因学注释
- 签名分析参考签名

### 4. 结构变体
**文件**：`CosmicStructExport.tsv.gz`

**内容**：
- 基因融合
- 结构断点
- 易位事件
- 大量删除/插入
- 复杂的重新排列

### 5. 拷贝数变化
**文件**：`CosmicCompleteCNA.tsv.gz`

**内容**：
- 复制数量得失
- 扩充和删除
- 段级数据
- 基因级注释

### 6. 基因表达
**文件**：`CosmicCompleteGeneExpression.tsv.gz`

**内容**：
- 表达过度/表达不足的数据
- 基因表达 Z 分数
- 组织特异性表达模式

### 7. 抗性突变
**文件**：`CosmicResistanceMutations.tsv.gz`

**内容**：
- 耐药突变
- 治疗协会
- 临床相关性

### 8. 细胞系项目
**文件**：各种细胞系特定文件

**内容**：
- 癌细胞系突变
- 复制细胞系的数字数据
- 细胞系中的融合基因
- 微卫星不稳定状态

### 9. 示例信息
**文件**：`CosmicSample.tsv.gz`

**内容**：
- 元数据示例
- 肿瘤部位/组织学
- 样本来源
- 研究参考资料

## 基因组组装

所有基因组数据可用于两个参考基因组：
- **GRCh37** (hg19) - 旧版组件
- **GRCh38** (hg38) - 当前组件（推荐）

文件路径使用以下模式：`{assembly}/cosmic/{version}/{filename}`

## 文件格式

### TSV/CSV 格式
- 制表符或逗号分隔值
- 包括列标题
- Gzip 压缩 (.gz)
- 可以使用 pandas、awk 或标准工具读取

### VCF 格式
- 标准变体调用格式
- 4.x 版规范
- 包括带有 COSMIC 注释的 INFO 字段
- Gzip 压缩和索引（.vcf.gz、.vcf.gz.tbi）

## 常用文件路径

使用 `latest` 作为最新版本：

```
# Coding mutations (TSV)
GRCh38/cosmic/latest/CosmicMutantExport.tsv.gz

# Coding mutations (VCF)
GRCh38/cosmic/latest/VCF/CosmicCodingMuts.vcf.gz

# Cancer Gene Census
GRCh38/cosmic/latest/cancer_gene_census.csv

# Structural variants
GRCh38/cosmic/latest/CosmicStructExport.tsv.gz

# Copy number alterations
GRCh38/cosmic/latest/CosmicCompleteCNA.tsv.gz

# Gene fusions
GRCh38/cosmic/latest/CosmicFusionExport.tsv.gz

# Gene expression
GRCh38/cosmic/latest/CosmicCompleteGeneExpression.tsv.gz

# Resistance mutations
GRCh38/cosmic/latest/CosmicResistanceMutations.tsv.gz

# Mutational signatures
signatures/signatures.tsv

# Sample information
GRCh38/cosmic/latest/CosmicSample.tsv.gz
```

## 关键数据字段

### 突变数据字段
- **基因名称** - HGNC 基因符号
- **登录号** - 转录本标识符
- **COSMIC ID** - 唯一突变标识符
- **CDS 突变** - 编码序列改变
- **AA突变** - 氨基酸变化
- **原发部位** - 肿瘤解剖位置
- **主要组织学** - 肿瘤类型分类
- **基因组坐标** - 染色体、位置、链
- **突变类型** - 替换、插入、删除等。
- **接合性** - 杂合/纯合状态
- **Pubmed ID** - 文献参考

### 癌症基因普查领域
- **基因符号** - 官方基因名称
- **Entrez GeneId** - NCBI 基因标识符
- **在癌症中的作用** - 癌基因、TSG、融合
- **突变类型** - 观察到的改变类型
- **易位伙伴** - 用于融合基因
- **层** - 证据分类（1 或 2）
- **标志** - 癌症标志协会
- **体细胞** - 是否记录体细胞突变
- **种系** - 是否记录了种系突变

## 数据更新

COSMIC 每季度更新一次新版本。每个版本包括：
- 来自文献和数据库的新突变数据
- 更新了癌症基因普查注释
- 修改突变特征（如果适用）
- 增强的示例注释

## 引文

使用 COSMIC 数据时，引用：
Tate JG、Bamford S、Jubb HC 等。 COSMIC：癌症体细胞突变目录。核酸研究。 2019；47（D1）：D941-D947。

## 其他资源

- **文档**：https://cancer.sanger.ac.uk/cosmic/help
- **发行说明**：https://cancer.sanger.ac.uk/cosmic/release_notes
- **联系方式**：cosmic@sanger.ac.uk
- **许可**：cosmic-translation@sanger.ac.uk