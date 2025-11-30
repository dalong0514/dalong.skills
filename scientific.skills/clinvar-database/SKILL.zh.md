<!-- 此文件由机器翻译自 SKILL.md -->

---
名称： clinvar-数据库
描述：“查询 NCBI ClinVar 以了解不同的临床意义。按基因/位置搜索，解释致病性分类，通过电子实用程序 API 或 FTP 访问，注释 VCF，用于基因组医学。”
---

# ClinVar 数据库

## 概述

ClinVar 是 NCBI 可免费访问的关于人类遗传变异和表型之间关系的报告档案，并附有支持证据。该数据库汇总了有关基因组变异及其与人类健康关系的信息，提供临床遗传学和研究中使用的标准化变异分类。

## 何时使用此技能

该技能应该在以下情况下使用：

- 按基因、状况或临床意义搜索变异
- 解释临床意义分类（致病性、良性、VUS）
- 通过电子实用程序 API 以编程方式访问 ClinVar 数据
- 从FTP下载和处理批量数据
- 了解评论状态和星级评定
- 解决相互冲突的变体解释
- 注释具有临床意义的变异调用集

## 核心能力

### 1.搜索查询ClinVar

#### Web 界面查询

使用网络界面搜索 ClinVar，网址为 https://www.ncbi.nlm.nih.gov/clinvar/

**常见搜索模式：**
- 按基因：`BRCA1[gene]`
- 按临床意义：`pathogenic[CLNSIG]`
- 按条件：`breast cancer[disorder]`
- 按变体：`NM_000059.3:c.1310_1313del[variant name]`
- 按染色体：`13[chr]`
- 组合：`BRCA1[gene] AND pathogenic[CLNSIG]`

#### 通过电子实用程序进行编程访问

使用 NCBI 的电子实用程序 API 以编程方式访问 ClinVar。请参阅 `references/api_reference.md` 以获取全面的 API 文档，包括：
- **研究** - 搜索匹配条件的变体
- **摘要** - 检索变体摘要
- **efetch** - 下载完整的 XML 记录
- **elink** - 在其他 NCBI 数据库中查找相关记录

**使用curl的快速示例：**
```bash
# Search for pathogenic BRCA1 variants
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=BRCA1[gene]+AND+pathogenic[CLNSIG]&retmode=json"
```

**最佳实践：**
- 在自动化之前在 Web 界面上测试查询
- 使用 API 密钥将速率限制从 3 个请求/秒提高到 10 个请求/秒
- 对速率限制错误实施指数退避
- 使用 Biopython 时设置 `Entrez.email`

### 2. 解释临床意义

#### 了解分类

ClinVar 使用标准化术语进行变异分类。有关详细解释指南，请参阅`references/clinical_significance.md`。

**关键种系分类术语 (ACMG/AMP)：**
- **致病性 (P)** - 变异导致疾病（~99% 的概率）
- **可能致病 (LP)** - 变异可能导致疾病（~90% 概率）
- **意义不确定 (VUS)** - 证据不足，无法分类
- **可能良性 (LB)** - 变异可能不会引起疾病
- **良性 (B)** - 变异不会引起疾病

**评论状态（星级）：**
- ★★★★ 练习指南 - 最高置信度
- ★★★ 专家小组评审（例如 ClinGen） - 高可信度
- ★★ 多个提交者，无冲突 - 中等信心
- ★ 有标准的单一提交者 - 标准重量
- ☆ 无断言标准 - 低置信度

**关键考虑因素：**
- 经常检查评论状态 - 更喜欢 ★★★ 或 ★★★★ 评级
- 相互矛盾的解释需要手动评估
- 随着新证据的出现，分类可能会发生变化
- VUS（意义不确定）变异缺乏足够的临床使用证据

### 3.从FTP下载批量数据

#### 访问 ClinVar FTP 站点

从`ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/`下载完整数据集

有关文件格式和处理的综合文档，请参阅`references/data_formats.md`。

**更新时间表：**
- 每月发布：每月第一个星期四（完整数据集，已存档）
- 每周更新：每周一（增量更新）

#### 可用格式

**XML 文件**（最全面）：
- VCV（变体）文件：`xml/clinvar_variation/` - 以变体为中心的聚合
- RCV（记录）文件：`xml/RCV/` - 变体条件对
- 包括完整的提交详细信息、证据和元数据

**VCF 文件**（用于基因组管道）：
- GRCh37：`vcf_GRCh37/clinvar.vcf.gz`
- GRCh38：`vcf_GRCh38/clinvar.vcf.gz`
- 限制：不包括 >10kb 的变体和复杂结构的变体

**制表符分隔文件**（用于快速分析）：
- `tab_delimited/variant_summary.txt.gz` - 所有变体的摘要
- `tab_delimited/var_citations.txt.gz` - PubMed 引文
- `tab_delimited/cross_references.txt.gz` - 数据库交叉引用

**示例下载：**
<<<代码块_1>>>

### 4. 处理和分析 ClinVar 数据

#### 使用 XML 文件

处理 XML 文件以提取变体详细信息、分类和证据。
**带有 xml.etree 的 Python 示例：**
<<<代码块_2>>>

#### 使用 VCF 文件

使用 bcftools 或 Python 对变异调用进行注释或按临床意义进行筛选。

**使用 bcftools：**
<<<代码块_3>>>

**在Python中使用PyVCF：**
<<<代码块_4>>>

#### 使用制表符分隔的文件

使用 pandas 或命令行工具进行快速过滤和分析。

**使用熊猫：**
<<<代码块_5>>>

**使用命令行工具：**
<<<代码块_6>>>

### 5. 处理相互冲突的解释

当多个提交者为同一变体提供不同的分类时，ClinVar 会报告“致病性的解释相互矛盾”。

**解决策略：**
1. 检查评论状态（星级）——评级越高，权重越大
2. 检查每个提交者的证据和断言标准
3. 考虑提交日期——较新的提交可能反映更新的证据
4. 查看人口频率数据（例如，gnomAD）以了解背景
5. 如果有的话，请咨询专家小组分类（★★★）
6. 对于临床使用，请务必遵循遗传学专业人士的意见

**搜索查询以排除冲突：**
```
TP53[gene] AND pathogenic[CLNSIG] NOT conflicting[RVSTAT]
```

### 6. 赛道分类更新

随着新证据的出现，变异分类可能会随着时间而改变。

**分类为何变化：**
- 新的功能研究或临床数据
- 更新了人口频率信息
- 修订了 ACMG/AMP 指南
- 来自其他家庭的隔离数据

**最佳实践：**
- 记录 ClinVar 版本和访问日期以确保可重复性
- 定期重新检查分类中的关键变体
- 订阅 ClinVar 邮件列表以获取主要更新
- 使用每月存档的版本来获得稳定的数据集

### 7. 向 ClinVar 提交数据

组织可以向 ClinVar 提交变体解释。

**提交方式：**
- 网络提交门户：https://submit.ncbi.nlm.nih.gov/subs/clinvar/
- API提交（需要服务帐户）：参见`references/api_reference.md`
- 通过Excel模板批量提交

**要求：**
- NCBI 的组织帐户
- 断言标准（最好是 ACMG/AMP 指南）
- 分类的支持证据

请联系：clinvar@ncbi.nlm.nih.gov 以设置提交帐户。

## 工作流程示例

### 示例 1：识别基因中高可信度的致病变异

**目标：** 通过专家小组审查找到 CFTR 基因的致病性变异。

**步骤：**
1. 使用网络界面或电子实用程序搜索：
   ```
   CFTR[gene] AND pathogenic[CLNSIG] AND (reviewed by expert panel[RVSTAT] OR practice guideline[RVSTAT])
   ```
2. 审核结果，注明审核状态（应为 ★★★ 或 ★★★★）
3.导出变体列表或通过efetch检索完整记录
4. 与临床表现交叉参考（如果适用）

### 示例 2：使用 ClinVar 分类注释 VCF

**目标：** 为变异调用添加临床意义注释。

**步骤：**
1. 下载适当的 ClinVar VCF（匹配基因组构建：GRCh37 或 GRCh38）：
   ```bash
   wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
   wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
   ```
2.使用bcftools进行注释：
   ```bash
   bcftools annotate -a clinvar.vcf.gz \
     -c INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT \
     -o annotated_variants.vcf \
     your_variants.vcf
   ```
3. 过滤带注释的 VCF 是否存在致病性变异：
   ```bash
   bcftools view -i 'INFO/CLNSIG~"Pathogenic"' annotated_variants.vcf
   ```

### 示例 3：分析特定疾病的变异

**目标：** 研究与遗传性乳腺癌相关的所有变异。

**步骤：**
1、按条件搜索：
   ```
   hereditary breast cancer[disorder] OR "Breast-ovarian cancer, familial"[disorder]
   ```
2. 下载 CSV 格式的结果或通过电子实用程序检索
3. 按审核状态过滤以优先考虑高置信度变体
4. 分析基因分布（BRCA1、BRCA2、PALB2 等）
5. 分别检查具有冲突解释的变体

### 示例4：批量下载和数据库构建

**目标：** 构建本地 ClinVar 数据库用于分析流程。

**步骤：**
1. 下载每月版本以实现可重复性：
   ```bash
   wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/clinvar_variation/ClinVarVariationRelease_YYYY-MM.xml.gz
   ```
2. 解析XML并加载到数据库（PostgreSQL、MySQL、MongoDB）
3. 按基因、位置、临床意义、审查状态分类的索引
4. 实现更新的版本跟踪
5. 安排每月从 FTP 站点更新

## 重要限制和注意事项

### 数据质量
- **并非所有提交的内容都具有相同的权重** - 检查审核状态（星级）
- **存在冲突的解释** - 需要手动评估
- **历史提交内容可能已过时** - 较新的数据可能更准确
- **VUS 分类不是临床诊断** - 意味着证据不足

### 范围限制
- **不适用于直接临床诊断** - 始终让遗传学专业人士参与
- **特定人群** - 变异频率因血统而异
- **覆盖不完整** - 并非所有基因或变异都经过充分研究
- **版本依赖性** - 跨分析协调基因组构建 (GRCh37/GRCh38)

### 技术限制
- **VCF 文件排除大变体** - >10kb 的变体不采用 VCF 格式
- **API 的速率限制** - 不带密钥时为 3 请求/秒，带 API 密钥时为 10 请求/秒
- **文件大小** - 完整的 XML 版本是多 GB 压缩文件
- **无实时更新** - 网站每周更新，FTP 每月/每周更新

## 资源

### 参考文档

该技能包括全面的参考文档：

- **`references/api_reference.md`** - 完整的电子实用程序 API 文档，包含 esearch、esummary、efetch 和 elink 示例；包括速率限制、身份验证和 Python/Biopython 代码示例

- **`references/clinical_significance.md`** - 解释临床意义分类、审查状态星级、冲突解决和变异解释最佳实践的详细指南

- **`references/data_formats.md`** - XML、VCF 和制表符分隔文件格式的文档； FTP目录结构、处理示例和格式选择指导

### 外部资源

- ClinVar 主页：https://www.ncbi.nlm.nih.gov/clinvar/
- ClinVar 文档：https://www.ncbi.nlm.nih.gov/clinvar/docs/
- 电子实用程序文档：https://www.ncbi.nlm.nih.gov/books/NBK25501/
- ACMG 变异解释指南：Richards 等人，2015（PMID：25741868）
- ClinGen 专家小组：https://clinicalgenome.org/

### 联系方式

有关 ClinVar 或数据提交的问题：clinvar@ncbi.nlm.nih.gov