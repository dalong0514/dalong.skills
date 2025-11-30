<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：gwas数据库
描述：“查询 NHGRI-EBI GWAS 目录中的 SNP 性状关联。通过 rs ID、疾病/性状、基因搜索变体，检索 p 值和汇总统计数据，以了解遗传流行病学和多基因风险评分。”
---

# GWAS 目录数据库

## 概述

GWAS 目录是由国家人类基因组研究所 (NHGRI) 和欧洲生物信息学研究所 (EBI) 维护的已发表的全基因组关联研究的综合存储库。该目录包含来自数千份 GWAS 出版物的精选 SNP 性状关联，包括遗传变异、相关性状和疾病、p 值、效应大小以及许多研究的完整摘要统计数据。

## 何时使用此技能

当查询涉及以下内容时应使用此技能：

- **遗传变异关联**：寻找与疾病或性状相关的 SNP
- **SNP 查找**：检索有关特定遗传变异 (rs ID) 的信息
- **性状/疾病搜索**：发现表型的遗传关联
- **基因关联**：寻找特定基因中或附近的变异
- **GWAS 摘要统计**：访问完整的全基因组关联数据
- **研究元数据**：检索出版物和队列信息
- **群体遗传学**：探索特定血统的关联
- **多基因风险评分**：识别风险预测模型的变体
- **功能基因组学**：了解变异效应和基因组背景
- **系统评论**：遗传关联的综合文献综合

## 核心能力

### 1.了解 GWAS 目录数据结构

GWAS 目录围绕四个核心实体组织：

- **研究**：带有元数据的 GWAS 出版物（PMID、作者、队列详细信息）
- **关联**：SNP 性状与统计证据的关联 (p ≤ 5×10⁻⁸)
- **变体**：具有基因组坐标和等位基因的遗传标记 (SNP)
- **性状**：表型和疾病（映射到 EFO 本体术语）

**关键标识符：**
- 研究种质：`GCST` ID（例如 GCST001234）
- 变体 ID：`rs` 数字（例如 rs7903146）或 `variant_id` 格式
- 特征 ID：EFO 术语（例如，EFO_0001360 表示 2 型糖尿病）
- 基因符号：HGNC 批准的名称（例如 TCF7L2）

### 2. Web 界面搜索

https://www.ebi.ac.uk/gwas/ 的网络界面支持多种搜索模式：

**按变体（rs ID）：**
```
rs7903146
```
返回此 SNP 的所有性状关联。

**按疾病/性状：**
<<<代码块_1>>>
返回所有相关的遗传变异。

**按基因：**
<<<代码块_2>>>
返回基因区域内或附近的变体。

**按染色体区域：**
<<<代码块_3>>>
返回指定基因组区间内的变异。

**按出版物：**
<<<代码块_4>>>
返回研究详细信息和所有报告的关联。

### 3.REST API 访问

GWAS 目录提供了两个用于编程访问的 REST API：

**基本网址：**
- GWAS 目录 API：`https://www.ebi.ac.uk/gwas/rest/api`
- 摘要统计 API：`https://www.ebi.ac.uk/gwas/summary-statistics/api`

**API文档：**
- 主要 API 文档：https://www.ebi.ac.uk/gwas/rest/docs/api
- 摘要统计文档：https://www.ebi.ac.uk/gwas/summary-statistics/docs/

**核心端点：**

1. **研究端点** - `/studies/{accessionID}`
   <<<代码块_5>>>

2. **关联端点** - `/associations`
   <<<代码块_6>>>

3. **变体端点** - `/singleNucleotidePolymorphisms/{rsID}`
   ```python
   # Get variant details
   url = "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/rs7903146"
   response = requests.get(url, headers={"Content-Type": "application/json"})
   variant_info = response.json()
   ```

4. **特征端点** - `/efoTraits/{efoID}`
   ```python
   # Get trait information
   url = "https://www.ebi.ac.uk/gwas/rest/api/efoTraits/EFO_0001360"
   response = requests.get(url, headers={"Content-Type": "application/json"})
   trait_info = response.json()
   ```

### 4. 查询示例和模式

**示例 1：查找与某种疾病的所有关联**
```python
import requests

trait = "EFO_0001360"  # Type 2 diabetes
base_url = "https://www.ebi.ac.uk/gwas/rest/api"

# Query associations for this trait
url = f"{base_url}/efoTraits/{trait}/associations"
response = requests.get(url, headers={"Content-Type": "application/json"})
associations = response.json()

# Process results
for assoc in associations.get('_embedded', {}).get('associations', []):
    variant = assoc.get('rsId')
    pvalue = assoc.get('pvalue')
    risk_allele = assoc.get('strongestAllele')
    print(f"{variant}: p={pvalue}, risk allele={risk_allele}")
```

**示例 2：获取变异信息和所有性状关联**
```python
import requests

variant = "rs7903146"
base_url = "https://www.ebi.ac.uk/gwas/rest/api"

# Get variant details
url = f"{base_url}/singleNucleotidePolymorphisms/{variant}"
response = requests.get(url, headers={"Content-Type": "application/json"})
variant_data = response.json()

# Get all associations for this variant
url = f"{base_url}/singleNucleotidePolymorphisms/{variant}/associations"
params = {"projection": "associationBySnp"}
response = requests.get(url, params=params, headers={"Content-Type": "application/json"})
associations = response.json()

# Extract trait names and p-values
for assoc in associations.get('_embedded', {}).get('associations', []):
    trait = assoc.get('efoTrait')
    pvalue = assoc.get('pvalue')
    print(f"Trait: {trait}, p-value: {pvalue}")
```

**示例 3：访问汇总统计数据**
```python
import requests

# Query summary statistics API
base_url = "https://www.ebi.ac.uk/gwas/summary-statistics/api"

# Find associations by trait with p-value threshold
trait = "EFO_0001360"  # Type 2 diabetes
p_upper = "0.000000001"  # p < 1e-9
url = f"{base_url}/traits/{trait}/associations"
params = {
    "p_upper": p_upper,
    "size": 100  # Number of results
}
response = requests.get(url, params=params)
results = response.json()

# Process genome-wide significant hits
for hit in results.get('_embedded', {}).get('associations', []):
    variant_id = hit.get('variant_id')
    chromosome = hit.get('chromosome')
    position = hit.get('base_pair_location')
    pvalue = hit.get('p_value')
    print(f"{chromosome}:{position} ({variant_id}): p={pvalue}")
```

**示例4：按染色体区域查询**
```python
import requests

# Find variants in a specific genomic region
chromosome = "10"
start_pos = 114000000
end_pos = 115000000

base_url = "https://www.ebi.ac.uk/gwas/rest/api"
url = f"{base_url}/singleNucleotidePolymorphisms/search/findByChromBpLocationRange"
params = {
    "chrom": chromosome,
    "bpStart": start_pos,
    "bpEnd": end_pos
}
response = requests.get(url, params=params, headers={"Content-Type": "application/json"})
variants_in_region = response.json()
```

### 5. 使用汇总统计数据

GWAS 目录包含许多研究的完整摘要统计数据，提供对所有测试变体（不仅仅是全基因组显着命中）的访问。

**访问方法：**
1. **FTP下载**：http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/
2. **REST API**：基于查询访问汇总统计信息
3. **Web界面**：通过网站浏览、下载

**摘要统计 API 功能：**
- 按染色体、位置、p 值过滤
- 查询研究中的特定变体
- 检索效应大小和等位基因频率
- 访问统一和标准化的数据

**示例：下载研究的摘要统计数据**
```python
import requests
import gzip

# Get available summary statistics
base_url = "https://www.ebi.ac.uk/gwas/summary-statistics/api"
url = f"{base_url}/studies/GCST001234"
response = requests.get(url)
study_info = response.json()

# Download link is provided in the response
# Alternatively, use FTP:
# ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCSTXXXXXX/
```

### 6. 数据集成和交叉引用

GWAS 目录提供外部资源的链接：

**基因组数据库：**
- Ensembl：基因注释和变异后果
- dbSNP：变异标识符和群体频率
- gnomAD：群体等位基因频率

**功能资源：**
- 开放目标：目标-疾病关联
- PGS 目录：多基因风险评分
- UCSC 基因组浏览器：基因组背景

**表型资源：**
- EFO（实验因素本体论）：标准化特征术语
- OMIM：疾病基因关系
- 疾病本体论：疾病层次结构

**API 响应中的以下链接：**
```python
import requests

# API responses include _links for related resources
response = requests.get("https://www.ebi.ac.uk/gwas/rest/api/studies/GCST001234")
study = response.json()

# Follow link to associations
associations_url = study['_links']['associations']['href']
associations_response = requests.get(associations_url)
```

## 查询工作流程

### 工作流程 1：探索疾病的遗传关联

1. **使用 EFO 术语或自由文本识别特征**：
   - 搜索网络界面查找疾病名称
   - 记下 EFO ID（例如，EFO_0001360 表示 2 型糖尿病）

2. **通过API查询关联：**
   ```python
   url = f"https://www.ebi.ac.uk/gwas/rest/api/efoTraits/{efo_id}/associations"
   ```

3. **按重要性和总体过滤：**
   - 检查 p 值（全基因组显着性：p ≤ 5×10⁻⁸）
   - 审查研究元数据中的祖先信息
   - 按样本大小或发现/复制状态过滤

4. **提取变体详细信息：**
   - 每个协会的 rs ID
   - 效应等位基因和方向
   - 效应大小（比值比、β 系数）
   - 人群等位基因频率

5. **与其他数据库的交叉引用：**
   - 在 Ensebl 中查找变量结果
   - 检查 gnomAD 中的人口频率
   - 探索基因功能和途径

### 工作流程 2：研究特定遗传变异

1. **查询变体：**
   ```python
   url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rs_id}"
   ```

2. **检索所有特征关联：**
   ```python
   url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rs_id}/associations"
   ```

3. **分析多效性：**
   - 识别与该变体相关的所有特征
   - 审查跨特征的效果方向
   - 寻找共享的生物途径

4. **检查基因组背景：**
   - 确定附近的基因
   - 识别变体是否位于编码/监管区域
   - 审查与其他变体的连锁不平衡

### 工作流程 3：以基因为中心的关联分析

1. **在网页界面中通过基因符号搜索**或：
   ```python
   url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/search/findByGene"
   params = {"geneName": gene_symbol}
   ```

2. **检索基因区域的变异：**
   - 获取基因的染色体坐标
   - 查询区域内的变体
   - 包括启动子和监管区域（扩展边界）

3. **分析关联模式：**
   - 识别与该基因变异相关的性状
   - 寻找跨研究的一致关联
   - 检查效果大小和方向

4. **功能解释：**
   - 确定变异后果（错义、监管等）
   - 检查表达QTL（eQTL）数据
   - 审查途径和网络环境

### 工作流程 4：遗传证据的系统审查

1. **定义研究问题：**
   - 感兴趣的特定特征或疾病
   - 人口因素
   - 研究设计要求

2. **全面的变异提取：**
   - 查询特征的所有关联
   - 设置显着性阈值
   - 笔记发现和复制研究

3. **质量评估：**
   - 审查研究样本量
   - 检查人口多样性
   - 评估研究之间的异质性
   - 识别潜在的偏见

4. **数据综合：**
   - 跨研究的聚合关联
   - 如果适用，进行荟萃分析
   - 创建汇总表
   - 生成曼哈顿或森林地块

5. **导出和文档：**
   - 下载完整的协会数据
   - 如果需要导出汇总统计数据
   - 文档检索策略和日期
   - 创建可重复的分析脚本

### 工作流程 5：访问和分析摘要统计数据

1. **通过汇总统计数据确定研究：**
   - 浏览摘要统计门户
   - 检查 FTP 目录列表
   - 查询可用研究的 API

2. **下载汇总统计数据：**
   ```bash
   # Via FTP
   wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCSTXXXXXX/harmonised/GCSTXXXXXX-harmonised.tsv.gz
   ```

3. **通过API查询特定变体：**
   ```python
   url = f"https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/{chrom}/associations"
   params = {"start": start_pos, "end": end_pos}
   ```

4. **处理和分析：**
   - 按 p 值阈值过滤
   - 提取效应大小和置信区间
   - 执行下游分析（精细映射、共定位等）

## 响应格式和数据字段

**关联记录中的关键字段：**
- `rsId`：变体标识符（rs 编号）
- `strongestAllele`：关联的风险等位基因
- `pvalue`：关联 p 值
- `pvalueText`：文本形式的 P 值（可能包括不等式）
- `orPerCopyNum`：优势比或 beta 系数
- `betaNum`：效应大小（针对数量性状）
- `betaUnit`：beta 的测量单位
- `range`：置信区间
- `efoTrait`：关联特征名称
- `mappedLabel`：EFO 映射特征术语

**研究元数据字段：**
- `accessionId`：GCST 研究标识符
- `pubmedId`：PubMed ID
- `author`：第一作者
- `publicationDate`：发布日期
- `ancestryInitial`：发现人群血统
- `ancestryReplication`：复制群体祖先
- `sampleSize`：样本总大小

**分页：**
结果分页（默认每页 20 个项目）。使用以下方式导航：
- `size` 参数：每页结果数
- `page` 参数：页码（0 索引）
- `_links` 响应：下一页/上一页的 URL

## 最佳实践

### 查询策略
- 从网络界面开始识别相关的 EFO 术语并研究添加物
- 使用API进行批量数据提取和自动分析
- 对大型结果集实施分页处理
- 缓存 API 响应以最大程度地减少冗余请求

### 数据解读
- 始终检查 p 值阈值（全基因组：5×10⁻⁸）
- 审查祖先信息的人群适用性
- 评估证据强度时考虑样本量
- 检查独立研究中的重复情况
- 注意效应大小估计中的赢家诅咒

### 速率限制和道德
- 遵守API使用指南（无过多请求）
- 使用摘要统计数据下载进行全基因组分析
- 在 API 调用之间实施适当的延迟
- 执行迭代分析时在本地缓存结果
- 在出版物中引用 GWAS 目录

### 数据质量注意事项
- GWAS 目录整理已发布的关联（可能包含不一致之处）
- 已公布的效应量报告（可能需要协调）
- 一些研究报告有条件或联合关联
- 合并结果时检查研究重叠
- 注意确定和选择偏差

## Python 集成示例

查询和分析 GWAS 数据的完整工作流程：

```python
import requests
import pandas as pd
from time import sleep

def query_gwas_catalog(trait_id, p_threshold=5e-8):
    """
    Query GWAS Catalog for trait associations

    Args:
        trait_id: EFO trait identifier (e.g., 'EFO_0001360')
        p_threshold: P-value threshold for filtering

    Returns:
        pandas DataFrame with association results
    """
    base_url = "https://www.ebi.ac.uk/gwas/rest/api"
    url = f"{base_url}/efoTraits/{trait_id}/associations"

    headers = {"Content-Type": "application/json"}
    results = []
    page = 0

    while True:
        params = {"page": page, "size": 100}
        response = requests.get(url, params=params, headers=headers)

        if response.status_code != 200:
            break

        data = response.json()
        associations = data.get('_embedded', {}).get('associations', [])

        if not associations:
            break

        for assoc in associations:
            pvalue = assoc.get('pvalue')
            if pvalue and float(pvalue) <= p_threshold:
                results.append({
                    'variant': assoc.get('rsId'),
                    'pvalue': pvalue,
                    'risk_allele': assoc.get('strongestAllele'),
                    'or_beta': assoc.get('orPerCopyNum') or assoc.get('betaNum'),
                    'trait': assoc.get('efoTrait'),
                    'pubmed_id': assoc.get('pubmedId')
                })

        page += 1
        sleep(0.1)  # Rate limiting

    return pd.DataFrame(results)

# Example usage
df = query_gwas_catalog('EFO_0001360')  # Type 2 diabetes
print(df.head())
print(f"\nTotal associations: {len(df)}")
print(f"Unique variants: {df['variant'].nunique()}")
```

## 资源

### 参考文献/api_reference.md

全面的 API 文档包括：
- 两个 API 的详细端点规范
- 查询参数和过滤器的完整列表
- 响应格式规范和字段说明
- 高级查询示例和模式
- 错误处理和故障排除
- 与外部数据库集成

在以下情况下请查阅此参考资料：
- 构建复杂的 API 查询
- 了解响应结构
- 实现分页或批量操作
- 排查 API 错误
- 探索高级过滤选项

### 培训材料

GWAS 目录团队提供研讨会材料：
- GitHub 存储库：https://github.com/EBISPOT/GWAS_Catalog-workshop
- 带有示例查询的 Jupyter 笔记本
- 用于云执行的 Google Colab 集成

## 重要提示

### 数据更新
- GWAS 目录定期更新新出版物
- 定期重新运行查询以实现全面覆盖
- 随着研究发布数据添加摘要统计数据
- EFO 映射可能会随着时间的推移而更新

### 引文要求
使用 GWAS 目录数据时，引用：
- 索利斯 E 等人。 (2023) NHGRI-EBI GWAS 目录：知识库和沉积资源。核酸研究。电话号码：37953337
- 包括访问日期和版本（如果可用）
- 在讨论具体发现时引用原始研究

### 限制
- 并非所有 GWAS 出版物均包含在内（适用管理标准）
- 可用于研究子集的完整摘要统计数据
- 效应大小可能需要跨研究协调
- 人口多样性正在增长，但历史上有限
- 一些关联代表条件或联合效应

### 数据访问
- 网页界面：免费，无需注册
- REST API：免费，无需 API 密钥
- FTP 下载：开放访问
- 速率限制适用于 API（请尊重）

## 其他资源

- **GWAS 目录网站**：https://www.ebi.ac.uk/gwas/
- **文档**：https://www.ebi.ac.uk/gwas/docs
- **API 文档**：https://www.ebi.ac.uk/gwas/rest/docs/api
- **摘要统计API**：https://www.ebi.ac.uk/gwas/summary-statistics/docs/
- **FTP 站点**：http://ftp.ebi.ac.uk/pub/databases/gwas/
- **培训材料**：https://github.com/EBISPOT/GWAS_Catalog-workshop
- **PGS 目录**（多基因分数）：https://www.pgscatalog.org/
- **帮助和支持**：gwas-info@ebi.ac.uk