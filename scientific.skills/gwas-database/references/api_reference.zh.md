<!-- 此文件由机器翻译自 api_reference.md -->

# GWAS 目录 API 参考

GWAS Catalog REST API 的综合参考，包括端点规范、查询参数、响应格式和高级使用模式。

## 目录

- [API 概述](#api-overview)
- [身份验证和速率限制](#authentication-and-rate-limiting)
- [GWAS 目录 REST API](#gwas-catalog-rest-api)
- [摘要统计API](#summary-statistics-api)
- [响应格式](#response-formats)
- [错误处理](#error-handling)
- [高级查询模式](#advanced-query-patterns)
- [集成示例](#integration-examples)

## API 概述

GWAS 目录提供了两个互补的 REST API：

1. **GWAS Catalog REST API**：访问精选的 SNP 性状关联、研究和元数据
2. **摘要统计 API**：访问完整的 GWAS 摘要统计（所有测试的变体）

这两个 API 都使用 RESTful 设计原则，并采用 HAL（超文本应用程序语言）格式的 JSON 响应，其中包括用于资源导航的 `_links`。

### 基本 URL

```
GWAS Catalog API:         https://www.ebi.ac.uk/gwas/rest/api
Summary Statistics API:   https://www.ebi.ac.uk/gwas/summary-statistics/api
```

### 版本信息

GWAS Catalog REST API v2.0 于 2024 年发布，具有重大改进：
- 新的终点（出版物、基因、基因组背景、祖先）
- 增强数据公开（群组、背景特征、许可证）
- 改进的查询功能
- 更好的性能和文档

为了向后兼容，之前的 API 版本在 2026 年 5 月之前仍然可用。

## 身份验证和速率限制

### 身份验证

**无需身份验证** - 两个 API 都是开放访问，不需要 API 密钥或注册。

### 速率限制

虽然没有记录明确的速率限制，但请遵循最佳实践：
- 在连续请求之间实现延迟（例如 0.1-0.5 秒）
- 对大型结果集使用分页
- 在本地缓存响应
- 使用批量下载 (FTP) 获取全基因组数据
- 避免用快速连续的请求来打击 API

**速率限制示例：**
<<<代码块_1>>>

## GWAS 目录 REST API

主要 API 提供对策划的 GWAS 关联、研究、变异和特征的访问。

### 核心端点

#### 1. 研究

**获取所有研究：**
<<<代码块_2>>>

**获取具体研究：**
<<<代码块_3>>>

**搜索研究：**
<<<代码块_4>>>

**查询参数：**
- `page`：页码（0 索引）
- `size`：每页结果（默认值：20）
- `sort`：排序字段（例如，`publicationDate,desc`）

**示例：**
<<<代码块_5>>>

**响应字段：**
- `accessionId`：研究标识符 (GCST ID)
- `title`：研究标题
- `publicationInfo`：发布详细信息，包括 PMID
- `initialSampleSize`：发现队列描述
- `replicationSampleSize`：复制队列描述
- `ancestries`：人口血统信息
- `genotypingTechnologies`：阵列或测序平台
- `_links`：相关资源的链接

#### 2. 协会

**获取所有关联：**
<<<代码块_6>>>

**获取具体关联：**
```
GET /associations/{associationId}
```

**获取特征的关联：**
```
GET /efoTraits/{efoId}/associations
```

**获取变体的关联：**
```
GET /singleNucleotidePolymorphisms/{rsId}/associations
```

**查询参数：**
- `projection`：响应投影（例如，`associationBySnp`）
- `page`、`size`、`sort`：分页控件

**示例：**
```python
import requests

# Find all associations for type 2 diabetes
trait_id = "EFO_0001360"
url = f"https://www.ebi.ac.uk/gwas/rest/api/efoTraits/{trait_id}/associations"
params = {"size": 100, "page": 0}
response = requests.get(url, params=params, headers={"Content-Type": "application/json"})
data = response.json()

associations = data.get('_embedded', {}).get('associations', [])
print(f"Found {len(associations)} associations")
```

**响应字段：**
- `rsId`：变体标识符
- `strongestAllele`：风险或影响等位基因
- `pvalue`：关联 p 值
- `pvalueText`：报告的 P 值（可能包括不等式）
- `pvalueMantissa`：p 值的尾数
- `pvalueExponent`：p 值的指数
- `orPerCopyNum`：每个等位基因拷贝的优势比
- `betaNum`：效应大小（数量性状）
- `betaUnit`：测量单位
- `range`：置信区间
- `standardError`：标准错误
- `efoTrait`：特征名称
- `mappedLabel`：EFO 标准化术语
- `studyId`：相关研究加入

#### 3. 变体（单核苷酸多态性）

**获取变体详细信息：**
```
GET /singleNucleotidePolymorphisms/{rsId}
```

**搜索变体：**
```
GET /singleNucleotidePolymorphisms/search/findByRsId?rsId={rsId}
GET /singleNucleotidePolymorphisms/search/findByChromBpLocationRange?chrom={chr}&bpStart={start}&bpEnd={end}
GET /singleNucleotidePolymorphisms/search/findByGene?geneName={gene}
```

**示例：**
```python
import requests

# Get variant information
rs_id = "rs7903146"
url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rs_id}"
response = requests.get(url, headers={"Content-Type": "application/json"})
variant = response.json()

print(f"rsID: {variant.get('rsId')}")
print(f"Location: chr{variant.get('locations', [{}])[0].get('chromosomeName')}:{variant.get('locations', [{}])[0].get('chromosomePosition')}")
```

**响应字段：**
- `rsId`: rs 号码
- `merged`：指示变体是否与另一个合并
- `functionalClass`：变体结果
- `locations`：基因组位置数组
  - `chromosomeName`：染色体编号
- `chromosomePosition`：碱基对位置
  - `region`：基因组区域信息
- `genomicContexts`：附近的基因
- `lastUpdateDate`：上次修改日期

#### 4. 特征（EFO 术语）

**获取特征信息：**
```
GET /efoTraits/{efoId}
```

**搜索特征：**
```
GET /efoTraits/search/findByEfoUri?uri={efoUri}
GET /efoTraits/search/findByTraitIgnoreCase?trait={traitName}
```

**示例：**
```python
import requests

# Get trait details
trait_id = "EFO_0001360"
url = f"https://www.ebi.ac.uk/gwas/rest/api/efoTraits/{trait_id}"
response = requests.get(url, headers={"Content-Type": "application/json"})
trait = response.json()

print(f"Trait: {trait.get('trait')}")
print(f"EFO URI: {trait.get('uri')}")
```

#### 5. 出版物

**获取出版信息：**
```
GET /publications
GET /publications/{publicationId}
GET /publications/search/findByPubmedId?pubmedId={pmid}
```

#### 6. 基因

**获取基因信息：**
```
GET /genes
GET /genes/{geneId}
GET /genes/search/findByGeneName?geneName={symbol}
```

### 分页和导航

所有列表端点都支持分页：

```python
import requests

def get_all_associations(trait_id):
    """Retrieve all associations for a trait with pagination"""
    base_url = "https://www.ebi.ac.uk/gwas/rest/api"
    url = f"{base_url}/efoTraits/{trait_id}/associations"
    all_associations = []
    page = 0

    while True:
        params = {"page": page, "size": 100}
        response = requests.get(url, params=params, headers={"Content-Type": "application/json"})

        if response.status_code != 200:
            break

        data = response.json()
        associations = data.get('_embedded', {}).get('associations', [])

        if not associations:
            break

        all_associations.extend(associations)
        page += 1

    return all_associations
```

### HAL 链接

响应包括用于资源导航的 `_links`：

```python
import requests

# Get study and follow links to associations
response = requests.get("https://www.ebi.ac.uk/gwas/rest/api/studies/GCST001795")
study = response.json()

# Follow link to associations
associations_url = study['_links']['associations']['href']
associations_response = requests.get(associations_url)
associations = associations_response.json()
```

## 摘要统计API

获取已存入完整数据的研究的完整 GWAS 摘要统计数据。

### 基本网址
```
https://www.ebi.ac.uk/gwas/summary-statistics/api
```

### 核心端点

#### 1. 研究

**获取所有研究的摘要统计数据：**
```
GET /studies
```

**获取具体研究：**
```
GET /studies/{gcstId}
```

#### 2. 特质

**获取特征信息：**
```
GET /traits/{efoId}
```

**获取特征的关联：**
```
GET /traits/{efoId}/associations
```

**查询参数：**
- `p_lower`：p 值阈值下限
- `p_upper`：p 值阈值上限
- `size`：结果数
- `page`：页码

**示例：**
```python
import requests

# Find highly significant associations for a trait
trait_id = "EFO_0001360"
base_url = "https://www.ebi.ac.uk/gwas/summary-statistics/api"
url = f"{base_url}/traits/{trait_id}/associations"
params = {
    "p_upper": "0.000000001",  # p < 1e-9
    "size": 100
}
response = requests.get(url, params=params)
results = response.json()
```

#### 3. 染色体

**通过染色体获取关联：**
```
GET /chromosomes/{chromosome}/associations
```

**按基因组区域查询：**
```
GET /chromosomes/{chromosome}/associations?start={start}&end={end}
```

**示例：**
```python
import requests

# Query variants in a specific region
chromosome = "10"
start_pos = 114000000
end_pos = 115000000

base_url = "https://www.ebi.ac.uk/gwas/summary-statistics/api"
url = f"{base_url}/chromosomes/{chromosome}/associations"
params = {
    "start": start_pos,
    "end": end_pos,
    "size": 1000
}
response = requests.get(url, params=params)
variants = response.json()
```

#### 4. 变体

**获取跨研究的特定变体：**
```
GET /variants/{variantId}
```

**按变体 ID 搜索：**
```
GET /variants/{variantId}/associations
```

### 响应字段

**关联领域：**
- `variant_id`：变体标识符
- `chromosome`：染色体编号
- `base_pair_location`：位置（bp）
- `effect_allele`：效果等位基因
- `other_allele`：参考等位基因
- `effect_allele_frequency`：等位基因频率
- `beta`：效果大小
- `standard_error`：标准错误
- `p_value`：P 值
- `ci_lower`：较低的置信区间
- `ci_upper`：置信区间上限
- `odds_ratio`：优势比（病例对照研究）
- `study_accession`：GCST ID

## 响应格式

### 内容类型

所有 API 请求都应包含标头：
```
Content-Type: application/json
```

### HAL 格式

响应遵循 HAL（超文本应用程序语言）规范：

```json
{
  "_embedded": {
    "associations": [
      {
        "rsId": "rs7903146",
        "pvalue": 1.2e-30,
        "efoTrait": "type 2 diabetes",
        "_links": {
          "self": {
            "href": "https://www.ebi.ac.uk/gwas/rest/api/associations/12345"
          }
        }
      }
    ]
  },
  "_links": {
    "self": {
      "href": "https://www.ebi.ac.uk/gwas/rest/api/efoTraits/EFO_0001360/associations?page=0"
    },
    "next": {
      "href": "https://www.ebi.ac.uk/gwas/rest/api/efoTraits/EFO_0001360/associations?page=1"
    }
  },
  "page": {
    "size": 20,
    "totalElements": 1523,
    "totalPages": 77,
    "number": 0
  }
}
```

### 页面元数据

分页响应包括页面信息：
- `size`：每页的项目数
- `totalElements`：结果总数
- `totalPages`：总页数
- `number`：当前页码（0索引）

## 错误处理

### HTTP 状态代码

- `200 OK`：请求成功
- `400 Bad Request`：无效参数
- `404 Not Found`：找不到资源
- `500 Internal Server Error`：服务器错误

### 错误响应格式

```json
{
  "timestamp": "2025-10-19T12:00:00.000+00:00",
  "status": 404,
  "error": "Not Found",
  "message": "No association found with id: 12345",
  "path": "/gwas/rest/api/associations/12345"
}
```

### 错误处理示例

```python
import requests

def safe_api_request(url, params=None):
    """Make API request with error handling"""
    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.HTTPError as e:
        print(f"HTTP Error: {e}")
        print(f"Response: {response.text}")
        return None
    except requests.exceptions.ConnectionError:
        print("Connection error - check network")
        return None
    except requests.exceptions.Timeout:
        print("Request timed out")
        return None
    except requests.exceptions.RequestException as e:
        print(f"Request error: {e}")
        return None
```

## 高级查询模式

### 1.交叉引用变体和性状

```python
import requests

def get_variant_pleiotropy(rs_id):
    """Get all traits associated with a variant"""
    base_url = "https://www.ebi.ac.uk/gwas/rest/api"
    url = f"{base_url}/singleNucleotidePolymorphisms/{rs_id}/associations"
    params = {"projection": "associationBySnp"}

    response = requests.get(url, params=params, headers={"Content-Type": "application/json"})
    data = response.json()

    traits = {}
    for assoc in data.get('_embedded', {}).get('associations', []):
        trait = assoc.get('efoTrait')
        pvalue = assoc.get('pvalue')
        if trait:
            if trait not in traits or float(pvalue) < float(traits[trait]):
                traits[trait] = pvalue

    return traits

# Example usage
pleiotropy = get_variant_pleiotropy('rs7903146')
for trait, pval in sorted(pleiotropy.items(), key=lambda x: float(x[1])):
    print(f"{trait}: p={pval}")
```

### 2. 按 P 值阈值过滤

```python
import requests

def get_significant_associations(trait_id, p_threshold=5e-8):
    """Get genome-wide significant associations"""
    base_url = "https://www.ebi.ac.uk/gwas/rest/api"
    url = f"{base_url}/efoTraits/{trait_id}/associations"

    results = []
    page = 0

    while True:
        params = {"page": page, "size": 100}
        response = requests.get(url, params=params, headers={"Content-Type": "application/json"})

        if response.status_code != 200:
            break

        data = response.json()
        associations = data.get('_embedded', {}).get('associations', [])

        if not associations:
            break

        for assoc in associations:
            pvalue = assoc.get('pvalue')
            if pvalue and float(pvalue) <= p_threshold:
                results.append(assoc)

        page += 1

    return results
```

### 3. 组合主要统计 API 和摘要统计 API

```python
import requests

def get_complete_variant_data(rs_id):
    """Get variant data from both APIs"""
    main_url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rs_id}"

    # Get basic variant info
    response = requests.get(main_url, headers={"Content-Type": "application/json"})
    variant_info = response.json()

    # Get associations
    assoc_url = f"{main_url}/associations"
    response = requests.get(assoc_url, headers={"Content-Type": "application/json"})
    associations = response.json()

    # Could also query summary statistics API for this variant
    # across all studies with summary data

    return {
        "variant": variant_info,
        "associations": associations
    }
```

### 4. 基因组区域查询

```python
import requests

def query_region(chromosome, start, end, p_threshold=None):
    """Query variants in genomic region"""
    # From main API
    base_url = "https://www.ebi.ac.uk/gwas/rest/api"
    url = f"{base_url}/singleNucleotidePolymorphisms/search/findByChromBpLocationRange"
    params = {
        "chrom": chromosome,
        "bpStart": start,
        "bpEnd": end,
        "size": 1000
    }

    response = requests.get(url, params=params, headers={"Content-Type": "application/json"})
    variants = response.json()

    # Can also query summary statistics API
    sumstats_url = f"https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/{chromosome}/associations"
    sumstats_params = {"start": start, "end": end, "size": 1000}
    if p_threshold:
        sumstats_params["p_upper"] = str(p_threshold)

    sumstats_response = requests.get(sumstats_url, params=sumstats_params)
    sumstats = sumstats_response.json()

    return {
        "catalog_variants": variants,
        "summary_stats": sumstats
    }
```

## 集成示例

### 完整工作流程：疾病遗传架构

```python
import requests
import pandas as pd
from time import sleep

class GWASCatalogQuery:
    def __init__(self):
        self.base_url = "https://www.ebi.ac.uk/gwas/rest/api"
        self.headers = {"Content-Type": "application/json"}

    def get_trait_associations(self, trait_id, p_threshold=5e-8):
        """Get all associations for a trait"""
        url = f"{self.base_url}/efoTraits/{trait_id}/associations"
        results = []
        page = 0

        while True:
            params = {"page": page, "size": 100}
            response = requests.get(url, params=params, headers=self.headers)

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
                        'rs_id': assoc.get('rsId'),
                        'pvalue': float(pvalue),
                        'risk_allele': assoc.get('strongestAllele'),
                        'or_beta': assoc.get('orPerCopyNum') or assoc.get('betaNum'),
                        'study': assoc.get('studyId'),
                        'pubmed_id': assoc.get('pubmedId')
                    })

            page += 1
            sleep(0.1)

        return pd.DataFrame(results)

    def get_variant_details(self, rs_id):
        """Get detailed variant information"""
        url = f"{self.base_url}/singleNucleotidePolymorphisms/{rs_id}"
        response = requests.get(url, headers=self.headers)

        if response.status_code == 200:
            return response.json()
        return None

    def get_gene_associations(self, gene_name):
        """Get variants associated with a gene"""
        url = f"{self.base_url}/singleNucleotidePolymorphisms/search/findByGene"
        params = {"geneName": gene_name}
        response = requests.get(url, params=params, headers=self.headers)

        if response.status_code == 200:
            return response.json()
        return None

# Example usage
gwas = GWASCatalogQuery()

# Query type 2 diabetes associations
df = gwas.get_trait_associations('EFO_0001360')
print(f"Found {len(df)} genome-wide significant associations")
print(f"Unique variants: {df['rs_id'].nunique()}")

# Get top variants
top_variants = df.nsmallest(10, 'pvalue')
print("\nTop 10 variants:")
print(top_variants[['rs_id', 'pvalue', 'risk_allele']])

# Get details for top variant
if len(top_variants) > 0:
    top_rs = top_variants.iloc[0]['rs_id']
    variant_info = gwas.get_variant_details(top_rs)
    if variant_info:
        loc = variant_info.get('locations', [{}])[0]
        print(f"\n{top_rs} location: chr{loc.get('chromosomeName')}:{loc.get('chromosomePosition')}")
```

### FTP 下载集成

```python
import requests
from pathlib import Path

def download_summary_statistics(gcst_id, output_dir="."):
    """Download summary statistics from FTP"""
    # FTP URL pattern
    ftp_base = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"

    # Try harmonised file first
    harmonised_url = f"{ftp_base}/{gcst_id}/harmonised/{gcst_id}-harmonised.tsv.gz"

    output_path = Path(output_dir) / f"{gcst_id}.tsv.gz"

    try:
        response = requests.get(harmonised_url, stream=True)
        response.raise_for_status()

        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

        print(f"Downloaded {gcst_id} to {output_path}")
        return output_path

    except requests.exceptions.HTTPError:
        print(f"Harmonised file not found for {gcst_id}")
        return None

# Example usage
download_summary_statistics("GCST001234", output_dir="./sumstats")
```

## 其他资源

- **交互式 API 文档**：https://www.ebi.ac.uk/gwas/rest/docs/api
- **统计 API 文档摘要**：https://www.ebi.ac.uk/gwas/summary-statistics/docs/
- **研讨会材料**：https://github.com/EBISPOT/GWAS_Catalog-workshop
- **关于 API v2 的博客文章**：https://ebispot.github.io/gwas-blog/rest-api-v2-release/
- **R 包 (gwasrapidd)**: https://cran.r-project.org/package=gwasrapidd