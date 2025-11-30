<!-- 此文件由机器翻译自 api_reference.md -->

# Reactome API 参考

本文档提供了 Reactome 的 REST API 的全面参考信息。

## 基本 URL

- **内容服务**：`https://reactome.org/ContentService`
- **分析服务**：`https://reactome.org/AnalysisService`

## 内容服务API

内容服务通过 REST 端点提供对 Reactome 精选路径数据的访问。

### 数据库信息

#### 获取数据库版本
```
GET /data/database/version
```

**响应：** 包含数据库版本号的纯文本

**示例：**
<<<代码块_1>>>

#### 获取数据库名称
<<<代码块_2>>>

**响应：** 包含数据库名称的纯文本

### 实体查询

#### 通过ID查询实体
<<<代码块_3>>>

**参数：**
- `id`（路径）：稳定标识符或数据库 ID（例如“R-HSA-69278”）

**响应：** JSON 对象包含完整的实体信息，包括：
- `stId`：稳定标识符
- `displayName`：人类可读的名称
- `schemaClass`：实体类型（路径、反应、复杂等）
- `species`：物种信息数组
- 附加类型特定字段

**示例：**
<<<代码块_4>>>

#### 查询实体属性
<<<代码块_5>>>

**参数：**
- `id`（路径）：实体标识符
- `attribute`（路径）：特定属性名称（例如，“displayName”、“compartment”）

**响应：** JSON 或纯文本，具体取决于属性类型

**示例：**
<<<代码块_6>>>

### 路径查询

#### 获取路径实体
```
GET /data/event/{id}/participatingPhysicalEntities
```

**参数：**
- `id`（路径）：路径或反应稳定标识符

**响应：** 参与该路径的物理实体（蛋白质、复合物、小分子）的 JSON 数组

**示例：**
```python
response = requests.get(
    "https://reactome.org/ContentService/data/event/R-HSA-69278/participatingPhysicalEntities"
)
entities = response.json()
for entity in entities:
    print(f"{entity['stId']}: {entity['displayName']} ({entity['schemaClass']})")
```

#### 收容事件
```
GET /data/pathway/{id}/containedEvents
```

**参数：**
- `id`（路径）：路径稳定标识符

**响应：** 路径中包含的事件（反应、子路径）的 JSON 数组

### 搜索查询

#### 按名称搜索
```
GET /data/query?name={query}
```

**参数：**
- `name`（查询）：搜索词

**响应：** 匹配实体的 JSON 数组

**示例：**
```python
response = requests.get(
    "https://reactome.org/ContentService/data/query",
    params={"name": "glycolysis"}
)
results = response.json()
```

## 分析服务API

分析服务执行通路富集和表达分析。

### 提交分析

#### 提交标识符（POST）
```
POST /identifiers/
POST /identifiers/projection/  # Map to human pathways only
```

**标题：**
- `Content-Type: text/plain`

**身体：**
- 对于过度代表性：标识符的纯文本列表（每行一个）
- 用于表达式分析：TSV 格式，标题以“#”开头

**表达数据格式：**
```
#Gene	Sample1	Sample2	Sample3
TP53	2.5	3.1	2.8
BRCA1	1.2	1.5	1.3
```

**响应：** JSON 对象包含：
```json
{
  "summary": {
    "token": "MzUxODM3NTQzMDAwMDA1ODI4MA==",
    "type": "OVERREPRESENTATION",
    "species": "9606",
    "sampleName": null,
    "fileName": null,
    "text": true
  },
  "pathways": [
    {
      "stId": "R-HSA-69278",
      "name": "Cell Cycle, Mitotic",
      "species": {
        "name": "Homo sapiens",
        "taxId": "9606"
      },
      "entities": {
        "found": 15,
        "total": 450,
        "pValue": 0.0000234,
        "fdr": 0.00156
      },
      "reactions": {
        "found": 12,
        "total": 342
      }
    }
  ],
  "resourceSummary": [
    {
      "resource": "TOTAL",
      "pathways": 25
    }
  ]
}
```

**示例：**
```python
import requests

# Overrepresentation analysis
identifiers = ["TP53", "BRCA1", "EGFR", "MYC", "CDK1"]
data = "\n".join(identifiers)

response = requests.post(
    "https://reactome.org/AnalysisService/identifiers/",
    headers={"Content-Type": "text/plain"},
    data=data
)

result = response.json()
token = result["summary"]["token"]

# Process pathways
for pathway in result["pathways"]:
    print(f"Pathway: {pathway['name']}")
    print(f"  Found: {pathway['entities']['found']}/{pathway['entities']['total']}")
    print(f"  p-value: {pathway['entities']['pValue']:.6f}")
    print(f"  FDR: {pathway['entities']['fdr']:.6f}")
```

#### 提交文件（表单上传）
```
POST /identifiers/form/
```

**内容类型：** `multipart/form-data`

**参数：**
- `file`：包含标识符或表达式数据的文件

#### 提交网址
```
POST /identifiers/url/
```

**参数：**
- `url`：指向数据文件的 URL

### 检索分析结果

#### 通过令牌获取结果
```
GET /token/{token}
GET /token/{token}/projection/  # With species projection
```

**参数：**
- `token`（路径）：提交返回的分析令牌

**响应：** 与初始分析响应结构相同

**示例：**
```python
token = "MzUxODM3NTQzMDAwMDA1ODI4MA=="
response = requests.get(f"https://reactome.org/AnalysisService/token/{token}")
results = response.json()
```

**注意：** 令牌的有效期为 7 天

#### 过滤结果
```
GET /token/{token}/filter/pathways?resource={resource}
```

**参数：**
- `token`（路径）：分析标记
- `resource`（查询）：资源过滤器（例如“TOTAL”、“UNIPROT”、“ENSEMBL”）

### 下载结果

#### 下载为 CSV
```
GET /download/{token}/pathways/{resource}/result.csv
```

#### 下载地图
```
GET /download/{token}/entities/found/{resource}/mapping.tsv
```

## 支持的标识符

Reactome 自动检测并处理各种标识符类型：

### 蛋白质和基因
- **UniProt**：P04637
- **基因符号**：TP53
- **合奏**：ENSG00000141510
- **EntrezGene**：7157
- **参考序列**：NM_000546
- **OMIM**：191170

### 小分子
- **ChEBI**：CHEBI：15377
- **KEGG 化合物**：C00031
- **PubChem**：702

### 其他
- **miRBase**：hsa-miR-21
- **InterPro**：IPR011616

## 响应格式

### JSON 对象

实体对象包含标准化字段：
```json
{
  "stId": "R-HSA-69278",
  "displayName": "Cell Cycle, Mitotic",
  "schemaClass": "Pathway",
  "species": [
    {
      "dbId": 48887,
      "displayName": "Homo sapiens",
      "taxId": "9606"
    }
  ],
  "isInDisease": false
}
```

### TSV 格式

对于批量查询，TSV 返回：
```
stId	displayName	schemaClass
R-HSA-69278	Cell Cycle, Mitotic	Pathway
R-HSA-69306	DNA Replication	Pathway
```

## 错误响应

### HTTP 状态代码
- `200`：成功
- `400`：错误请求（无效参数）
- `404`：未找到（无效 ID）
- `415`：不支持的媒体类型
- `500`：内部服务器错误

### 错误 JSON 结构
```json
{
  "code": 404,
  "reason": "NOT_FOUND",
  "messages": ["Pathway R-HSA-INVALID not found"]
}
```

## 速率限制
Reactome 目前不执行严格的速率限制，但请考虑：
- 在请求之间实施合理的延迟
- 可用时使用批处理操作
- 适当时缓存结果
- 遵守 7 天令牌有效期

## 最佳实践

### 1. 使用分析令牌
存储和重用分析标记以避免冗余计算：
```python
# Store token after analysis
token = result["summary"]["token"]
save_token(token)  # Save to file or database

# Retrieve results later
result = requests.get(f"https://reactome.org/AnalysisService/token/{token}")
```

### 2.批量查询
在单个请求中提交多个标识符，而不是单个查询：
```python
# Good: Single batch request
identifiers = ["TP53", "BRCA1", "EGFR"]
result = analyze_batch(identifiers)

# Avoid: Multiple individual requests
# for gene in genes:
#     result = analyze_single(gene)  # Don't do this
```

### 3. 适当处理物种
使用 `/projection/` 端点将非人类标识符映射到人类路径：
```python
# For mouse genes, project to human pathways
response = requests.post(
    "https://reactome.org/AnalysisService/identifiers/projection/",
    headers={"Content-Type": "text/plain"},
    data=mouse_genes
)
```

### 4. 处理大型结果集
对于返回许多路径的分析，请按重要性过滤：
```python
significant_pathways = [
    p for p in result["pathways"]
    if p["entities"]["fdr"] < 0.05
]
```

## 集成示例

### 完整的分析工作流程
```python
import requests
import json

def analyze_gene_list(genes, output_file="analysis_results.json"):
    """
    Perform pathway enrichment analysis on a list of genes
    """
    # Submit analysis
    data = "\n".join(genes)
    response = requests.post(
        "https://reactome.org/AnalysisService/identifiers/",
        headers={"Content-Type": "text/plain"},
        data=data
    )

    if response.status_code != 200:
        raise Exception(f"Analysis failed: {response.text}")

    result = response.json()
    token = result["summary"]["token"]

    # Filter significant pathways (FDR < 0.05)
    significant = [
        p for p in result["pathways"]
        if p["entities"]["fdr"] < 0.05
    ]

    # Save results
    with open(output_file, "w") as f:
        json.dump({
            "token": token,
            "total_pathways": len(result["pathways"]),
            "significant_pathways": len(significant),
            "pathways": significant
        }, f, indent=2)

    # Generate browser URL for top pathway
    if significant:
        top_pathway = significant[0]
        url = f"https://reactome.org/PathwayBrowser/#{top_pathway['stId']}&DTAB=AN&ANALYSIS={token}"
        print(f"View top result: {url}")

    return result

# Usage
genes = ["TP53", "BRCA1", "BRCA2", "CDK1", "CDK2"]
result = analyze_gene_list(genes)
```

## 其他资源

- **交互式 API 文档**：https://reactome.org/dev/content-service
- **分析服务文档**：https://reactome.org/dev/analysis
- **用户指南**：https://reactome.org/userguide
- **数据下载**：https://reactome.org/download-data