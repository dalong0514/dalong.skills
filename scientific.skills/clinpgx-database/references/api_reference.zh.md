<!-- 此文件由机器翻译自 api_reference.md -->

# ClinPGx API 参考

ClinPGx REST API 的完整参考文档。

## 基本网址

```
https://api.clinpgx.org/v1/
```

## 速率限制

- **最大速率**：每秒 2 个请求
- **强制**：超过限制的请求将收到 HTTP 429（请求过多）
- **最佳实践**：在请求之间实现 500 毫秒的延迟（0.5 秒）
- **建议**：对于大量 API 使用，请联系 api@clinpgx.org

## 身份验证

基本 API 访问不需要身份验证。所有端点均可公开访问。

## 数据许可

通过 API 访问的所有数据均须遵守：
- 知识共享署名-相同方式共享 4.0 国际许可证
- ClinPGx 数据使用政策

## 响应格式

所有成功的响应都会返回带有适当 HTTP 状态代码的 JSON：
- `200 OK`：请求成功
- `404 Not Found`：资源不存在
- `429 Too Many Requests`：超出速率限制
- `500 Internal Server Error`：服务器错误

## 核心端点

### 1. 基因端点

检索药物基因信息，包括功能、变异和临床意义。

#### 通过符号获取基因

<<<代码块_1>>>

**参数：**
- `gene_symbol`（路径，必填）：基因符号（例如 CYP2D6、TPMT、DPYD）

**请求示例：**
<<<代码块_2>>>

**响应示例：**
<<<代码块_3>>>

#### 搜索基因

<<<代码块_4>>>

**参数：**
- `q`（查询，可选）：基因名称或符号的搜索词

**示例：**
<<<代码块_5>>>

### 2. 化学/药物终点

访问药物和化合物信息，包括药物基因组注释。

#### 通过 ID 获取药物

<<<代码块_6>>>

**参数：**
- `drug_id`（路径，必需）：ClinPGx 药物标识符（例如 PA448515）

**请求示例：**
```bash
curl "https://api.clinpgx.org/v1/chemical/PA448515"
```

#### 按名称搜索药物

```http
GET /v1/chemical?name={drug_name}
```

**参数：**
- `name`（查询，可选）：药物名称或同义词

**示例：**
```bash
curl "https://api.clinpgx.org/v1/chemical?name=warfarin"
```

**响应示例：**
```json
[
  {
    "id": "PA448515",
    "name": "warfarin",
    "genericNames": ["warfarin sodium"],
    "tradeNames": ["Coumadin", "Jantoven"],
    "drugClasses": ["Anticoagulants"],
    "indication": "Prevention of thrombosis",
    "relatedGenes": ["CYP2C9", "VKORC1", "CYP4F2"]
  }
]
```

### 3. 基因-药物对终点

通过临床注释查询策划的基因-药物相互作用关系。

#### 获取基因-药物对

```http
GET /v1/geneDrugPair?gene={gene}&drug={drug}
```

**参数：**
- `gene`（查询，可选）：基因符号
- `drug`（查询，可选）：药品名称
- `cpicLevel`（查询，可选）：按 CPIC 推荐级别（A、B、C、D）过滤

**请求示例：**
```bash
# Get all pairs for a gene
curl "https://api.clinpgx.org/v1/geneDrugPair?gene=CYP2D6"

# Get specific gene-drug pair
curl "https://api.clinpgx.org/v1/geneDrugPair?gene=CYP2D6&drug=codeine"

# Get all CPIC Level A pairs
curl "https://api.clinpgx.org/v1/geneDrugPair?cpicLevel=A"
```

**响应示例：**
```json
[
  {
    "gene": "CYP2D6",
    "drug": "codeine",
    "sources": ["CPIC", "FDA", "DPWG"],
    "cpicLevel": "A",
    "evidenceLevel": "1A",
    "clinicalAnnotationCount": 45,
    "hasGuideline": true,
    "guidelineUrl": "https://www.clinpgx.org/guideline/..."
  }
]
```

### 4. 指导端点

获取来自 CPIC、DPWG 和其他来源的临床实践指南。

#### 获取指南

```http
GET /v1/guideline?source={source}&gene={gene}&drug={drug}
```

**参数：**
- `source`（查询，可选）：指南来源（CPIC、DPWG、FDA）
- `gene`（查询，可选）：基因符号
- `drug`（查询，可选）：药品名称

**请求示例：**
```bash
# Get all CPIC guidelines
curl "https://api.clinpgx.org/v1/guideline?source=CPIC"

# Get guideline for specific gene-drug
curl "https://api.clinpgx.org/v1/guideline?gene=CYP2C19&drug=clopidogrel"
```

#### 通过 ID 获取指南

```http
GET /v1/guideline/{guideline_id}
```

**示例：**
```bash
curl "https://api.clinpgx.org/v1/guideline/PA166104939"
```

**响应示例：**
```json
{
  "id": "PA166104939",
  "name": "CPIC Guideline for CYP2C19 and Clopidogrel",
  "source": "CPIC",
  "genes": ["CYP2C19"],
  "drugs": ["clopidogrel"],
  "recommendationLevel": "A",
  "lastUpdated": "2023-08-01",
  "summary": "Alternative antiplatelet therapy recommended for...",
  "recommendations": [...],
  "pdfUrl": "https://www.clinpgx.org/...",
  "pmid": "23400754"
}
```

### 5. 等位基因端点

查询等位基因定义、功能和群体频率。

#### 获取一个基因的所有等位基因

```http
GET /v1/allele?gene={gene_symbol}
```

**参数：**
- `gene`（查询，必需）：基因符号

**请求示例：**
```bash
curl "https://api.clinpgx.org/v1/allele?gene=CYP2D6"
```

**响应示例：**
```json
[
  {
    "name": "CYP2D6*1",
    "gene": "CYP2D6",
    "function": "Normal function",
    "activityScore": 1.0,
    "frequencies": {
      "European": 0.42,
      "African": 0.37,
      "East Asian": 0.50,
      "Latino": 0.44
    },
    "definingVariants": ["Reference allele"],
    "pharmVarId": "PV00001"
  },
  {
    "name": "CYP2D6*4",
    "gene": "CYP2D6",
    "function": "No function",
    "activityScore": 0.0,
    "frequencies": {
      "European": 0.20,
      "African": 0.05,
      "East Asian": 0.01,
      "Latino": 0.10
    },
    "definingVariants": ["rs3892097"],
    "pharmVarId": "PV00004"
  }
]
```

#### 获取特定等位基因

```http
GET /v1/allele/{allele_name}
```

**参数：**
- `allele_name`（路径，必需）：带有星号命名法的等位基因名称（例如，CYP2D6*4）

**示例：**
```bash
curl "https://api.clinpgx.org/v1/allele/CYP2D6*4"
```

### 6. 变体端点

搜索遗传变异及其药物基因组注释。

#### 通过 rsID 获取变体

```http
GET /v1/variant/{rsid}
```

**参数：**
- `rsid`（路径，必需）：dbSNP 参考 SNP ID

**请求示例：**
```bash
curl "https://api.clinpgx.org/v1/variant/rs4244285"
```

**响应示例：**
```json
{
  "rsid": "rs4244285",
  "chromosome": "10",
  "position": 94781859,
  "gene": "CYP2C19",
  "alleles": ["CYP2C19*2"],
  "consequence": "Splice site variant",
  "clinicalSignificance": "Pathogenic - reduced enzyme activity",
  "frequencies": {
    "European": 0.15,
    "African": 0.18,
    "East Asian": 0.29,
    "Latino": 0.12
  },
  "references": [...]
}
```

#### 按位置搜索变体

```http
GET /v1/variant?chromosome={chr}&position={pos}
```

**参数：**
- `chromosome`（查询，可选）：染色体编号（1-22，X，Y）
- `position`（查询，可选）：基因组位置 (GRCh38)

**示例：**
```bash
curl "https://api.clinpgx.org/v1/variant?chromosome=10&position=94781859"
```

### 7. 临床注释端点

访问基因-药物-表型关系的精选文献注释。

#### 获取临床注释

```http
GET /v1/clinicalAnnotation?gene={gene}&drug={drug}&evidenceLevel={level}
```

**参数：**
- `gene`（查询，可选）：基因符号
- `drug`（查询，可选）：药品名称
- `evidenceLevel`（查询，可选）：证据级别（1A、1B、2A、2B、3、4）
- `phenotype`（查询，可选）：表型或结果

**请求示例：**
```bash
# Get all annotations for a gene
curl "https://api.clinpgx.org/v1/clinicalAnnotation?gene=CYP2D6"

# Get high-quality evidence only
curl "https://api.clinpgx.org/v1/clinicalAnnotation?evidenceLevel=1A"

# Get annotations for specific gene-drug pair
curl "https://api.clinpgx.org/v1/clinicalAnnotation?gene=TPMT&drug=azathioprine"
```

**响应示例：**
```json
[
  {
    "id": "PA166153683",
    "gene": "CYP2D6",
    "drug": "codeine",
    "phenotype": "Reduced analgesic effect",
    "evidenceLevel": "1A",
    "annotation": "Poor metabolizers have reduced conversion...",
    "pmid": "24618998",
    "studyType": "Clinical trial",
    "population": "European",
    "sources": ["CPIC"]
  }
]
```
**证据级别：**
- **1A**：来自指南的高质量证据（CPIC、FDA、DPWG）
- **1B**：高质量证据尚未成为指南
- **2A**：来自精心设计的研究的适度证据
- **2B**：中等证据，但有一些局限性
- **3**：证据有限或相互矛盾
- **4**：病例报告或证据薄弱

### 8. 药品标签端点

检索包含药物基因组内容的监管药物标签信息。

#### 获取药品标签

```http
GET /v1/drugLabel?drug={drug_name}&source={source}
```

**参数：**
- `drug`（查询，必填）：药品名称
- `source`（查询，可选）：监管来源（FDA、EMA、PMDA、加拿大卫生部）

**请求示例：**
```bash
# Get all labels for warfarin
curl "https://api.clinpgx.org/v1/drugLabel?drug=warfarin"

# Get only FDA labels
curl "https://api.clinpgx.org/v1/drugLabel?drug=warfarin&source=FDA"
```

**响应示例：**
```json
[
  {
    "id": "DL001234",
    "drug": "warfarin",
    "source": "FDA",
    "sections": {
      "testing": "Consider CYP2C9 and VKORC1 genotyping...",
      "dosing": "Dose adjustment based on genotype...",
      "warnings": "Risk of bleeding in certain genotypes"
    },
    "biomarkers": ["CYP2C9", "VKORC1"],
    "testingRecommended": true,
    "labelUrl": "https://dailymed.nlm.nih.gov/...",
    "lastUpdated": "2024-01-15"
  }
]
```

### 9. 路径端点

访问药代动力学和药效学途径图和信息。

#### 通过 ID 获取 Pathway

```http
GET /v1/pathway/{pathway_id}
```

**参数：**
- `pathway_id`（路径，必需）：ClinPGx 通路标识符

**示例：**
```bash
curl "https://api.clinpgx.org/v1/pathway/PA146123006"
```

#### 搜索路径

```http
GET /v1/pathway?drug={drug_name}&gene={gene}
```

**参数：**
- `drug`（查询，可选）：药品名称
- `gene`（查询，可选）：基因符号

**示例：**
```bash
curl "https://api.clinpgx.org/v1/pathway?drug=warfarin"
```

**响应示例：**
```json
{
  "id": "PA146123006",
  "name": "Warfarin Pharmacokinetics and Pharmacodynamics",
  "drugs": ["warfarin"],
  "genes": ["CYP2C9", "VKORC1", "CYP4F2", "GGCX"],
  "description": "Warfarin is metabolized primarily by CYP2C9...",
  "diagramUrl": "https://www.clinpgx.org/pathway/...",
  "steps": [
    {
      "step": 1,
      "process": "Absorption",
      "genes": []
    },
    {
      "step": 2,
      "process": "Metabolism",
      "genes": ["CYP2C9", "CYP2C19"]
    },
    {
      "step": 3,
      "process": "Target interaction",
      "genes": ["VKORC1"]
    }
  ]
}
```

## 查询模式和示例

### 常见查询模式

#### 1. 患者用药审查

查询患者药物的所有基因-药物对：

```python
import requests

patient_meds = ["clopidogrel", "simvastatin", "codeine"]
patient_genes = {"CYP2C19": "*1/*2", "CYP2D6": "*1/*1", "SLCO1B1": "*1/*5"}

for med in patient_meds:
    for gene in patient_genes:
        response = requests.get(
            "https://api.clinpgx.org/v1/geneDrugPair",
            params={"gene": gene, "drug": med}
        )
        pairs = response.json()
        # Check for interactions
```

#### 2. 可操作的基因面板

查找所有具有 CPIC A 级建议的基因：

```python
response = requests.get(
    "https://api.clinpgx.org/v1/geneDrugPair",
    params={"cpicLevel": "A"}
)
actionable_pairs = response.json()

genes = set(pair['gene'] for pair in actionable_pairs)
print(f"Panel should include: {sorted(genes)}")
```

#### 3. 人口频率分析

比较不同人群的等位基因频率：

```python
alleles = requests.get(
    "https://api.clinpgx.org/v1/allele",
    params={"gene": "CYP2D6"}
).json()

# Calculate phenotype frequencies
pm_freq = {}  # Poor metabolizer frequencies
for allele in alleles:
    if allele['function'] == 'No function':
        for pop, freq in allele['frequencies'].items():
            pm_freq[pop] = pm_freq.get(pop, 0) + freq
```

#### 4. 药物安全筛查

检查高风险基因药物关联：

```python
# Screen for HLA-B*57:01 before abacavir
response = requests.get(
    "https://api.clinpgx.org/v1/geneDrugPair",
    params={"gene": "HLA-B", "drug": "abacavir"}
)
pair = response.json()[0]

if pair['cpicLevel'] == 'A':
    print("CRITICAL: Do not use if HLA-B*57:01 positive")
```

## 错误处理

### 常见错误响应

#### 404 未找到
```json
{
  "error": "Resource not found",
  "message": "Gene 'INVALID' does not exist"
}
```

#### 429 请求过多
```json
{
  "error": "Rate limit exceeded",
  "message": "Maximum 2 requests per second allowed"
}
```

### 推荐的错误处理模式

```python
import requests
import time

def safe_query(url, params=None, max_retries=3):
    for attempt in range(max_retries):
        try:
            response = requests.get(url, params=params, timeout=10)

            if response.status_code == 200:
                time.sleep(0.5)  # Rate limiting
                return response.json()
            elif response.status_code == 429:
                wait = 2 ** attempt
                print(f"Rate limited. Waiting {wait}s...")
                time.sleep(wait)
            elif response.status_code == 404:
                print("Resource not found")
                return None
            else:
                response.raise_for_status()

        except requests.RequestException as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            if attempt == max_retries - 1:
                raise

    return None
```

## 最佳实践

### 速率限制
- 在请求之间实现 500 毫秒的延迟（最多 2 个请求/秒）
- 对速率限制错误使用指数退避
- 考虑缓存经常访问的数据的结果
- 对于批量操作，请联系 api@clinpgx.org

### 缓存策略
```python
import json
from pathlib import Path

def cached_query(cache_file, query_func, *args, **kwargs):
    cache_path = Path(cache_file)

    if cache_path.exists():
        with open(cache_path) as f:
            return json.load(f)

    result = query_func(*args, **kwargs)

    if result:
        with open(cache_path, 'w') as f:
            json.dump(result, f)

    return result
```

### 批处理
```python
import time

def batch_gene_query(genes, delay=0.5):
    results = {}
    for gene in genes:
        response = requests.get(f"https://api.clinpgx.org/v1/gene/{gene}")
        if response.status_code == 200:
            results[gene] = response.json()
        time.sleep(delay)
    return results
```

## 数据模式定义

### 基因对象
```typescript
{
  id: string;              // ClinPGx gene ID
  symbol: string;          // HGNC gene symbol
  name: string;            // Full gene name
  chromosome: string;      // Chromosome location
  function: string;        // Pharmacogenomic function
  clinicalAnnotations: number;  // Count of annotations
  relatedDrugs: string[];  // Associated drugs
}
```

### 药物对象
```typescript
{
  id: string;              // ClinPGx drug ID
  name: string;            // Generic name
  tradeNames: string[];    // Brand names
  drugClasses: string[];   // Therapeutic classes
  indication: string;      // Primary indication
  relatedGenes: string[];  // Pharmacogenes
}
```

### 基因-药物对对象
```typescript
{
  gene: string;            // Gene symbol
  drug: string;            // Drug name
  sources: string[];       // CPIC, FDA, DPWG, etc.
  cpicLevel: string;       // A, B, C, D
  evidenceLevel: string;   // 1A, 1B, 2A, 2B, 3, 4
  hasGuideline: boolean;   // Has clinical guideline
}
```

### 等位基因对象
```typescript
{
  name: string;            // Allele name (e.g., CYP2D6*4)
  gene: string;            // Gene symbol
  function: string;        // Normal/decreased/no/increased/uncertain
  activityScore: number;   // 0.0 to 2.0+
  frequencies: {           // Population frequencies
    [population: string]: number;
  };
  definingVariants: string[];  // rsIDs or descriptions
}
```

## API 稳定性和版本控制

### 当前状态
- API版本：v1
- 稳定性：Beta - 端点稳定，参数可能会改变
- 监控：https://blog.clinpgx.org/ 的更新

### 从 PharmGKB 迁移
截至 2025 年 7 月，PharmGKB URL 重定向至 ClinPGx。更新参考：
- 旧：`https://api.pharmgkb.org/`
- 新：`https://api.clinpgx.org/`

### 未来的变化
- 关注 API v2 公告
- 重大变更将在 ClinPGx 博客上公布
- 考虑生产应用程序的版本固定

## 支持和联系

- **API 问题**：api@clinpgx.org
- **文档**：https://api.clinpgx.org/
- **一般问题**：https://www.clinpgx.org/page/faqs
- **博客**：https://blog.clinpgx.org/
- **太保指引**：https://cpicpgx.org/

## 相关资源

- **PharmCAT**：药物基因组变异调用和注释工具
- **PharmVar**：Pharmacogene 等位基因命名数据库
- **CPIC**：临床药物遗传学实施联盟
- **DPWG**：荷兰药物遗传学工作组
- **ClinGen**：临床基因组资源