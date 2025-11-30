<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：clinpgx-数据库
描述：“访问 ClinPGx 药物基因组学数据（PharmGKB 的后继者）。查询基因-药物相互作用、CPIC 指南、等位基因功能，以实现精准医疗和基因型指导的剂量决策。”
---

# ClinPGx 数据库

## 概述

ClinPGx（临床药物基因组学数据库）是临床药物基因组学信息的综合资源，是 PharmGKB 的继承者。它整合了 PharmGKB、CPIC 和 PharmCAT 的数据，提供有关遗传变异如何影响药物反应的精选信息。获取基因药物对、临床指南、等位基因功能和药物标签，以实现精准医学应用。

## 何时使用此技能

该技能应该在以下情况下使用：

- **基因-药物相互作用**：查询基因变异如何影响药物代谢、功效或毒性
- **CPIC 指南**：获取基于证据的药物遗传学临床实践指南
- **等位基因信息**：检索等位基因功能、频率和表型数据
- **药物标签**：探索 FDA 和其他监管药物基因组药物标签
- **药物基因组注释**：访问有关基因-药物-疾病关系的精选文献
- **临床决策支持**：使用 PharmDOG 工具进行表型转化和定制基因型解释
- **精准医学**：在临床实践中实施药物基因组学测试
- **药物代谢**：了解 CYP450 和其他药物基因功能
- **个性化剂量**：寻找基因型指导的剂量建议
- **药物不良反应**：识别药物毒性的遗传风险因素

## 安装和设置

### Python API 访问

ClinPGx REST API 提供对所有数据库资源的编程访问。基本设置：

```bash
uv pip install requests
```

### API 端点

<<<代码块_1>>>

**速率限制**：
- 每秒最多 2 个请求
- 过多的请求将导致 HTTP 429（请求过多）响应

**身份验证**：基本访问不需要

**数据许可证**：知识共享署名-相同方式共享 4.0 国际许可证

对于大量 API 使用，请通过 api@clinpgx.org 通知 ClinPGx 团队

## 核心能力

### 1. 基因查询

**检索基因信息**，包括功能、临床注释和药物基因组学意义：

<<<代码块_2>>>

**关键药效基因**：
- **CYP450 酶**：CYP2D6、CYP2C19、CYP2C9、CYP3A4、CYP3A5
- **运输机**：SLCO1B1、ABCB1、ABCG2
- **其他代谢物**：TPMT、DPYD、NUDT15、UGT1A1
- **受体**：OPRM1、HTR2A、ADRB1
- **HLA 基因**：HLA-B、HLA-A

### 2. 药物和化学品查询

**检索药物信息**，包括药物基因组注释和机制：

<<<代码块_3>>>

**具有药物基因组学意义的药物类别**：
- 抗凝剂（华法林、氯吡格雷）
- 抗抑郁药（SSRI、TCA）
- 免疫抑制剂（他克莫司、硫唑嘌呤）
- 肿瘤药物（5-氟尿嘧啶、伊立替康、他莫昔芬）
- 心血管药物（他汀类药物、β受体阻滞剂）
- 止痛药（可待因、曲马多）
- 抗病毒药（阿巴卡韦）

### 3. 基因-药物对查询

**通过临床注释访问精心策划的基因-药物关系**：

<<<代码块_4>>>

**临床注释来源**：
- CPIC（临床药物遗传学实施联盟）
- DPWG（荷兰药物遗传学工作组）
- FDA（食品和药物管理局）标签
- 同行评审的文献摘要注释

### 4.太保指引

**获取基于证据的临床实践指南**：

<<<代码块_5>>>

**CPIC 指南组成部分**：
- 涵盖基因-药物对
- 按表型的临床建议
- 证据级别和强度评级
- 支持文献
- 可下载的 PDF 和补充材料
- 实施注意事项

**指南示例**：
- CYP2D6-可待因（避免用于超快速代谢者）
- CYP2C19-氯吡格雷（代谢不良者的替代疗法）
- TPMT-硫唑嘌呤（中间/弱代谢者减少剂量）
- DPYD-氟嘧啶（根据活性调整剂量）
- HLA-B*57:01-阿巴卡韦（如果呈阳性则避免）

### 5. 等位基因和变异信息

**查询等位基因功能和频率数据**：

<<<代码块_6>>>

**等位基因信息包括**：
- 功能状态（正常、下降、无功能、增加、不确定）
- 各族裔群体的人口频率
- 定义变异（SNP、插入缺失、CNV）
- 表型分配
- 对 PharmVar 和其他命名系统的引用

**表型类别**：
- **超快速代谢器** (UM)：增加酶活性
- **正常代谢者** (NM)：正常酶活性
- **中间代谢物** (IM)：酶活性降低
- **代谢能力差** (PM)：酶活性很少甚至没有

### 6. 变体注释

**访问特定遗传变异的临床注释**：

```python
# Get variant information
response = requests.get("https://api.clinpgx.org/v1/variant/rs4244285")
variant_data = response.json()

# Search variants by position (if supported)
response = requests.get("https://api.clinpgx.org/v1/variant",
                       params={"chromosome": "10", "position": "94781859"})
variants = response.json()
```

**变体数据包括**：
- rsID和基因组坐标
- 基因和功能后果
- 等位基因关联
- 临床意义
- 人口频率
- 文献参考

### 7. 临床注释

**检索精选文献注释**（以前称为 PharmGKB 临床注释）：

```python
# Get clinical annotations
response = requests.get("https://api.clinpgx.org/v1/clinicalAnnotation",
                       params={"gene": "CYP2D6"})
annotations = response.json()

# Filter by evidence level
response = requests.get("https://api.clinpgx.org/v1/clinicalAnnotation",
                       params={"evidenceLevel": "1A"})
high_evidence = response.json()
```

**证据级别**（从最高到最低）：
- **1A级**：高质量证据，CPIC/FDA/DPWG指南
- **1B 级**：高质量证据，尚未成为指南
- **2A 级**：来自精心设计的研究的适度证据
- **2B 级**：中等证据，但有一些局限性
- **3级**：证据有限或相互矛盾
- **4级**：病例报告或证据薄弱

### 8. 药品标签

**从药物标签获取药物基因组信息**：

```python
# Get drug labels with PGx information
response = requests.get("https://api.clinpgx.org/v1/drugLabel",
                       params={"drug": "warfarin"})
labels = response.json()

# Filter by regulatory source
response = requests.get("https://api.clinpgx.org/v1/drugLabel",
                       params={"source": "FDA"})
fda_labels = response.json()
```

**标签信息包括**：
- 测试建议
- 按基因型的剂量指导
- 警告和预防措施
- 生物标志物信息
- 监管来源（FDA、EMA、PMDA 等）

### 9. 途径

**探索药代动力学和药效学途径**：

```python
# Get pathway information
response = requests.get("https://api.clinpgx.org/v1/pathway/PA146123006")  # Warfarin pathway
pathway_data = response.json()

# Search pathways by drug
response = requests.get("https://api.clinpgx.org/v1/pathway",
                       params={"drug": "warfarin"})
pathways = response.json()
```

**路径图**显示：
- 药物代谢步骤
- 涉及酶和转运蛋白
- 影响每一步的基因变异
- 对功效/毒性的下游影响
- 与其他途径的相互作用

## 查询工作流程

### 工作流程 1：药物处方的临床决策支持

1. **确定患者基因型**的相关药物基因：
   ```python
   # Example: Patient is CYP2C19 *1/*2 (intermediate metabolizer)
   response = requests.get("https://api.clinpgx.org/v1/allele/CYP2C19*2")
   allele_function = response.json()
   ```

2. **查询基因-药物对**以获取感兴趣的药物：
   ```python
   response = requests.get("https://api.clinpgx.org/v1/geneDrugPair",
                          params={"gene": "CYP2C19", "drug": "clopidogrel"})
   pair_info = response.json()
   ```

3. **检索 CPIC 指南**以获取剂量建议：
   ```python
   response = requests.get("https://api.clinpgx.org/v1/guideline",
                          params={"gene": "CYP2C19", "drug": "clopidogrel"})
   guideline = response.json()
   # Recommendation: Alternative antiplatelet therapy for IM/PM
   ```

4. **检查药品标签**以获取监管指导：
   ```python
   response = requests.get("https://api.clinpgx.org/v1/drugLabel",
                          params={"drug": "clopidogrel"})
   label = response.json()
   ```

### 工作流程 2：基因组分析

1. **获取临床面板中的药物基因列表**：
   ```python
   pgx_panel = ["CYP2C19", "CYP2D6", "CYP2C9", "TPMT", "DPYD", "SLCO1B1"]
   ```

2. **对于每个基因，检索所有药物相互作用**：
   ```python
   all_interactions = {}
   for gene in pgx_panel:
       response = requests.get("https://api.clinpgx.org/v1/geneDrugPair",
                              params={"gene": gene})
       all_interactions[gene] = response.json()
   ```

3. **筛选 CPIC 指南级证据**：
   ```python
   for gene, pairs in all_interactions.items():
       for pair in pairs:
           if pair.get('cpicLevel'):  # Has CPIC guideline
               print(f"{gene} - {pair['drug']}: {pair['cpicLevel']}")
   ```

4. **生成具有可操作的药物基因组学发现的患者报告**。

### 工作流程 3：药物安全评估

1. **查询药物的 PGx 关联**：
   ```python
   response = requests.get("https://api.clinpgx.org/v1/chemical",
                          params={"name": "abacavir"})
   drug_id = response.json()[0]['id']
   ```

2. **获取临床注释**：
   ```python
   response = requests.get("https://api.clinpgx.org/v1/clinicalAnnotation",
                          params={"drug": drug_id})
   annotations = response.json()
   ```

3. **检查 HLA 关联**和毒性风险：
   ```python
   for annotation in annotations:
       if 'HLA' in annotation.get('genes', []):
           print(f"Toxicity risk: {annotation['phenotype']}")
           print(f"Evidence level: {annotation['evidenceLevel']}")
   ```

4. **从指南和标签中检索筛查建议**。

### 工作流程 4：研究分析 - 群体药物基因组学

1. **获取等位基因频率**以进行群体比较：
   ```python
   response = requests.get("https://api.clinpgx.org/v1/allele",
                          params={"gene": "CYP2D6"})
   alleles = response.json()
   ```

2. **提取特定人群的频率**：
   ```python
   populations = ['European', 'African', 'East Asian', 'Latino']
   frequency_data = {}
   for allele in alleles:
       allele_name = allele['name']
       frequency_data[allele_name] = {
           pop: allele.get(f'{pop}_frequency', 'N/A')
           for pop in populations
       }
   ```

3. **按人群计算表型分布**：
   ```python
   # Combine allele frequencies with function to predict phenotypes
   phenotype_dist = calculate_phenotype_frequencies(frequency_data)
   ```

4. **分析不同人群药物剂量的影响**。

### 工作流程 5：文献证据审查

1. **搜索基因-药物对**：
   ```python
   response = requests.get("https://api.clinpgx.org/v1/geneDrugPair",
                          params={"gene": "TPMT", "drug": "azathioprine"})
   pair = response.json()
   ```

2. **检索所有临床注释**：
   ```python
   response = requests.get("https://api.clinpgx.org/v1/clinicalAnnotation",
                          params={"gene": "TPMT", "drug": "azathioprine"})
   annotations = response.json()
   ```

3. **按证据级别和发布日期过滤**：
   ```python
   high_quality = [a for a in annotations
                   if a['evidenceLevel'] in ['1A', '1B', '2A']]
   ```

4. **提取 PMID** 并检索完整参考文献：
   ```python
   pmids = [a['pmid'] for a in high_quality if 'pmid' in a]
   # Use PubMed skill to retrieve full citations
   ```

## 速率限制和最佳实践

### 速率限制合规性

```python
import time

def rate_limited_request(url, params=None, delay=0.5):
    """Make API request with rate limiting (2 req/sec max)"""
    response = requests.get(url, params=params)
    time.sleep(delay)  # Wait 0.5 seconds between requests
    return response

# Use in loops
genes = ["CYP2D6", "CYP2C19", "CYP2C9"]
for gene in genes:
    response = rate_limited_request(
        "https://api.clinpgx.org/v1/gene/" + gene
    )
    data = response.json()
```

### 错误处理

```python
def safe_api_call(url, params=None, max_retries=3):
    """API call with error handling and retries"""
    for attempt in range(max_retries):
        try:
            response = requests.get(url, params=params, timeout=10)

            if response.status_code == 200:
                return response.json()
            elif response.status_code == 429:
                # Rate limit exceeded
                wait_time = 2 ** attempt  # Exponential backoff
                print(f"Rate limit hit. Waiting {wait_time}s...")
                time.sleep(wait_time)
            else:
                response.raise_for_status()

        except requests.exceptions.RequestException as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            if attempt == max_retries - 1:
                raise
            time.sleep(1)
```

### 缓存结果

```python
import json
from pathlib import Path

def cached_query(cache_file, api_func, *args, **kwargs):
    """Cache API results to avoid repeated queries"""
    cache_path = Path(cache_file)

    if cache_path.exists():
        with open(cache_path) as f:
            return json.load(f)

    result = api_func(*args, **kwargs)

    with open(cache_path, 'w') as f:
        json.dump(result, f, indent=2)

    return result

# Usage
gene_data = cached_query(
    'cyp2d6_cache.json',
    rate_limited_request,
    "https://api.clinpgx.org/v1/gene/CYP2D6"
)
```

## PharmDOG 工具

PharmDOG（以前称为 DDRx）是 ClinPGx 的临床决策支持工具，用于解释药物基因组测试结果：

**主要特点**：
- **表型转化计算器**：调整影响 CYP2D6 的药物相互作用的表型预测
- **自定义基因型**：输入患者基因型以获得表型预测
- **二维码共享**：生成可共享的患者报告
- **灵活的指导来源**：选择要应用的指南（CPIC、DPWG、FDA）
- **多药物分析**：同时评估多种药物

**访问**：可通过 https://www.clinpgx.org/pharmacogenomic-decision-support 获取

**用例**：
- PGx 小组结果的临床解释
- 已知基因型患者的药物审查
- 患者教育材料
- 即时护理决策支持

## 资源

### 脚本/query_clinpgx.py

Python 脚本，具有用于常见 ClinPGx 查询的即用型函数：

- `get_gene_info(gene_symbol)` - 检索基因详细信息
- `get_drug_info(drug_name)` - 获取药品信息
- `get_gene_drug_pairs(gene, drug)` - 查询基因-药物相互作用
- `get_cpic_guidelines(gene, drug)` - 检索 CPIC 指南
- `get_alleles(gene)` - 获取基因的所有等位基因
- `get_clinical_annotations(gene, drug, evidence_level)` - 查询文献注释
- `get_drug_labels(drug)` - 检索药物基因组药物标签
- `search_variants(rsid)` - 按变体 rsID 搜索
- `export_to_dataframe(data)` - 将结果转换为 pandas DataFrame

请参阅此脚本以获取具有适当速率限制和错误处理的实现示例。

### 参考文献/api_reference.md

全面的 API 文档包括：

- 带参数的完整端点列表
- 请求/响应格式规范
- 每个端点的查询示例
- 过滤运算符和搜索模式
- 数据模式定义
- 速率限制细节
- 身份验证要求（如果有）
- 排除常见错误

当需要详细的 API 信息或构建复杂查询时，请参阅此文档。

## 重要提示

### 数据源和集成

ClinPGx 整合了多个权威来源：
- **PharmGKB**：精心策划的药物基因组学知识库（现在是 ClinPGx 的一部分）
- **CPIC**：循证临床实施指南
- **PharmCAT**：等位基因调用和表型解释工具
- **DPWG**：荷兰药物遗传学指南
- **FDA/EMA 标签**：监管药物基因组信息

截至 2025 年 7 月，所有 PharmGKB URL 都重定向到相应的 ClinPGx 页面。

### 临床实施注意事项

- **证据水平**：临床应用前务必检查证据强度
- **人群差异**：等位基因频率在不同人群中差异显着
- **表型转化**：考虑影响酶活性的药物间相互作用
- **多基因效应**：一些药物受多种药物基因影响
- **非遗传因素**：年龄、器官功能、药物相互作用也会影响反应
- **测试限制**：并非所有检测都检测到所有临床相关等位基因

### 数据更新

- ClinPGx 不断更新新的证据和指南
- 检查临床注释的出版日期
- 监控 ClinPGx 博客 (https://blog.clinpgx.org/) 以获取公告
- 随着新证据的出现，太保指南更新
- PharmVar 提供等位基因定义的命名法更新

### API 稳定性

- API端点相对稳定，但在开发过程中可能会发生变化
- 参数和响应格式可能会被修改
- 监控 API 变更日志和 ClinPGx 博客以获取更新
- 考虑生产应用程序的版本固定
- 在生产部署之前测试开发中的 API 更改

## 常见用例

### 预防性药物基因组测试

查询所有临床上可行的基因药物对以指导面板选择：

```python
# Get all CPIC guideline pairs
response = requests.get("https://api.clinpgx.org/v1/geneDrugPair",
                       params={"cpicLevel": "A"})  # Level A recommendations
actionable_pairs = response.json()
```

### 药物治疗管理

根据已知基因型审查患者用药：

```python
patient_genes = {"CYP2C19": "*1/*2", "CYP2D6": "*1/*1", "SLCO1B1": "*1/*5"}
medications = ["clopidogrel", "simvastatin", "escitalopram"]

for med in medications:
    for gene in patient_genes:
        response = requests.get("https://api.clinpgx.org/v1/geneDrugPair",
                               params={"gene": gene, "drug": med})
        # Check for interactions and dosing guidance
```

### 临床试验资格

筛查药物基因组禁忌症：

```python
# Check for HLA-B*57:01 before abacavir trial
response = requests.get("https://api.clinpgx.org/v1/geneDrugPair",
                       params={"gene": "HLA-B", "drug": "abacavir"})
pair_info = response.json()
# CPIC: Do not use if HLA-B*57:01 positive
```

## 其他资源

- **ClinPGx 网站**：https://www.clinpgx.org/
- **ClinPGx 博客**：https://blog.clinpgx.org/
- **API 文档**：https://api.clinpgx.org/
- **太保网站**：https://cpicpgx.org/
- **PharmCAT**：https://pharmcat.clinpgx.org/
- **ClinGen**：https://clinicalgenome.org/
- **联系方式**：api@clinpgx.org（用于大量 API 使用）