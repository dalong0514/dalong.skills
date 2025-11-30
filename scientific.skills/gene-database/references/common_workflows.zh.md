<!-- 此文件由机器翻译自 common_workflows.md -->

# 通用基因数据库工作流程

本文档提供了使用 NCBI 基因数据库的常见工作流程和用例的示例。

## 目录

1.【疾病基因发现】(#disease-gene-discovery)
2. [基因注释管道](#gene-annotation-pipeline)
3. [跨物种基因比较](#cross-species-gene-comparison)
4. [通路分析](#pathway-analysis)
5. [变异分析](#variant-analysis)
6.【出版物挖矿】(#publication-mining)

---

## 疾病基因发现

### 用例

识别与特定疾病或表型相关的基因。

### 工作流程

1. **按疾病名称搜索**

```bash
# Find genes associated with Alzheimer's disease
python scripts/query_gene.py --search "Alzheimer disease[disease]" --organism human --max-results 50
```

2. **按染色体位置过滤**

<<<代码块_1>>>

3. **检索详细信息**

<<<代码块_2>>>

### 预期输出

- 与疾病相关的基因列表
- 基因符号、描述和染色体位置
- 相关出版物和临床注释

---

## 基因注释管道

### 用例

使用全面的元数据注释基因标识符列表。

### 工作流程

1. **准备基因列表**

创建一个带有基因符号的文件 `genes.txt`（每行一个）：
<<<代码块_3>>>

2. **批量查找**

<<<代码块_4>>>

3. **解析结果**

<<<代码块_5>>>

4. **丰富序列数据**

<<<代码块_6>>>

### 用例

- 为出版物创建基因注释表
- 分析前验证基因列表
- 建立基因参考数据库
- 基因组管道的质量控制

---

## 跨物种基因比较

### 用例

查找直向同源物或比较不同物种的相同基因。

### 工作流程

1. **在多种生物体中寻找基因**

```bash
# Find TP53 in human
python scripts/fetch_gene_data.py --symbol TP53 --taxon human

# Find TP53 in mouse
python scripts/fetch_gene_data.py --symbol TP53 --taxon mouse

# Find TP53 in zebrafish
python scripts/fetch_gene_data.py --symbol TP53 --taxon zebrafish
```

2. **比较不同物种的基因 ID**

```python
# Compare gene information across species
species = {
    'human': '9606',
    'mouse': '10090',
    'rat': '10116'
}

gene_symbol = 'TP53'

for organism, taxon_id in species.items():
    # Fetch gene data
    # ... (use fetch_gene_by_symbol)
    print(f"{organism}: {gene_data}")
```

3. **使用 ELink 查找直向同源物**

```bash
# Get HomoloGene links for a gene
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=homologene&id=7157&retmode=json"
```

### 应用程序

- 进化研究
- 模式生物研究
- 比较基因组学
- 跨物种实验设计

---

## 通路分析

### 用例

识别参与特定生物途径或过程的基因。

### 工作流程

1. **按基因本体 (GO) 术语搜索**

```bash
# Find genes involved in apoptosis
python scripts/query_gene.py --search "GO:0006915[biological process]" --organism human --max-results 100
```

2. **按途径名称搜索**

```bash
# Find genes in insulin signaling pathway
python scripts/query_gene.py --search "insulin signaling pathway[pathway]" --organism human
```

3. **获取通路相关基因**

```python
# Example: Get all genes in a specific pathway
import urllib.request
import json

# Search for pathway genes
query = "MAPK signaling pathway[pathway] AND human[organism]"
url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={query}&retmode=json&retmax=200"

with urllib.request.urlopen(url) as response:
    data = json.loads(response.read().decode())
    gene_ids = data['esearchresult']['idlist']

print(f"Found {len(gene_ids)} genes in MAPK signaling pathway")
```

4. **批量检索基因详细信息**

```bash
# Get details for all pathway genes
python scripts/batch_gene_lookup.py --ids 5594,5595,5603,5604 --output mapk_genes.json
```

### 应用程序

- 通路富集分析
- 基因集分析
- 系统生物学研究
- 药物靶点识别

---

## 变异分析

### 用例

查找具有临床相关变异或疾病相关突变的基因。

### 工作流程

1. **搜索具有临床变异的基因**

```bash
# Find genes with pathogenic variants
python scripts/query_gene.py --search "pathogenic[clinical significance]" --organism human --max-results 50
```

2. **链接到 ClinVar 数据库**

```bash
# Get ClinVar records for a gene
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=clinvar&id=672&retmode=json"
```

3. **寻找药物基因组基因**

```bash
# Find genes associated with drug response
python scripts/query_gene.py --search "pharmacogenomic[property]" --organism human
```

4. **获取变体汇总数据**

```python
# Example: Get genes with known variants
from scripts.query_gene import esearch, efetch

# Search for genes with variants
gene_ids = esearch("has variants[filter] AND human[organism]", retmax=100)

# Fetch detailed records
for gene_id in gene_ids[:10]:  # First 10
    data = efetch([gene_id], retmode='xml')
    # Parse XML for variant information
    print(f"Gene {gene_id} variant data...")
```

### 应用程序

- 临床遗传学
- 精准医疗
- 药物基因组学
- 遗传咨询

---

## 出版物挖掘

### 用例

查找最近出版物中提到的基因或将基因与文献联系起来。

### 工作流程

1. **搜索特定出版物中提到的基因**

```bash
# Find genes mentioned in papers about CRISPR
python scripts/query_gene.py --search "CRISPR[text word]" --organism human --max-results 100
```

2. **获取基因的 PubMed 文章**

```bash
# Get all publications for BRCA1
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=pubmed&id=672&retmode=json"
```

3. **按作者或期刊搜索**

```bash
# Find genes studied by specific research group
python scripts/query_gene.py --search "Smith J[author] AND 2024[pdat]" --organism human
```

4. **提取基因-发表关系**

```python
# Example: Build gene-publication network
from scripts.query_gene import esearch, esummary
import urllib.request
import json

# Get gene
gene_id = '672'

# Get publications for gene
url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=pubmed&id={gene_id}&retmode=json"

with urllib.request.urlopen(url) as response:
    data = json.loads(response.read().decode())

# Extract PMIDs
pmids = []
for linkset in data.get('linksets', []):
    for linksetdb in linkset.get('linksetdbs', []):
        pmids.extend(linksetdb.get('links', []))

print(f"Gene {gene_id} has {len(pmids)} publications")
```

### 应用程序

- 文献综述
- 拨款写作
- 知识库建设
- 基因组学研究趋势分析

---

## 高级模式

### 组合多个搜索

```python
# Example: Find genes at intersection of multiple criteria
def find_genes_multi_criteria(organism='human'):
    # Criteria 1: Disease association
    disease_genes = set(esearch("diabetes[disease] AND human[organism]"))

    # Criteria 2: Chromosome location
    chr_genes = set(esearch("11[chromosome] AND human[organism]"))

    # Criteria 3: Gene type
    coding_genes = set(esearch("protein coding[gene type] AND human[organism]"))

    # Intersection
    candidates = disease_genes & chr_genes & coding_genes

    return list(candidates)
```

### 限速批处理

```python
import time

def process_genes_with_rate_limit(gene_ids, batch_size=200, delay=0.1):
    results = []

    for i in range(0, len(gene_ids), batch_size):
        batch = gene_ids[i:i + batch_size]

        # Process batch
        batch_results = esummary(batch)
        results.append(batch_results)

        # Rate limit
        time.sleep(delay)

    return results
```

### 错误处理和重试

```python
import time

def robust_gene_fetch(gene_id, max_retries=3):
    for attempt in range(max_retries):
        try:
            data = fetch_gene_by_id(gene_id)
            return data
        except Exception as e:
            if attempt < max_retries - 1:
                wait = 2 ** attempt  # Exponential backoff
                time.sleep(wait)
            else:
                print(f"Failed to fetch gene {gene_id}: {e}")
                return None
```

---

## 提示和最佳实践

1. **开始具体，然后扩大**：从精确查询开始，并根据需要进行扩展
2. **使用生物体过滤器**：始终指定生物体进行基因符号搜索
3. **验证结果**：检查基因 ID 和符号的准确性
4. **缓存常用数据**：将常用查询存储在本地
5. **监控速率限制**：使用 API 密钥并实施延迟
6. **组合 API**：使用电子实用程序进行搜索，使用数据集 API 获取详细数据
7. **处理歧义**：基因符号可能指不同物种的不同基因
8. **检查数据流通**：基因注释定期更新
9. **使用批量操作**：尽可能一起处理多个基因
10. **记录您的查询**：保留搜索词和参数的记录