<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：pubchem-数据库
描述：“通过 PUG-REST API/PubChemPy（110M+ 化合物）查询 PubChem。按名称/CID/SMILES 搜索，检索属性、相似性/子结构搜索、生物活性，用于化学信息学。”
---

# PubChem 数据库

## 概述

PubChem 是世界上最大的免费化学数据库，包含 1.1 亿多种化合物和 2.7 亿多种生物活性。按名称、CID 或 SMILES 查询化学结构，检索分子属性，执行相似性和子结构搜索，使用 PUG-REST API 和 PubChemPy 访问生物活性数据。

## 何时使用此技能

该技能应该在以下情况下使用：
- 按名称、结构 (SMILES/InChI) 或分子式搜索化合物
- 检索分子特性（MW、LogP、TPSA、氢键描述符）
- 执行相似性搜索以查找结构相关的化合物
- 对特定化学基序进行子结构搜索
- 从筛选分析中获取生物活性数据
- 化学标识符格式之间的转换（CID、SMILES、InChI）
- 批量处理多种化合物以进行药物相似性筛选或性质分析

## 核心能力

### 1.化学结构搜索

使用多种标识符类型搜索化合物：

**按化学名称**：
```python
import pubchempy as pcp
compounds = pcp.get_compounds('aspirin', 'name')
compound = compounds[0]
```

**按 CID（化合物 ID）**：
<<<代码块_1>>>

**通过微笑**：
<<<代码块_2>>>

**作者：InChI**：
<<<代码块_3>>>

**按分子式**：
<<<代码块_4>>>

### 2. 财产检索

使用高级或低级方法检索化合物的分子特性：

**使用 PubChemPy（推荐）**：
<<<代码块_5>>>

**获取特定属性**：
<<<代码块_6>>>

**批量属性检索**：
```python
import pandas as pd

compound_names = ['aspirin', 'ibuprofen', 'paracetamol']
all_properties = []

for name in compound_names:
    props = pcp.get_properties(
        ['MolecularFormula', 'MolecularWeight', 'XLogP'],
        name,
        'name'
    )
    all_properties.extend(props)

df = pd.DataFrame(all_properties)
```

**可用属性**：MolecularFormula、MolecularWeight、CanonicalSMILES、IsomericSMILES、InChI、InChIKey、IUPACName、XLogP、TPSA、HBondDonorCount、HBondAcceptorCount、RotatableBondCount、Complexity、Charge 等（完整列表请参阅 `references/api_reference.md`）。

### 3.相似性搜索

使用 Tanimoto 相似性查找结构相似的化合物：

```python
import pubchempy as pcp

# Start with a query compound
query_compound = pcp.get_compounds('gefitinib', 'name')[0]
query_smiles = query_compound.canonical_smiles

# Perform similarity search
similar_compounds = pcp.get_compounds(
    query_smiles,
    'smiles',
    searchtype='similarity',
    Threshold=85,  # Similarity threshold (0-100)
    MaxRecords=50
)

# Process results
for compound in similar_compounds[:10]:
    print(f"CID {compound.cid}: {compound.iupac_name}")
    print(f"  MW: {compound.molecular_weight}")
```

**注意**：对于大型查询，相似性搜索是异步的，可能需要 15-30 秒才能完成。 PubChemPy 自动处理异步模式。

### 4. 子结构搜索

查找包含特定结构基序的化合物：

```python
import pubchempy as pcp

# Search for compounds containing pyridine ring
pyridine_smiles = 'c1ccncc1'

matches = pcp.get_compounds(
    pyridine_smiles,
    'smiles',
    searchtype='substructure',
    MaxRecords=100
)

print(f"Found {len(matches)} compounds containing pyridine")
```

**常见子结构**：
- 苯环：`c1ccccc1`
- 吡啶：`c1ccncc1`
- 苯酚：`c1ccc(O)cc1`
- 羧酸：`C(=O)O`

### 5.格式转换

不同化学结构格式之间的转换：

```python
import pubchempy as pcp

compound = pcp.get_compounds('aspirin', 'name')[0]

# Convert to different formats
smiles = compound.canonical_smiles
inchi = compound.inchi
inchikey = compound.inchikey
cid = compound.cid

# Download structure files
pcp.download('SDF', 'aspirin', 'name', 'aspirin.sdf', overwrite=True)
pcp.download('JSON', '2244', 'cid', 'aspirin.json', overwrite=True)
```

### 6. 结构可视化

生成 2D 结构图像：

```python
import pubchempy as pcp

# Download compound structure as PNG
pcp.download('PNG', 'caffeine', 'name', 'caffeine.png', overwrite=True)

# Using direct URL (via requests)
import requests

cid = 2244  # Aspirin
url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG?image_size=large"
response = requests.get(url)

with open('structure.png', 'wb') as f:
    f.write(response.content)
```

### 7. 同义词检索

获取化合物的所有已知名称和同义词：

```python
import pubchempy as pcp

synonyms_data = pcp.get_synonyms('aspirin', 'name')

if synonyms_data:
    cid = synonyms_data[0]['CID']
    synonyms = synonyms_data[0]['Synonym']

    print(f"CID {cid} has {len(synonyms)} synonyms:")
    for syn in synonyms[:10]:  # First 10
        print(f"  - {syn}")
```

### 8. 生物活性数据访问

从测定中检索生物活性数据：

```python
import requests
import json

# Get bioassay summary for a compound
cid = 2244  # Aspirin
url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"

response = requests.get(url)
if response.status_code == 200:
    data = response.json()
    # Process bioassay information
    table = data.get('Table', {})
    rows = table.get('Row', [])
    print(f"Found {len(rows)} bioassay records")
```

**对于更复杂的生物活性查询**，请使用 `scripts/bioactivity_query.py` 帮助程序脚本，它提供：
- 具有活动结果过滤的生物测定摘要
- 检测目标识别
- 按生物靶标搜索化合物
- 用于特定测定的活性化合物列表

### 9. 综合复合注释

通过 PUG-View 访问详细的化合物信息：

```python
import requests

cid = 2244
url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"

response = requests.get(url)
if response.status_code == 200:
    annotations = response.json()
    # Contains extensive data including:
    # - Chemical and Physical Properties
    # - Drug and Medication Information
    # - Pharmacology and Biochemistry
    # - Safety and Hazards
    # - Toxicity
    # - Literature references
    # - Patents
```

**获取特定部分**：
```python
# Get only drug information
url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON?heading=Drug and Medication Information"
```

## 安装要求

安装 PubChemPy 以进行基于 Python 的访问：

```bash
uv pip install pubchempy
```

对于直接 API 访问和生物活性查询：

```bash
uv pip install requests
```

数据分析可选：

```bash
uv pip install pandas
```

## 帮助脚本

此技能包括用于常见 PubChem 任务的 Python 脚本：

### 脚本/compound_search.py

提供用于搜索和检索化合物信息的实用函数：

**主要功能**：
- `search_by_name(name, max_results=10)`：按名称搜索化合物
- `search_by_smiles(smiles)`：按 SMILES 字符串搜索
- `get_compound_by_cid(cid)`：通过 CID 检索化合物
- `get_compound_properties(identifier, namespace, properties)`：获取特定属性
- `similarity_search(smiles, threshold, max_records)`：执行相似性搜索
- `substructure_search(smiles, max_records)`：执行子结构搜索
- `get_synonyms(identifier, namespace)`：获取所有同义词
- `batch_search(identifiers, namespace, properties)`：批量搜索多个化合物
- `download_structure(identifier, namespace, format, filename)`：下载结构
- `print_compound_info(compound)`：打印格式化的化合物信息

**用法**：
```python
from scripts.compound_search import search_by_name, get_compound_properties

# Search for a compound
compounds = search_by_name('ibuprofen')

# Get specific properties
props = get_compound_properties('aspirin', 'name', ['MolecularWeight', 'XLogP'])
```

### 脚本/bioactivity_query.py

提供检索生物活性数据的功能：

**主要功能**：
- `get_bioassay_summary(cid)`：获取化合物的生物测定摘要
- `get_compound_bioactivities(cid, activity_outcome)`：获取过滤的生物活性
- `get_assay_description(aid)`：获取详细的检测信息
- `get_assay_targets(aid)`：获取用于分析的生物目标
- `search_assays_by_target(target_name, max_results)`：按目标查找检测
- `get_active_compounds_in_assay(aid, max_results)`：获取活性化合物
- `get_compound_annotations(cid, section)`：获取 PUG-View 注释
- `summarize_bioactivities(cid)`：生成生物活性汇总统计数据
- `find_compounds_by_bioactivity(target, threshold, max_compounds)`：按目标查找化合物

**用法**：
```python
from scripts.bioactivity_query import get_bioassay_summary, summarize_bioactivities

# Get bioactivity summary
summary = summarize_bioactivities(2244)  # Aspirin
print(f"Total assays: {summary['total_assays']}")
print(f"Active: {summary['active']}, Inactive: {summary['inactive']}")
```

## API 速率限制和最佳实践

**速率限制**：
- 每秒最多 5 个请求
- 每分钟最多 400 个请求
- 每分钟最长运行时间 300 秒

**最佳实践**：
1. **使用CID进行重复查询**：CID比名称或结构更高效
2. **本地缓存结果**：存储经常访问的数据
3. **批量请求**：尽可能组合多个查询
4. **实施延迟**：在请求之间添加 0.2-0.3 秒的延迟
5. **优雅地处理错误**：检查HTTP错误和丢失的数据
6. **使用 PubChemPy**：更高级别的抽象处理许多边缘情况
7. **利用异步模式**：用于大型相似性/子结构搜索
8. **指定MaxRecords**：限制结果以避免超时

**错误处理**：
```python
from pubchempy import BadRequestError, NotFoundError, TimeoutError

try:
    compound = pcp.get_compounds('query', 'name')[0]
except NotFoundError:
    print("Compound not found")
except BadRequestError:
    print("Invalid request format")
except TimeoutError:
    print("Request timed out - try reducing scope")
except IndexError:
    print("No results returned")
```

## 常见工作流程

### 工作流程 1：化学标识符转换管道

不同化学标识符之间的转换：

```python
import pubchempy as pcp

# Start with any identifier type
compound = pcp.get_compounds('caffeine', 'name')[0]

# Extract all identifier formats
identifiers = {
    'CID': compound.cid,
    'Name': compound.iupac_name,
    'SMILES': compound.canonical_smiles,
    'InChI': compound.inchi,
    'InChIKey': compound.inchikey,
    'Formula': compound.molecular_formula
}
```

### 工作流程 2：类药特性筛选

使用 Lipinski 五法则筛选化合物：

```python
import pubchempy as pcp

def check_drug_likeness(compound_name):
    compound = pcp.get_compounds(compound_name, 'name')[0]

    # Lipinski's Rule of Five
    rules = {
        'MW <= 500': compound.molecular_weight <= 500,
        'LogP <= 5': compound.xlogp <= 5 if compound.xlogp else None,
        'HBD <= 5': compound.h_bond_donor_count <= 5,
        'HBA <= 10': compound.h_bond_acceptor_count <= 10
    }

    violations = sum(1 for v in rules.values() if v is False)
    return rules, violations

rules, violations = check_drug_likeness('aspirin')
print(f"Lipinski violations: {violations}")
```

### 工作流程 3：寻找相似的候选药物

识别与已知药物结构相似的化合物：

```python
import pubchempy as pcp

# Start with known drug
reference_drug = pcp.get_compounds('imatinib', 'name')[0]
reference_smiles = reference_drug.canonical_smiles

# Find similar compounds
similar = pcp.get_compounds(
    reference_smiles,
    'smiles',
    searchtype='similarity',
    Threshold=85,
    MaxRecords=20
)

# Filter by drug-like properties
candidates = []
for comp in similar:
    if comp.molecular_weight and 200 <= comp.molecular_weight <= 600:
        if comp.xlogp and -1 <= comp.xlogp <= 5:
            candidates.append(comp)

print(f"Found {len(candidates)} drug-like candidates")
```

### 工作流程 4：批量化合物特性比较

比较多种化合物的特性：

```python
import pubchempy as pcp
import pandas as pd

compound_list = ['aspirin', 'ibuprofen', 'naproxen', 'celecoxib']

properties_list = []
for name in compound_list:
    try:
        compound = pcp.get_compounds(name, 'name')[0]
        properties_list.append({
            'Name': name,
            'CID': compound.cid,
            'Formula': compound.molecular_formula,
            'MW': compound.molecular_weight,
            'LogP': compound.xlogp,
            'TPSA': compound.tpsa,
            'HBD': compound.h_bond_donor_count,
            'HBA': compound.h_bond_acceptor_count
        })
    except Exception as e:
        print(f"Error processing {name}: {e}")

df = pd.DataFrame(properties_list)
print(df.to_string(index=False))
```

### 工作流程 5：基于子结构的虚拟筛选

筛选含有特定药效基团的化合物：

```python
import pubchempy as pcp

# Define pharmacophore (e.g., sulfonamide group)
pharmacophore_smiles = 'S(=O)(=O)N'

# Search for compounds containing this substructure
hits = pcp.get_compounds(
    pharmacophore_smiles,
    'smiles',
    searchtype='substructure',
    MaxRecords=100
)

# Further filter by properties
filtered_hits = [
    comp for comp in hits
    if comp.molecular_weight and comp.molecular_weight < 500
]

print(f"Found {len(filtered_hits)} compounds with desired substructure")
```

## 参考文档

有关详细的 API 文档，包括完整的属性列表、URL 模式、高级查询选项和更多示例，请参阅 `references/api_reference.md`。该综合参考包括：

- 完整的 PUG-REST API 端点文档
- 可用分子特性的完整列表
- 异步请求处理模式
- PubChemPy API 参考
- 用于注释的 PUG-View API
- 常见的工作流程和用例
- PubChem 官方文档的链接

## 故障排除

**未找到化合物**：
- 尝试其他名称或同义词
- 如果已知，请使用 CID
- 检查拼写和化学名称格式

**超时错误**：
- 减少MaxRecords参数
- 在请求之间添加延迟
- 使用 CID 代替名称以加快查询速度

**空属性值**：
- 并非所有化合物都具有所有属性
- 访问前检查属性是否存在：`if compound.xlogp:`
- 某些属性仅适用于某些化合物类型

**超出速率限制**：
- 在请求之间实现延迟（0.2-0.3 秒）
- 尽可能使用批处理操作
- 考虑在本地缓存结果

**相似性/子结构搜索挂起**：
- 这些是异步操作，可能需要 15-30 秒
- PubChemPy 自动处理轮询
- 如果超时则减少 MaxRecords

## 其他资源

- PubChem 主页：https://pubchem.ncbi.nlm.nih.gov/
- PUG-REST 文档：https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
- PUG-REST 教程：https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest-tutorial
- PubChemPy 文档：https://pubchempy.readthedocs.io/
- PubChemPy GitHub：https://github.com/mcs07/PubChemPy