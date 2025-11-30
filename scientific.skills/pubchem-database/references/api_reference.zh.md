<!-- 此文件由机器翻译自 api_reference.md -->

# PubChem API 参考

## 概述

PubChem 是世界上最大的免费化学数据库，由国家生物技术信息中心 (NCBI) 维护。它包含来自 770 多个数据源的超过 1.1 亿个独特的化学结构和超过 2.7 亿个生物活性。

## 数据库结构

PubChem 包含三个主要子数据库：

1. **化合物数据库**：具有计算属性的经过验证的独特化学结构
2. **物质数据库**：存放来自数据源的化学物质记录
3. **BioAssay数据库**：化合物的生物活性测试结果

## PubChem PUG-REST API

### 基本 URL 结构

```
https://pubchem.ncbi.nlm.nih.gov/rest/pug/<input>/<operation>/<output>
```

组件：
- `<input>`：化合物/cid、物质/sid、测定/辅助或搜索规范
- `<operation>`：可选操作，如属性、同义词、分类等。
- `<output>`：JSON、XML、CSV、PNG、SDF 等格式。

### 常见请求模式

#### 1. 通过标识符检索

通过CID（化合物ID）获取化合物：
<<<代码块_1>>>

按名称获取化合物：
<<<代码块_2>>>

通过 SMILES 获得复合：
<<<代码块_3>>>

通过 InChI 得到化合物：
<<<代码块_4>>>

#### 2. 可用属性

可检索的常见分子特性：
- `MolecularFormula`
- `MolecularWeight`
- `CanonicalSMILES`
- `IsomericSMILES`
- `InChI`
- `InChIKey`
- `IUPACName`
- `XLogP`
- `ExactMass`
- `MonoisotopicMass`
- `TPSA`（拓扑极表面积）
- `Complexity`
- `Charge`
- `HBondDonorCount`
- `HBondAcceptorCount`
- `RotatableBondCount`
- `HeavyAtomCount`
- `IsotopeAtomCount`
- `AtomStereoCount`
- `BondStereoCount`
- `CovalentUnitCount`
- `Volume3D`
- `XStericQuadrupole3D`
- `YStericQuadrupole3D`
- `ZStericQuadrupole3D`
- `FeatureCount3D`

要检索多个属性，请用逗号分隔它们：
<<<代码块_5>>>

#### 3.结构搜索操作

**相似性搜索**：
<<<代码块_6>>>

**子结构搜索**：
```
POST /rest/pug/compound/substructure/smiles/{smiles}/cids/JSON
```

**上层建筑搜索**：
```
POST /rest/pug/compound/superstructure/smiles/{smiles}/cids/JSON
```

#### 4. 图像生成

获取二维结构图像：
```
GET /rest/pug/compound/cid/{cid}/PNG
Optional parameters: image_size=small|large
```

#### 5. 格式转换

获取 SDF（结构数据文件）形式的化合物：
```
GET /rest/pug/compound/cid/{cid}/SDF
```

得到化合物 MOL：
```
GET /rest/pug/compound/cid/{cid}/record/SDF
```

#### 6.同义词检索

获取化合物的所有同义词：
```
GET /rest/pug/compound/cid/{cid}/synonyms/JSON
```

#### 7. 生物测定数据

获取化合物的生物测定数据：
```
GET /rest/pug/compound/cid/{cid}/assaysummary/JSON
```

获取具体的检测信息：
```
GET /rest/pug/assay/aid/{aid}/description/JSON
```

### 异步请求

对于大型查询（相似性/子结构搜索），PUG-REST 使用异步模式：

1.提交查询（返回ListKey）
2. 使用ListKey检查状态
3. 准备好后检索结果

工作流程示例：
```python
# Step 1: Submit similarity search
response = requests.post(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/similarity/smiles/{smiles}/cids/JSON",
    data={"Threshold": 90}
)
listkey = response.json()["Waiting"]["ListKey"]

# Step 2: Check status
status_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/{listkey}/cids/JSON"

# Step 3: Poll until ready (with timeout)
# Step 4: Retrieve results from the same URL
```

### 使用限制

**速率限制**：
- 每秒最多 5 个请求
- 每分钟最多 400 个请求
- 每分钟最长运行时间 300 秒

**最佳实践**：
- 尽可能使用批量请求
- 实施指数退避重试
- 适当时缓存结果
- 对大型查询使用异步模式

## PubChemPy Python 库

PubChemPy 是一个 Python 包装器，可简化 PUG-REST API 访问。

### 安装

```bash
pip install pubchempy
```

### 重点课程

#### 复合类

表示化合物的主类：

```python
import pubchempy as pcp

# Get by CID
compound = pcp.Compound.from_cid(2244)

# Access properties
compound.molecular_formula  # 'C9H8O4'
compound.molecular_weight   # 180.16
compound.iupac_name        # '2-acetyloxybenzoic acid'
compound.canonical_smiles   # 'CC(=O)OC1=CC=CC=C1C(=O)O'
compound.isomeric_smiles    # Same as canonical for non-stereoisomers
compound.inchi             # InChI string
compound.inchikey          # InChI Key
compound.xlogp             # Partition coefficient
compound.tpsa              # Topological polar surface area
```

#### 搜索方法

**按姓名**：
```python
compounds = pcp.get_compounds('aspirin', 'name')
# Returns list of Compound objects
```

**通过微笑**：
```python
compound = pcp.get_compounds('CC(=O)OC1=CC=CC=C1C(=O)O', 'smiles')[0]
```

**作者：InChI**：
```python
compound = pcp.get_compounds('InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)', 'inchi')[0]
```

**按公式**：
```python
compounds = pcp.get_compounds('C9H8O4', 'formula')
# Returns all compounds with this formula
```

**相似性搜索**：
```python
results = pcp.get_compounds('CC(=O)OC1=CC=CC=C1C(=O)O', 'smiles',
                           searchtype='similarity',
                           Threshold=90)
```

**子结构搜索**：
```python
results = pcp.get_compounds('c1ccccc1', 'smiles',
                           searchtype='substructure')
# Returns all compounds containing benzene ring
```

#### 财产检索

获取多种化合物的特定属性：
```python
properties = pcp.get_properties(
    ['MolecularFormula', 'MolecularWeight', 'CanonicalSMILES'],
    'aspirin',
    'name'
)
# Returns list of dictionaries
```

获取 pandas DataFrame 的属性：
```python
import pandas as pd
df = pd.DataFrame(properties)
```

#### 同义词

获取化合物的所有同义词：
```python
synonyms = pcp.get_synonyms('aspirin', 'name')
# Returns list of dictionaries with CID and synonym lists
```

#### 下载格式

下载各种格式的化合物：
```python
# Get as SDF
sdf_data = pcp.download('SDF', 'aspirin', 'name', overwrite=True)

# Get as JSON
json_data = pcp.download('JSON', '2244', 'cid')

# Get as PNG image
pcp.download('PNG', '2244', 'cid', 'aspirin.png', overwrite=True)
```

### 错误处理

```python
from pubchempy import BadRequestError, NotFoundError, TimeoutError

try:
    compound = pcp.get_compounds('nonexistent', 'name')
except NotFoundError:
    print("Compound not found")
except BadRequestError:
    print("Invalid request")
except TimeoutError:
    print("Request timed out")
```

## PUG-查看 API

PUG-View 提供对完整文本注释和专业报告的访问。

### 关键端点

获取复合注释：
```
GET /rest/pug_view/data/compound/{cid}/JSON
```

获取具体注释部分：
```
GET /rest/pug_view/data/compound/{cid}/JSON?heading={section_name}
```

可用部分包括：
- 化学和物理特性
- 药物和药物信息
- 药理学和生物化学
- 安全和危险
- 毒性
- 文学
- 专利
- 生物分子相互作用和途径

## 常见工作流程
### 1. 化学标识符转换

从名称转换为 SMILES 到 InChI：
```python
import pubchempy as pcp

compound = pcp.get_compounds('caffeine', 'name')[0]
smiles = compound.canonical_smiles
inchi = compound.inchi
inchikey = compound.inchikey
cid = compound.cid
```

### 2.批量属性检索

获取多种化合物的属性：
```python
compound_names = ['aspirin', 'ibuprofen', 'paracetamol']
properties = []

for name in compound_names:
    props = pcp.get_properties(
        ['MolecularFormula', 'MolecularWeight', 'XLogP'],
        name,
        'name'
    )
    properties.extend(props)

import pandas as pd
df = pd.DataFrame(properties)
```

### 3. 寻找相似的化合物

查找与查询结构相似的化合物：
```python
# Start with a known compound
query_compound = pcp.get_compounds('gefitinib', 'name')[0]
query_smiles = query_compound.canonical_smiles

# Perform similarity search
similar = pcp.get_compounds(
    query_smiles,
    'smiles',
    searchtype='similarity',
    Threshold=85
)

# Get properties for similar compounds
for compound in similar[:10]:  # First 10 results
    print(f"{compound.cid}: {compound.iupac_name}, MW: {compound.molecular_weight}")
```

### 4. 下部结构筛选

查找包含特定子结构的所有化合物：
```python
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

### 5.生物活性数据检索

```python
import requests

cid = 2244  # Aspirin
url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"

response = requests.get(url)
if response.status_code == 200:
    bioassay_data = response.json()
    # Process bioassay information
```

## 提示和最佳实践

1. **使用CID进行重复查询**：CID比名称或结构更高效
2. **缓存结果**：将经常访问的数据存储在本地
3. **批量请求**：尽可能组合多个查询
4. **处理速率限制**：在请求之间实现延迟
5. **使用适当的搜索类型**：相关化合物的相似性、基序查找的子结构
6. **利用 PubChemPy**：更高级别的抽象简化了常见任务
7. **处理缺失数据**：并非所有属性都适用于所有化合物
8. **使用异步模式**：用于大量相似性/子结构搜索
9. **指定输出格式**：选择 JSON 进行编程访问，选择 SDF 进行化学信息学工具
10. **阅读文档**：完整的 PUG-REST 文档位于 https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest

## 其他资源

- PubChem 主页：https://pubchem.ncbi.nlm.nih.gov/
- PUG-REST 文档：https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
- PUG-REST 教程：https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest-tutorial
- PubChemPy 文档：https://pubchempy.readthedocs.io/
- PubChemPy GitHub：https://github.com/mcs07/PubChemPy
- IUPAC 教程：https://iupac.github.io/WFChemCookbook/datasources/pubchem_pugrest.html