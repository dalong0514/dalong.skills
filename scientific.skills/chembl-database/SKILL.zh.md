<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：chembl 数据库
描述：“查询 ChEMBL 的生物活性分子和药物发现数据。按结构/性质搜索化合物，检索生物活性数据（IC50、Ki）、查找抑制剂、执行 SAR 研究，用于药物化学。”
---

# ChEMBL 数据库

## 概述

ChEMBL 是由欧洲生物信息学研究所 (EBI) 维护的手动管理的生物活性分子数据库，包含超过 200 万种化合物、1900 万个生物活性测量值、13,000 多个药物靶点以及已批准药物和临床候选药物的数据。使用 ChEMBL Python 客户端以编程方式访问和查询这些数据，以进行药物发现和药物化学研究。

## 何时使用此技能

该技能应该在以下情况下使用：

- **化合物搜索**：按名称、结构或属性查找分子
- **目标信息**：检索有关蛋白质、酶或生物目标的数据
- **生物活性数据**：查询 IC50、Ki、EC50 或其他活性测量值
- **药物信息**：查找批准的药物、机制或适应症
- **结构搜索**：执行相似性或子结构搜索
- **化学信息学**：分析分子特性和药物相似性
- **靶标-配体关系**：探索化合物-靶标相互作用
- **药物发现**：识别抑制剂、激动剂或生物活性分子

## 安装和设置

### Python 客户端

编程访问需要 ChEMBL Python 客户端：

```bash
uv pip install chembl_webresource_client
```

### 基本使用模式

<<<代码块_1>>>

## 核心能力

### 1. 分子查询

**通过 ChEMBL ID 检索：**
<<<代码块_2>>>

**按名称搜索：**
<<<代码块_3>>>

**按属性过滤：**
<<<代码块_4>>>

### 2. 目标查询

**检索目标信息：**
<<<代码块_5>>>

**搜索特定目标类型：**
<<<代码块_6>>>

### 3. 生物活性数据

**查询目标的活动：**
```python
activity = new_client.activity
# Find potent EGFR inhibitors
results = activity.filter(
    target_chembl_id='CHEMBL203',
    standard_type='IC50',
    standard_value__lte=100,
    standard_units='nM'
)
```

**获取化合物的所有活性：**
```python
compound_activities = activity.filter(
    molecule_chembl_id='CHEMBL25',
    pchembl_value__isnull=False
)
```

### 4. 基于结构的搜索

**相似性搜索：**
```python
similarity = new_client.similarity
# Find compounds similar to aspirin
similar = similarity.filter(
    smiles='CC(=O)Oc1ccccc1C(=O)O',
    similarity=85  # 85% similarity threshold
)
```

**子结构搜索：**
```python
substructure = new_client.substructure
# Find compounds containing benzene ring
results = substructure.filter(smiles='c1ccccc1')
```

### 5. 药品信息

**检索药品数据：**
```python
drug = new_client.drug
drug_info = drug.get('CHEMBL25')
```

**获取作用机制：**
```python
mechanism = new_client.mechanism
mechanisms = mechanism.filter(molecule_chembl_id='CHEMBL25')
```

**查询药品适应症：**
```python
drug_indication = new_client.drug_indication
indications = drug_indication.filter(molecule_chembl_id='CHEMBL25')
```

## 查询工作流程

### 工作流程 1：寻找目标抑制剂

1. **通过名称搜索来识别目标**：
   ```python
   targets = new_client.target.filter(pref_name__icontains='EGFR')
   target_id = targets[0]['target_chembl_id']
   ```

2. **查询该目标的生物活性数据**：
   ```python
   activities = new_client.activity.filter(
       target_chembl_id=target_id,
       standard_type='IC50',
       standard_value__lte=100
   )
   ```

3. **提取化合物 ID** 并检索详细信息：
   ```python
   compound_ids = [act['molecule_chembl_id'] for act in activities]
   compounds = [new_client.molecule.get(cid) for cid in compound_ids]
   ```

### 工作流程 2：分析已知药物

1. **获取药品信息**：
   ```python
   drug_info = new_client.drug.get('CHEMBL1234')
   ```

2. **检索机制**：
   ```python
   mechanisms = new_client.mechanism.filter(molecule_chembl_id='CHEMBL1234')
   ```

3. **查找所有生物活性**：
   ```python
   activities = new_client.activity.filter(molecule_chembl_id='CHEMBL1234')
   ```

### 工作流程 3：构效关系 (SAR) 研究

1. **找到相似的化合物**：
   ```python
   similar = new_client.similarity.filter(smiles='query_smiles', similarity=80)
   ```

2. **获取每个化合物的活性**：
   ```python
   for compound in similar:
       activities = new_client.activity.filter(
           molecule_chembl_id=compound['molecule_chembl_id']
       )
   ```

3. **利用结果中的分子特性分析特性-活性关系**。

## 过滤运算符

ChEMBL 支持 Django 风格的查询过滤器：

- `__exact` - 完全匹配
- `__iexact` - 不区分大小写的精确匹配
- `__contains` / `__icontains` - 子字符串匹配
- `__startswith` / `__endswith` - 前缀/后缀匹配
- `__gt`、`__gte`、`__lt`、`__lte` - 数字比较
- `__range` - 范围内的值
- `__in` - 列表中的值
- `__isnull` - 空/非空检查

## 数据导出与分析

将结果转换为 pandas DataFrame 进行分析：

```python
import pandas as pd

activities = new_client.activity.filter(target_chembl_id='CHEMBL203')
df = pd.DataFrame(list(activities))

# Analyze results
print(df['standard_value'].describe())
print(df.groupby('standard_type').size())
```

## 性能优化

### 缓存

客户端自动缓存结果24小时。配置缓存：

```python
from chembl_webresource_client.settings import Settings

# Disable caching
Settings.Instance().CACHING = False

# Adjust cache expiration (seconds)
Settings.Instance().CACHE_EXPIRE = 86400
```

### 惰性评估

仅当访问数据时才执行查询。转换为列表强制执行：

```python
# Query is not executed yet
results = molecule.filter(pref_name__icontains='aspirin')

# Force execution
results_list = list(results)
```

### 分页

结果自动分页。迭代所有结果：

```python
for activity in new_client.activity.filter(target_chembl_id='CHEMBL203'):
    # Process each activity
    print(activity['molecule_chembl_id'])
```

## 常见用例

### 寻找激酶抑制剂

```python
# Identify kinase targets
kinases = new_client.target.filter(
    target_type='SINGLE PROTEIN',
    pref_name__icontains='kinase'
)

# Get potent inhibitors
for kinase in kinases[:5]:  # First 5 kinases
    activities = new_client.activity.filter(
        target_chembl_id=kinase['target_chembl_id'],
        standard_type='IC50',
        standard_value__lte=50
    )
```

### 探索药物再利用

```python
# Get approved drugs
drugs = new_client.drug.filter()

# For each drug, find all targets
for drug in drugs[:10]:
    mechanisms = new_client.mechanism.filter(
        molecule_chembl_id=drug['molecule_chembl_id']
    )
```

### 虚拟筛选

```python
# Find compounds with desired properties
candidates = new_client.molecule.filter(
    molecule_properties__mw_freebase__range=[300, 500],
    molecule_properties__alogp__lte=5,
    molecule_properties__hba__lte=10,
    molecule_properties__hbd__lte=5
)
```

## 资源

### 脚本/example_queries.py

演示常见 ChEMBL 查询模式的即用型 Python 函数：

- `get_molecule_info()` - 通过 ID 检索分子详细信息
- `search_molecules_by_name()` - 基于名称的分子搜索
- `find_molecules_by_properties()` - 基于属性的过滤
- `get_bioactivity_data()` - 查询目标的生物活性
- `find_similar_compounds()` - 相似性搜索
- `substructure_search()` - 子结构匹配
- `get_drug_info()` - 检索药品信息
- `find_kinase_inhibitors()` - 专门的激酶抑制剂搜索
- `export_to_dataframe()` - 将结果转换为 pandas DataFrame

请参阅此脚本以获取实现细节和使用示例。

### 参考资料/api_reference.md

全面的 API 文档包括：

- 完整的终点列表（分子、靶标、活性、测定、药物等）
- 所有过滤运算符和查询模式
- 分子特性和生物活性领域
- 高级查询示例
- 配置和性能调整
- 错误处理和速率限制

当需要详细的 API 信息或对查询进行故障排除时，请参阅此文档。

## 重要提示

### 数据可靠性

- ChEMBL 数据是手动整理的，但可能包含不一致的内容
- 始终检查活动记录中的 `data_validity_comment` 字段
- 注意 `potential_duplicate` 标志

### 单位和标准

- 生物活性值使用标准单位（nM、uM等）
- `pchembl_value` 提供标准化活动（-对数刻度）
- 检查`standard_type`以了解测量类型（IC50、Ki、EC50等）

### 速率限制

- 尊重 ChEMBL 的公平使用政策
- 使用缓存来最大程度地减少重复请求
- 考虑批量下载大型数据集
- 避免用快速连续的请求来打击 API

### 化学结构格式

- SMILES 字符串是主要结构格式
- InChI 键可用于化合物
- 可以通过图像端点生成SVG图像

## 其他资源

- ChEMBL 网站：https://www.ebi.ac.uk/chembl/
- API 文档：https://www.ebi.ac.uk/chembl/api/data/docs
- Python 客户端 GitHub：https://github.com/chembl/chembl_webresource_client
- 接口文档：https://chembl.gitbook.io/chembl-interface-documentation/
- 示例笔记本：https://github.com/chembl/notebooks