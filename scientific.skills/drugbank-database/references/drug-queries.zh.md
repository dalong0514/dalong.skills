<!-- 此文件由机器翻译自 drug-queries.md -->

# 药品信息查询

## 概述
DrugBank 提供全面的药物信息，每个条目有 200 多个数据字段，包括化学性质、药理学、作用机制和临床数据。

## 数据库内容

### 药物类别
- **FDA 批准的小分子**：~2,037 种药物
- **生物技术/生物药物**：约 241 个条目
- **营养保健品**：~96 种化合物
- **实验药物**：约 6,000 多种化合物
- **撤回/停产**：具有安全数据的历史药物

### 数据字段（每个条目 200 多个）
- **标识符**：DrugBank ID、CAS 编号、UNII、PubChem CID
- **名称**：通用名、品牌、同义词、IUPAC
- **化学**：结构（SMILES、InChI）、分子式、分子量
- **药理学**：适应症、作用机制、药效学
- **药代动力学**：吸收、分布、代谢、排泄 (ADME)
- **毒性**：LD50、不良反应、禁忌症
- **临床**：剂型、给药途径、半衰期
- **目标**：蛋白质、酶、转运蛋白、载体
- **相互作用**：药物-药物、药物-食物相互作用
- **参考文献**：引用文献和临床研究

## XML 结构导航

### 基本 XML 结构
```xml
<drugbank>
  <drug type="small molecule" created="..." updated="...">
    <drugbank-id primary="true">DB00001</drugbank-id>
    <name>Lepirudin</name>
    <description>...</description>
    <cas-number>...</cas-number>
    <synthesis-reference>...</synthesis-reference>
    <indication>...</indication>
    <pharmacodynamics>...</pharmacodynamics>
    <mechanism-of-action>...</mechanism-of-action>
    <toxicity>...</toxicity>
    <metabolism>...</metabolism>
    <absorption>...</absorption>
    <half-life>...</half-life>
    <protein-binding>...</protein-binding>
    <route-of-elimination>...</route-of-elimination>
    <calculated-properties>...</calculated-properties>
    <experimental-properties>...</experimental-properties>
    <targets>...</targets>
    <enzymes>...</enzymes>
    <transporters>...</transporters>
    <drug-interactions>...</drug-interactions>
  </drug>
</drugbank>
```

### 命名空间
DrugBank XML 使用名称空间。正确处理它们：
<<<代码块_1>>>

## 通过药品标识符查询

### 通过 DrugBank ID 查询
<<<代码块_2>>>

### 按名称查询
<<<代码块_3>>>

### 按CAS号查询
<<<代码块_4>>>

## 提取具体信息

### 药物基本信息
<<<代码块_5>>>

### 化学性质
<<<代码块_6>>>

### 药理学信息
```python
def extract_pharmacology(drug):
    """Extract pharmacological information"""
    ns = {'db': 'http://www.drugbank.ca'}

    pharm = {
        'indication': get_text_safe(drug.find('db:indication', ns)),
        'pharmacodynamics': get_text_safe(drug.find('db:pharmacodynamics', ns)),
        'mechanism_of_action': get_text_safe(drug.find('db:mechanism-of-action', ns)),
        'toxicity': get_text_safe(drug.find('db:toxicity', ns)),
        'metabolism': get_text_safe(drug.find('db:metabolism', ns)),
        'absorption': get_text_safe(drug.find('db:absorption', ns)),
        'half_life': get_text_safe(drug.find('db:half-life', ns)),
        'protein_binding': get_text_safe(drug.find('db:protein-binding', ns)),
        'route_of_elimination': get_text_safe(drug.find('db:route-of-elimination', ns)),
        'volume_of_distribution': get_text_safe(drug.find('db:volume-of-distribution', ns)),
        'clearance': get_text_safe(drug.find('db:clearance', ns)),
    }
    return pharm
```

### 外部标识符
```python
def extract_external_identifiers(drug):
    """Extract cross-references to other databases"""
    ns = {'db': 'http://www.drugbank.ca'}

    identifiers = {}

    external_ids = drug.find('db:external-identifiers', ns)
    if external_ids is not None:
        for ext_id in external_ids.findall('db:external-identifier', ns):
            resource = ext_id.find('db:resource', ns).text
            identifier = ext_id.find('db:identifier', ns).text
            identifiers[resource] = identifier

    return identifiers

# Common external databases:
# - PubChem Compound
# - PubChem Substance
# - ChEMBL
# - ChEBI
# - UniProtKB
# - KEGG Drug
# - PharmGKB
# - RxCUI (RxNorm)
# - ZINC
```

## 构建药物数据集

### 创建药物词典
```python
def build_drug_database():
    """Build searchable dictionary of all drugs"""
    root = get_drugbank_root()
    ns = {'db': 'http://www.drugbank.ca'}

    drug_db = {}

    for drug in root.findall('db:drug', ns):
        db_id = drug.find('db:drugbank-id[@primary="true"]', ns).text

        drug_info = {
            'id': db_id,
            'name': get_text_safe(drug.find('db:name', ns)),
            'type': drug.get('type'),
            'description': get_text_safe(drug.find('db:description', ns)),
            'cas': get_text_safe(drug.find('db:cas-number', ns)),
            'indication': get_text_safe(drug.find('db:indication', ns)),
        }

        drug_db[db_id] = drug_info

    return drug_db

# Create searchable database
drugs = build_drug_database()
print(f"Total drugs: {len(drugs)}")
```

### 导出到数据框
```python
import pandas as pd

def create_drug_dataframe():
    """Create pandas DataFrame of drug information"""
    root = get_drugbank_root()
    ns = {'db': 'http://www.drugbank.ca'}

    drugs_data = []

    for drug in root.findall('db:drug', ns):
        drug_dict = {
            'drugbank_id': drug.find('db:drugbank-id[@primary="true"]', ns).text,
            'name': get_text_safe(drug.find('db:name', ns)),
            'type': drug.get('type'),
            'cas_number': get_text_safe(drug.find('db:cas-number', ns)),
            'description': get_text_safe(drug.find('db:description', ns)),
            'indication': get_text_safe(drug.find('db:indication', ns)),
        }
        drugs_data.append(drug_dict)

    df = pd.DataFrame(drugs_data)
    return df

# Usage
df = create_drug_dataframe()
df.to_csv('drugbank_drugs.csv', index=False)
```

### 按药物类型过滤
```python
def filter_by_type(drug_type='small molecule'):
    """Get drugs of specific type"""
    root = get_drugbank_root()
    ns = {'db': 'http://www.drugbank.ca'}

    filtered_drugs = []

    for drug in root.findall('db:drug', ns):
        if drug.get('type') == drug_type:
            db_id = drug.find('db:drugbank-id[@primary="true"]', ns).text
            name = get_text_safe(drug.find('db:name', ns))
            filtered_drugs.append({'id': db_id, 'name': name})

    return filtered_drugs

# Get all biotech drugs
biotech_drugs = filter_by_type('biotech')
```

### 按关键字搜索
```python
def search_drugs_by_keyword(keyword, field='indication'):
    """Search drugs by keyword in specific field"""
    root = get_drugbank_root()
    ns = {'db': 'http://www.drugbank.ca'}

    results = []
    keyword_lower = keyword.lower()

    for drug in root.findall('db:drug', ns):
        field_elem = drug.find(f'db:{field}', ns)
        if field_elem is not None and field_elem.text:
            if keyword_lower in field_elem.text.lower():
                db_id = drug.find('db:drugbank-id[@primary="true"]', ns).text
                name = get_text_safe(drug.find('db:name', ns))
                results.append({
                    'id': db_id,
                    'name': name,
                    field: field_elem.text[:200]  # First 200 chars
                })

    return results

# Example: Find drugs for cancer treatment
cancer_drugs = search_drugs_by_keyword('cancer', 'indication')
```

## 性能优化

### 索引以加快查询速度
```python
def build_indexes():
    """Build indexes for faster lookups"""
    root = get_drugbank_root()
    ns = {'db': 'http://www.drugbank.ca'}

    # Index by ID, name, and CAS
    id_index = {}
    name_index = {}
    cas_index = {}

    for drug in root.findall('db:drug', ns):
        db_id = drug.find('db:drugbank-id[@primary="true"]', ns).text
        id_index[db_id] = drug

        name = get_text_safe(drug.find('db:name', ns))
        if name:
            name_index[name.lower()] = drug

        cas = get_text_safe(drug.find('db:cas-number', ns))
        if cas:
            cas_index[cas] = drug

    return {'id': id_index, 'name': name_index, 'cas': cas_index}

# Build once, query many times
indexes = build_indexes()
drug = indexes['name'].get('aspirin')
```