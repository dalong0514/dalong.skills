<!-- 此文件由机器翻译自 chemical-analysis.md -->

# 化学性质及相似性分析

## 概述
DrugBank 提供广泛的化学性质数据，包括分子结构、理化性质和计算描述符。该信息支持基于结构的分析、相似性搜索和 QSAR 建模。

## 化学标识符和结构

### 可用的结构格式
- **SMILES**：简化的分子输入线输入系统
- **InChI**：国际化学标识符
- **InChIKey**：用于数据库搜索的散列 InChI
- **分子式**：化学式（例如，C9H8O4）
- **IUPAC 名称**：系统化学名称
- **传统名称**：常用名称和同义词

### 提取化学结构
```python
from drugbank_downloader import get_drugbank_root

def get_drug_structures(drugbank_id):
    """Extract chemical structure representations"""
    root = get_drugbank_root()
    ns = {'db': 'http://www.drugbank.ca'}

    for drug in root.findall('db:drug', ns):
        primary_id = drug.find('db:drugbank-id[@primary="true"]', ns)
        if primary_id is not None and primary_id.text == drugbank_id:
            structures = {}

            # Get calculated properties
            calc_props = drug.find('db:calculated-properties', ns)
            if calc_props is not None:
                for prop in calc_props.findall('db:property', ns):
                    kind = prop.find('db:kind', ns).text
                    value = prop.find('db:value', ns).text

                    if kind in ['SMILES', 'InChI', 'InChIKey', 'Molecular Formula', 'IUPAC Name']:
                        structures[kind] = value

            return structures
    return {}

# Usage
structures = get_drug_structures('DB00001')
print(f"SMILES: {structures.get('SMILES')}")
print(f"InChI: {structures.get('InChI')}")
```

## 理化性质

### 计算属性
从结构计算的属性：
- **分子量**：精确质量（道尔顿）
- **logP**：分配系数（亲脂性）
- **logS**：水溶性
- **极表面积 (PSA)**：拓扑极表面积
- **氢键供体**：氢键供体的数量
- **H-Bond Acceptors**：氢键受体的数量
- **可旋转债券**：可旋转债券的数量
- **折射率**：摩尔折射率
- **极化率**：分子极化率

### 实验特性
文献中测量的特性：
- **熔点**：物理熔点
- **水溶性**：实验溶解度数据
- **pKa**：酸解离常数
- **疏水性**：实验 logP/logD 值

### 提取所有属性
<<<代码块_1>>>

## Lipinski 的五分析法则

### 五格规则
<<<代码块_2>>>

### Veber 的规则
<<<代码块_3>>>

## 化学相似性分析

### 与 RDKit 的基于结构的相似性
<<<代码块_4>>>

### 查找类似药物
<<<代码块_5>>>

### 批量相似度矩阵
<<<代码块_6>>>

## 分子指纹

### 生成不同的指纹类型
```python
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.Fingerprints import FingerprintMols

def generate_fingerprints(smiles):
    """Generate multiple types of molecular fingerprints"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    fingerprints = {
        'morgan_fp': AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048),
        'maccs_keys': MACCSkeys.GenMACCSKeys(mol),
        'topological_fp': FingerprintMols.FingerprintMol(mol),
        'atom_pairs': Pairs.GetAtomPairFingerprint(mol)
    }

    return fingerprints

# Generate fingerprints for a drug
structures = get_drug_structures('DB00001')
fps = generate_fingerprints(structures.get('SMILES'))
```

### 子结构搜索
```python
from rdkit.Chem import Fragments

def search_substructure(substructure_smarts):
    """Find drugs containing a specific substructure"""
    root = get_drugbank_root()
    ns = {'db': 'http://www.drugbank.ca'}

    pattern = Chem.MolFromSmarts(substructure_smarts)
    if pattern is None:
        print("Invalid SMARTS pattern")
        return []

    matching_drugs = []

    for drug in root.findall('db:drug', ns):
        drug_id = drug.find('db:drugbank-id[@primary="true"]', ns).text
        structures = get_drug_structures(drug_id)
        smiles = structures.get('SMILES')

        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol and mol.HasSubstructMatch(pattern):
                drug_name = drug.find('db:name', ns).text
                matching_drugs.append({
                    'drug_id': drug_id,
                    'drug_name': drug_name
                })

    return matching_drugs

# Example: Find drugs with benzene ring
benzene_drugs = search_substructure('c1ccccc1')
print(f"Found {len(benzene_drugs)} drugs with benzene ring")
```

## ADMET 属性预测

### 预测吸收
```python
def predict_oral_absorption(drugbank_id):
    """Predict oral absorption based on physicochemical properties"""
    props = get_all_properties(drugbank_id)
    calc_props = props.get('calculated', {})

    mw = float(calc_props.get('Molecular Weight', {}).get('value', 0))
    logp = float(calc_props.get('logP', {}).get('value', 0))
    psa = float(calc_props.get('Polar Surface Area (PSA)', {}).get('value', 0))
    h_donors = int(calc_props.get('H Bond Donor Count', {}).get('value', 0))

    # Simple absorption prediction
    good_absorption = (
        mw <= 500 and
        -0.5 <= logp <= 5.0 and
        psa <= 140 and
        h_donors <= 5
    )

    absorption_score = 0
    if mw <= 500:
        absorption_score += 25
    if -0.5 <= logp <= 5.0:
        absorption_score += 25
    if psa <= 140:
        absorption_score += 25
    if h_donors <= 5:
        absorption_score += 25

    return {
        'predicted_absorption': 'good' if good_absorption else 'poor',
        'absorption_score': absorption_score,
        'properties': {
            'molecular_weight': mw,
            'logP': logp,
            'psa': psa,
            'h_donors': h_donors
        }
    }
```

### BBB 渗透率预测
```python
def predict_bbb_permeability(drugbank_id):
    """Predict blood-brain barrier permeability"""
    props = get_all_properties(drugbank_id)
    calc_props = props.get('calculated', {})

    mw = float(calc_props.get('Molecular Weight', {}).get('value', 0))
    logp = float(calc_props.get('logP', {}).get('value', 0))
    psa = float(calc_props.get('Polar Surface Area (PSA)', {}).get('value', 0))
    h_donors = int(calc_props.get('H Bond Donor Count', {}).get('value', 0))

    # BBB permeability criteria (simplified)
    bbb_permeable = (
        mw <= 450 and
        logp <= 5.0 and
        psa <= 90 and
        h_donors <= 3
    )

    return {
        'bbb_permeable': bbb_permeable,
        'properties': {
            'molecular_weight': mw,
            'logP': logp,
            'psa': psa,
            'h_donors': h_donors
        }
    }
```

## 化学空间分析

### 主成分分析
```python
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def perform_chemical_space_pca(drug_ids):
    """Perform PCA on chemical descriptor space"""
    # Extract properties for all drugs
    properties_list = []
    valid_ids = []

    for drug_id in drug_ids:
        props = get_all_properties(drug_id)
        calc_props = props.get('calculated', {})

        try:
            prop_vector = [
                float(calc_props.get('Molecular Weight', {}).get('value', 0)),
                float(calc_props.get('logP', {}).get('value', 0)),
                float(calc_props.get('Polar Surface Area (PSA)', {}).get('value', 0)),
                int(calc_props.get('H Bond Donor Count', {}).get('value', 0)),
                int(calc_props.get('H Bond Acceptor Count', {}).get('value', 0)),
                int(calc_props.get('Rotatable Bond Count', {}).get('value', 0)),
            ]
            properties_list.append(prop_vector)
            valid_ids.append(drug_id)
        except (ValueError, TypeError):
            continue

    # Perform PCA
    X = np.array(properties_list)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    # Create DataFrame
    df = pd.DataFrame({
        'drug_id': valid_ids,
        'PC1': X_pca[:, 0],
        'PC2': X_pca[:, 1]
    })

    return df, pca

# Visualize chemical space
# drug_list = [all approved drugs]
# pca_df, pca_model = perform_chemical_space_pca(drug_list)
```

### 按化学性质聚类
```python
from sklearn.cluster import KMeans

def cluster_drugs_by_properties(drug_ids, n_clusters=10):
    """Cluster drugs based on chemical properties"""
    properties_list = []
    valid_ids = []

    for drug_id in drug_ids:
        props = get_all_properties(drug_id)
        calc_props = props.get('calculated', {})

        try:
            prop_vector = [
                float(calc_props.get('Molecular Weight', {}).get('value', 0)),
                float(calc_props.get('logP', {}).get('value', 0)),
                float(calc_props.get('Polar Surface Area (PSA)', {}).get('value', 0)),
                int(calc_props.get('H Bond Donor Count', {}).get('value', 0)),
                int(calc_props.get('H Bond Acceptor Count', {}).get('value', 0)),
            ]
            properties_list.append(prop_vector)
            valid_ids.append(drug_id)
        except (ValueError, TypeError):
            continue

    X = np.array(properties_list)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # K-means clustering
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(X_scaled)

    df = pd.DataFrame({
        'drug_id': valid_ids,
        'cluster': clusters
    })

    return df, kmeans
```

## 导出化学数据

### 创建化学性质数据库
```python
def export_chemical_properties(output_file='drugbank_chemical_properties.csv'):
    """Export all chemical properties to CSV"""
    root = get_drugbank_root()
    ns = {'db': 'http://www.drugbank.ca'}

    all_properties = []

    for drug in root.findall('db:drug', ns):
        drug_id = drug.find('db:drugbank-id[@primary="true"]', ns).text
        drug_name = drug.find('db:name', ns).text

        props = get_all_properties(drug_id)
        calc_props = props.get('calculated', {})

        property_dict = {
            'drug_id': drug_id,
            'drug_name': drug_name,
            'smiles': calc_props.get('SMILES', {}).get('value'),
            'inchi': calc_props.get('InChI', {}).get('value'),
            'inchikey': calc_props.get('InChIKey', {}).get('value'),
            'molecular_weight': calc_props.get('Molecular Weight', {}).get('value'),
            'logP': calc_props.get('logP', {}).get('value'),
            'psa': calc_props.get('Polar Surface Area (PSA)', {}).get('value'),
            'h_donors': calc_props.get('H Bond Donor Count', {}).get('value'),
            'h_acceptors': calc_props.get('H Bond Acceptor Count', {}).get('value'),
            'rotatable_bonds': calc_props.get('Rotatable Bond Count', {}).get('value'),
        }

        all_properties.append(property_dict)

    df = pd.DataFrame(all_properties)
    df.to_csv(output_file, index=False)
    print(f"Exported {len(all_properties)} drug properties to {output_file}")

# Usage
export_chemical_properties()
```

## 最佳实践

1. **结构验证**：在与 RDKit 一起使用之前始终验证 SMILES/InChI
2. **多个描述符**：使用多种指纹类型进行综合相似度
3. **阈值选择**：Tanimoto >0.85 = 非常相似，0.7-0.85 = 相似，<0.7 = 不同
4. **规则应用**：Lipinski 的 Ro5 和 Veber 的规则是指导方针，而不是绝对的界限
5. **ADMET预测**：使用计算预测作为筛选，通过实验验证
6. **化学空间**：可视化化学空间以了解药物多样性
7. **标准化**：比较前对分子进行标准化（中和、去除盐）
8. **性能**：缓存计算指纹以进行大规模相似性搜索