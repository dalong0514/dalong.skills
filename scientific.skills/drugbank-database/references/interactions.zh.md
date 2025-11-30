<!-- 此文件由机器翻译自 interactions.md -->

# 药物间相互作用

## 概述
DrugBank 提供全面的药物相互作用 (DDI) 数据，包括机制、严重程度和临床意义。这些信息对于药物警戒、临床决策支持和药物安全研究至关重要。

## 交互数据结构

### XML 结构
```xml
<drug-interactions>
  <drug-interaction>
    <drugbank-id>DB00001</drugbank-id>
    <name>Warfarin</name>
    <description>The risk or severity of adverse effects can be increased...</description>
  </drug-interaction>
  <drug-interaction>
    <drugbank-id>DB00002</drugbank-id>
    <name>Aspirin</name>
    <description>May increase the anticoagulant activities...</description>
  </drug-interaction>
</drug-interactions>
```

### 交互组件
- **drugbank-id**：相互作用药物的 DrugBank ID
- **名称**：相互作用药物的名称
- **描述**：相互作用机制和临床意义的详细描述

## 提取药物相互作用

### 基本交互提取
<<<代码块_1>>>

### 双向交互映射
<<<代码块_2>>>

## 分析交互模式

### 计算每种药物的相互作用
<<<代码块_3>>>

### 寻找共同的互动伙伴
<<<代码块_4>>>

### 检查特定药物对
<<<代码块_5>>>

## 交互分类

### 解析交互描述
<<<代码块_6>>>

### 对交互进行分类
```python
def categorize_drug_interactions(drugbank_id):
    """Categorize interactions by severity and mechanism"""
    interactions = get_drug_interactions(drugbank_id)

    categorized = {
        'major': [],
        'moderate': [],
        'minor': [],
        'unknown': []
    }

    for interaction in interactions:
        severity = classify_interaction_severity(interaction['description'])
        interaction['severity'] = severity
        interaction['mechanisms'] = classify_interaction_mechanism(interaction['description'])
        categorized[severity].append(interaction)

    return categorized

# Usage
categorized = categorize_drug_interactions('DB00001')
print(f"Major: {len(categorized['major'])}")
print(f"Moderate: {len(categorized['moderate'])}")
print(f"Minor: {len(categorized['minor'])}")
```

## 构建交互矩阵

### 创建配对交互矩阵
```python
import pandas as pd
import numpy as np

def create_interaction_matrix(drug_ids):
    """Create binary interaction matrix for specified drugs"""
    n = len(drug_ids)
    matrix = np.zeros((n, n), dtype=int)

    # Build index mapping
    id_to_idx = {drug_id: idx for idx, drug_id in enumerate(drug_ids)}

    # Fill matrix
    for i, drug_id in enumerate(drug_ids):
        interactions = get_drug_interactions(drug_id)
        for interaction in interactions:
            partner_id = interaction['partner_id']
            if partner_id in id_to_idx:
                j = id_to_idx[partner_id]
                matrix[i, j] = 1
                matrix[j, i] = 1  # Symmetric

    df = pd.DataFrame(matrix, index=drug_ids, columns=drug_ids)
    return df

# Example: Create matrix for top 100 drugs
top_100_drugs = [drug['id'] for drug in rank_drugs_by_interactions()[:100]]
interaction_matrix = create_interaction_matrix(top_100_drugs)
```

### 导出互动网络
```python
def export_interaction_network_csv(output_file='drugbank_interactions.csv'):
    """Export all interactions as edge list (CSV)"""
    root = get_drugbank_root()
    ns = {'db': 'http://www.drugbank.ca'}

    edges = []

    for drug in root.findall('db:drug', ns):
        drug_id = drug.find('db:drugbank-id[@primary="true"]', ns).text
        drug_name = drug.find('db:name', ns).text

        ddi_elem = drug.find('db:drug-interactions', ns)
        if ddi_elem is not None:
            for interaction in ddi_elem.findall('db:drug-interaction', ns):
                partner_id = interaction.find('db:drugbank-id', ns).text
                partner_name = interaction.find('db:name', ns).text
                description = interaction.find('db:description', ns).text

                edges.append({
                    'drug1_id': drug_id,
                    'drug1_name': drug_name,
                    'drug2_id': partner_id,
                    'drug2_name': partner_name,
                    'description': description
                })

    df = pd.DataFrame(edges)
    df.to_csv(output_file, index=False)
    print(f"Exported {len(edges)} interactions to {output_file}")

# Usage
export_interaction_network_csv()
```

## 网络分析

### 图形表示
```python
import networkx as nx

def build_interaction_graph():
    """Build NetworkX graph of drug interactions"""
    network = build_interaction_network()

    G = nx.Graph()

    # Add nodes and edges
    for drug_id, partners in network.items():
        G.add_node(drug_id)
        for partner_id in partners:
            G.add_edge(drug_id, partner_id)

    return G

# Build graph
G = build_interaction_graph()
print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

# Network statistics
density = nx.density(G)
print(f"Network density: {density:.4f}")

# Find highly connected drugs (hubs)
degree_dict = dict(G.degree())
top_hubs = sorted(degree_dict.items(), key=lambda x: x[1], reverse=True)[:10]
print("Top 10 hubs:", top_hubs)
```

### 社区检测
```python
def detect_interaction_communities():
    """Detect communities in interaction network"""
    G = build_interaction_graph()

    # Louvain community detection
    from networkx.algorithms import community
    communities = community.louvain_communities(G)

    print(f"Detected {len(communities)} communities")

    # Analyze communities
    for i, comm in enumerate(communities[:5], 1):  # Top 5 communities
        print(f"Community {i}: {len(comm)} drugs")

    return communities

# Usage
communities = detect_interaction_communities()
```

## 临床应用

### 多药分析
```python
def check_polypharmacy_interactions(drug_list):
    """Check for interactions in a drug regimen"""
    print(f"Checking interactions for {len(drug_list)} drugs...")

    all_interactions = []

    # Check all pairs
    for i, drug1 in enumerate(drug_list):
        for drug2 in drug_list[i+1:]:
            interaction = check_interaction(drug1, drug2)
            if interaction:
                interaction['drug1'] = drug1
                interaction['drug2'] = drug2
                all_interactions.append(interaction)

    return all_interactions

# Example: Check patient drug regimen
patient_drugs = ['DB00001', 'DB00002', 'DB00005', 'DB00009']
interactions = check_polypharmacy_interactions(patient_drugs)

print(f"\nFound {len(interactions)} interactions:")
for interaction in interactions:
    print(f"\n{interaction['drug1']} + {interaction['drug2']}")
    print(f"  {interaction['description'][:100]}...")
```

### 互动风险评分
```python
def calculate_interaction_risk_score(drug_list):
    """Calculate overall interaction risk for drug combination"""
    interactions = check_polypharmacy_interactions(drug_list)

    severity_weights = {'major': 3, 'moderate': 2, 'minor': 1, 'unknown': 1}

    total_score = 0
    for interaction in interactions:
        severity = classify_interaction_severity(interaction['description'])
        total_score += severity_weights[severity]

    return {
        'total_interactions': len(interactions),
        'risk_score': total_score,
        'average_severity': total_score / len(interactions) if interactions else 0
    }

# Usage
risk = calculate_interaction_risk_score(patient_drugs)
print(f"Risk Score: {risk['risk_score']}, Avg Severity: {risk['average_severity']:.2f}")
```

## 最佳实践

1. **双向检查**：始终检查两个方向的交互（A→B 和 B→A）
2. **背景很重要**：在解释相互作用的意义时考虑临床背景
3. **最新数据**：使用最新的 DrugBank 版本获取最新的相互作用数据
4. **严重程度分类**：根据您的临床需求实施自定义分类
5. **网络分析**：使用图形分析来识别高风险药物组合
6. **临床验证**：与临床指南和文献的交叉引用
7. **文档**：记录 DrugBank 版本和分析方法的重现性