<!-- 此文件由机器翻译自 targets-pathways.md -->

# 药物靶点和途径

## 概述
DrugBank 提供有关药物-蛋白质相互作用的全面信息，包括靶标、酶、转运蛋白和载体。这些数据对于了解药物机制、确定重新利用机会和预测脱靶效应至关重要。

## 蛋白质相互作用类别

### 目标蛋白质
药物结合产生治疗效果的主要蛋白质：
- **受体**：G蛋白偶联受体、核受体、离子通道
- **酶**：激酶、蛋白酶、磷酸酶
- **传输器**：用作目标（不仅仅是 ADME）
- **其他**：结构蛋白、DNA/RNA

### 代谢酶
参与药物代谢的酶：
- **细胞色素 P450 酶**：CYP3A4、CYP2D6、CYP2C9 等。
- **第二阶段酶**：UGT、SULT、GST
- **酯酶和肽酶**

### 运输机
参与药物跨膜转运的蛋白质：
- **摄取转运蛋白**：OATP、OCT
- **外排转运蛋白**：P-糖蛋白、BCRP、MRP
- **其他**：SLC 和 ABC 运输机系列

### 运营商
结合和转运药物的血浆蛋白：
- **白蛋白**：血液中的主要药物载体
- **α-1-酸性糖蛋白**
- **脂蛋白**
- **特异性结合蛋白**：SHBG、CBG 等。

## XML 数据结构

### 目标元素结构
```xml
<targets>
  <target>
    <id>BE0000001</id>
    <name>Prothrombin</name>
    <organism>Humans</organism>
    <actions>
      <action>inhibitor</action>
    </actions>
    <known-action>yes</known-action>
    <polypeptide id="P00734" source="Swiss-Prot">
      <name>Prothrombin</name>
      <general-function>Serine-type endopeptidase activity</general-function>
      <specific-function>Thrombin plays a role in...</specific-function>
      <gene-name>F2</gene-name>
      <organism>Homo sapiens</organism>
      <external-identifiers>
        <external-identifier>
          <resource>UniProtKB</resource>
          <identifier>P00734</identifier>
        </external-identifier>
      </external-identifiers>
      <amino-acid-sequence>MAHVRGLQLP...</amino-acid-sequence>
      <pfams>...</pfams>
      <go-classifiers>...</go-classifiers>
    </polypeptide>
  </target>
</targets>
```

## 提取目标信息

### 获取药物靶点
<<<代码块_1>>>

### 获取所有蛋白质相互作用
<<<代码块_2>>>

## 建立目标药物网络

### 创建目标药物矩阵
<<<代码块_3>>>

### 寻找针对特定蛋白质的药物
<<<代码块_4>>>

### 寻找具有共同目标的药物
<<<代码块_5>>>

## 通路分析

### 提取途径信息
<<<代码块_6>>>

### 建立通路网络
```python
def build_pathway_drug_network():
    """Build network of pathways and drugs"""
    root = get_drugbank_root()
    ns = {'db': 'http://www.drugbank.ca'}

    pathway_network = {}

    for drug in root.findall('db:drug', ns):
        drug_id = drug.find('db:drugbank-id[@primary="true"]', ns).text

        pathways_elem = drug.find('db:pathways', ns)
        if pathways_elem is not None:
            for pathway in pathways_elem.findall('db:pathway', ns):
                pathway_id = pathway.find('db:smpdb-id', ns).text
                pathway_name = pathway.find('db:name', ns).text

                if pathway_id not in pathway_network:
                    pathway_network[pathway_id] = {
                        'name': pathway_name,
                        'drugs': []
                    }

                pathway_network[pathway_id]['drugs'].append(drug_id)

    return pathway_network
```

## 基于目标的药物再利用

### 寻找具有相似靶点的药物
```python
def find_similar_target_profiles(drugbank_id, min_shared_targets=2):
    """Find drugs with similar target profiles for repurposing"""
    reference_targets = get_drug_targets(drugbank_id)
    reference_target_ids = set(t.get('uniprot_id') or t['name'] for t in reference_targets)

    root = get_drugbank_root()
    ns = {'db': 'http://www.drugbank.ca'}

    similar_drugs = []

    for drug in root.findall('db:drug', ns):
        drug_id = drug.find('db:drugbank-id[@primary="true"]', ns).text

        if drug_id == drugbank_id:
            continue

        drug_targets = get_drug_targets(drug_id)
        drug_target_ids = set(t.get('uniprot_id') or t['name'] for t in drug_targets)

        shared = reference_target_ids.intersection(drug_target_ids)

        if len(shared) >= min_shared_targets:
            drug_name = drug.find('db:name', ns).text
            indication = get_text_safe(drug.find('db:indication', ns))

            similar_drugs.append({
                'drug_id': drug_id,
                'drug_name': drug_name,
                'shared_targets': len(shared),
                'total_targets': len(drug_target_ids),
                'overlap_ratio': len(shared) / len(drug_target_ids) if drug_target_ids else 0,
                'indication': indication,
                'shared_target_names': list(shared)
            })

    # Sort by overlap ratio
    similar_drugs.sort(key=lambda x: x['overlap_ratio'], reverse=True)
    return similar_drugs

# Example: Find repurposing candidates
candidates = find_similar_target_profiles('DB00001', min_shared_targets=2)
for drug in candidates[:5]:
    print(f"{drug['drug_name']}: {drug['shared_targets']} shared targets")
```

### 多药理分析
```python
def analyze_polypharmacology(drugbank_id):
    """Analyze on-target and off-target effects"""
    targets = get_drug_targets(drugbank_id)

    analysis = {
        'total_targets': len(targets),
        'known_action_targets': [],
        'unknown_action_targets': [],
        'target_classes': {},
        'organisms': {}
    }

    for target in targets:
        if target.get('known_action') == 'yes':
            analysis['known_action_targets'].append(target)
        else:
            analysis['unknown_action_targets'].append(target)

        # Count by organism
        organism = target.get('organism', 'Unknown')
        analysis['organisms'][organism] = analysis['organisms'].get(organism, 0) + 1

    return analysis

# Usage
poly_analysis = analyze_polypharmacology('DB00001')
print(f"Total targets: {poly_analysis['total_targets']}")
print(f"Known action: {len(poly_analysis['known_action_targets'])}")
print(f"Unknown action: {len(poly_analysis['unknown_action_targets'])}")
```

## 酶和转运蛋白分析

### CYP450 相互作用分析
```python
def analyze_cyp450_metabolism(drugbank_id):
    """Analyze CYP450 enzyme involvement"""
    interactions = get_all_protein_interactions(drugbank_id)
    enzymes = interactions['enzymes']

    cyp_enzymes = []
    for enzyme in enzymes:
        gene_name = enzyme.get('gene_name', '')
        if gene_name and gene_name.startswith('CYP'):
            cyp_enzymes.append({
                'gene': gene_name,
                'name': enzyme['name'],
                'actions': enzyme.get('actions', [])
            })

    return cyp_enzymes

# Check CYP involvement
cyp_data = analyze_cyp450_metabolism('DB00001')
for cyp in cyp_data:
    print(f"{cyp['gene']}: {cyp['actions']}")
```

### 转运蛋白底物分析
```python
def analyze_transporter_substrates(drugbank_id):
    """Identify transporter involvement for ADME"""
    interactions = get_all_protein_interactions(drugbank_id)
    transporters = interactions['transporters']

    transporter_info = {
        'efflux': [],
        'uptake': [],
        'other': []
    }

    for transporter in transporters:
        name = transporter['name'].lower()
        gene = transporter.get('gene_name', '').upper()

        if 'p-glycoprotein' in name or gene == 'ABCB1':
            transporter_info['efflux'].append(transporter)
        elif 'oatp' in name or 'slco' in gene.lower():
            transporter_info['uptake'].append(transporter)
        else:
            transporter_info['other'].append(transporter)

    return transporter_info
```

## GO 术语和蛋白质功能分析

### 提取 GO 术语
```python
def get_target_go_terms(drugbank_id):
    """Extract Gene Ontology terms for drug targets"""
    root = get_drugbank_root()
    ns = {'db': 'http://www.drugbank.ca'}

    for drug in root.findall('db:drug', ns):
        primary_id = drug.find('db:drugbank-id[@primary="true"]', ns)
        if primary_id is not None and primary_id.text == drugbank_id:
            go_terms = []

            targets_elem = drug.find('db:targets', ns)
            if targets_elem is not None:
                for target in targets_elem.findall('db:target', ns):
                    polypeptide = target.find('db:polypeptide', ns)
                    if polypeptide is not None:
                        go_classifiers = polypeptide.find('db:go-classifiers', ns)
                        if go_classifiers is not None:
                            for go_class in go_classifiers.findall('db:go-classifier', ns):
                                go_term = {
                                    'category': go_class.find('db:category', ns).text,
                                    'description': go_class.find('db:description', ns).text,
                                }
                                go_terms.append(go_term)

            return go_terms
    return []
```

## 最佳实践

1. **UniProt Cross-Reference**：使用 UniProt ID 进行跨数据库的准确蛋白质匹配
2. **作用分类**：注意作用类型（抑制剂、激动剂、拮抗剂等）
3. **已知与未知**：区分已验证的目标和预测/未知的相互作用
4. **生物体特异性**：分析目标数据时考虑生物体
5. **多药理学**：预测药物作用时考虑多个目标
6. **通路背景**：使用通路数据了解系统效应
7. **CYP450 分析**：对于预测药物间相互作用至关重要
8. **转运蛋白分析**：对于了解生物利用度和组织分布至关重要