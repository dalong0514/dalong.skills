<!-- 此文件由机器翻译自 workflows.md -->

# gget 工作流程示例

扩展工作流程示例演示如何组合多个 gget 模块来执行常见的生物信息学任务。

## 目录
1.【完整的基因分析流程】(#complete-gene-analysis-pipeline)
2.【比较结构生物学】(#comparative-structural-biology)
3. [癌症基因组学分析](#cancer-genomics-analysis)
4. [单细胞表达分析](#single-cell-express-analysis)
5. [构建参考转录组](#building-reference-transcriptomes)
6.【突变影响评估】(#mutation-impact-assessment)
7. [药物靶点发现](#drug-target-discovery)

---

## 完整的基因分析流程

对基因从发现到功能注释的全面分析。

```python
import gget
import pandas as pd

# Step 1: Search for genes of interest
print("Step 1: Searching for GABA receptor genes...")
search_results = gget.search(["GABA", "receptor", "alpha"],
                             species="homo_sapiens",
                             andor="and")
print(f"Found {len(search_results)} genes")

# Step 2: Get detailed information
print("\nStep 2: Getting detailed information...")
gene_ids = search_results["ensembl_id"].tolist()[:5]  # Top 5 genes
gene_info = gget.info(gene_ids, pdb=True)
print(gene_info[["ensembl_id", "gene_name", "uniprot_id", "description"]])

# Step 3: Retrieve sequences
print("\nStep 3: Retrieving sequences...")
nucleotide_seqs = gget.seq(gene_ids)
protein_seqs = gget.seq(gene_ids, translate=True)

# Save sequences
with open("gaba_receptors_nt.fasta", "w") as f:
    f.write(nucleotide_seqs)
with open("gaba_receptors_aa.fasta", "w") as f:
    f.write(protein_seqs)

# Step 4: Get expression data
print("\nStep 4: Getting tissue expression...")
for gene_id, gene_name in zip(gene_ids, gene_info["gene_name"]):
    expr_data = gget.archs4(gene_name, which="tissue")
    print(f"\n{gene_name} expression:")
    print(expr_data.head())

# Step 5: Find correlated genes
print("\nStep 5: Finding correlated genes...")
correlated = gget.archs4(gene_info["gene_name"].iloc[0], which="correlation")
correlated_top = correlated.head(20)
print(correlated_top)

# Step 6: Enrichment analysis on correlated genes
print("\nStep 6: Performing enrichment analysis...")
gene_list = correlated_top["gene_symbol"].tolist()
enrichment = gget.enrichr(gene_list, database="ontology", plot=True)
print(enrichment.head(10))

# Step 7: Get disease associations
print("\nStep 7: Getting disease associations...")
for gene_id, gene_name in zip(gene_ids[:3], gene_info["gene_name"][:3]):
    diseases = gget.opentargets(gene_id, resource="diseases", limit=5)
    print(f"\n{gene_name} disease associations:")
    print(diseases)

# Step 8: Check for orthologs
print("\nStep 8: Finding orthologs...")
orthologs = gget.bgee(gene_ids[0], type="orthologs")
print(orthologs)

print("\nComplete gene analysis pipeline finished!")
```

---

## 比较结构生物学

比较不同物种的蛋白质结构并分析功能基序。

<<<代码块_1>>>

---

## 癌症基因组学分析

分析癌症相关基因及其突变。

<<<代码块_2>>>

---

## 单细胞表达分析

分析特定细胞类型和组织的单细胞 RNA-seq 数据。

<<<代码块_3>>>

---

## 构建参考转录组

为 RNA-seq 分析流程准备参考数据。

<<<代码块_4>>>

<<<代码块_5>>>

---

## 突变影响评估

分析基因突变对蛋白质结构和功能的影响。

<<<代码块_6>>>

---

## 药物靶点发现

识别并验证特定疾病的潜在药物靶点。

```python
import gget
import pandas as pd

print("Drug Target Discovery Workflow")
print("=" * 50)

# Step 1: Search for disease-related genes
disease = "alzheimer"
print(f"\n1. Searching for {disease} disease genes...")
genes = gget.search([disease], species="homo_sapiens", limit=50)
print(f"Found {len(genes)} potential genes")

# Step 2: Get detailed information
print("\n2. Getting detailed gene information...")
gene_ids = genes["ensembl_id"].tolist()[:20]  # Top 20
gene_info = gget.info(gene_ids[:10])  # Limit to avoid timeout

# Step 3: Get disease associations from OpenTargets
print("\n3. Getting disease associations...")
disease_scores = []
for gene_id, gene_name in zip(gene_info["ensembl_id"], gene_info["gene_name"]):
    diseases = gget.opentargets(gene_id, resource="diseases", limit=10)

    # Filter for Alzheimer's disease
    alzheimer = diseases[diseases["disease_name"].str.contains("Alzheimer", case=False, na=False)]

    if len(alzheimer) > 0:
        disease_scores.append({
            "ensembl_id": gene_id,
            "gene_name": gene_name,
            "disease_score": alzheimer["overall_score"].max()
        })

disease_df = pd.DataFrame(disease_scores).sort_values("disease_score", ascending=False)
print("\nTop disease-associated genes:")
print(disease_df.head(10))

# Step 4: Get tractability information
print("\n4. Assessing target tractability...")
top_targets = disease_df.head(5)
for _, row in top_targets.iterrows():
    tractability = gget.opentargets(
        row["ensembl_id"],
        resource="tractability"
    )
    print(f"\n{row['gene_name']} tractability:")
    print(tractability)

# Step 5: Get expression data
print("\n5. Getting tissue expression data...")
for _, row in top_targets.iterrows():
    # Brain expression from OpenTargets
    expression = gget.opentargets(
        row["ensembl_id"],
        resource="expression",
        filter_tissue="brain"
    )
    print(f"\n{row['gene_name']} brain expression:")
    print(expression)

    # Tissue expression from ARCHS4
    tissue_expr = gget.archs4(row["gene_name"], which="tissue")
    brain_expr = tissue_expr[tissue_expr["tissue"].str.contains("brain", case=False, na=False)]
    print(f"ARCHS4 brain expression:")
    print(brain_expr)

# Step 6: Check for existing drugs
print("\n6. Checking for existing drugs...")
for _, row in top_targets.iterrows():
    drugs = gget.opentargets(row["ensembl_id"], resource="drugs", limit=5)
    print(f"\n{row['gene_name']} drug associations:")
    if len(drugs) > 0:
        print(drugs[["drug_name", "drug_type", "max_phase_for_all_diseases"]])
    else:
        print("No drugs found")

# Step 7: Get protein-protein interactions
print("\n7. Getting protein-protein interactions...")
for _, row in top_targets.iterrows():
    interactions = gget.opentargets(
        row["ensembl_id"],
        resource="interactions",
        limit=10
    )
    print(f"\n{row['gene_name']} interacts with:")
    if len(interactions) > 0:
        print(interactions[["gene_b_symbol", "interaction_score"]])

# Step 8: Enrichment analysis
print("\n8. Performing pathway enrichment...")
gene_list = top_targets["gene_name"].tolist()
enrichment = gget.enrichr(gene_list, database="pathway", plot=True)
print("\nTop enriched pathways:")
print(enrichment.head(10))

# Step 9: Get structure information
print("\n9. Getting structure information...")
for _, row in top_targets.iterrows():
    info = gget.info([row["ensembl_id"]], pdb=True)

    if "pdb_id" in info.columns and pd.notna(info["pdb_id"].iloc[0]):
        pdb_ids = info["pdb_id"].iloc[0].split(";")
        print(f"\n{row['gene_name']} PDB structures: {', '.join(pdb_ids[:3])}")
    else:
        print(f"\n{row['gene_name']}: No PDB structure available")
        # Could predict with AlphaFold
        print(f"  Consider AlphaFold prediction")

# Step 10: Generate target summary report
print("\n10. Generating target summary report...")
report = []
for _, row in top_targets.iterrows():
    report.append({
        "Gene": row["gene_name"],
        "Ensembl ID": row["ensembl_id"],
        "Disease Score": row["disease_score"],
        "Target Status": "High Priority"
    })

report_df = pd.DataFrame(report)
report_df.to_csv("drug_targets_report.csv", index=False)
print("\nTarget report saved to drug_targets_report.csv")

print("\nDrug target discovery workflow completed!")
```

---

## 工作流程开发技巧

### 错误处理
```python
import gget

def safe_gget_call(func, *args, **kwargs):
    """Wrapper for gget calls with error handling"""
    try:
        result = func(*args, **kwargs)
        return result
    except Exception as e:
        print(f"Error in {func.__name__}: {str(e)}")
        return None

# Usage
result = safe_gget_call(gget.search, ["ACE2"], species="homo_sapiens")
if result is not None:
    print(result)
```

### 速率限制
```python
import time
import gget

def rate_limited_queries(gene_ids, delay=1):
    """Query multiple genes with rate limiting"""
    results = []
    for i, gene_id in enumerate(gene_ids):
        print(f"Querying {i+1}/{len(gene_ids)}: {gene_id}")
        result = gget.info([gene_id])
        results.append(result)

        if i < len(gene_ids) - 1:  # Don't sleep after last query
            time.sleep(delay)

    return pd.concat(results, ignore_index=True)
```

### 缓存结果
```python
import os
import pickle
import gget

def cached_gget(cache_file, func, *args, **kwargs):
    """Cache gget results to avoid repeated queries"""
    if os.path.exists(cache_file):
        print(f"Loading from cache: {cache_file}")
        with open(cache_file, "rb") as f:
            return pickle.load(f)

    result = func(*args, **kwargs)

    with open(cache_file, "wb") as f:
        pickle.dump(result, f)
    print(f"Saved to cache: {cache_file}")

    return result

# Usage
result = cached_gget("ace2_info.pkl", gget.info, ["ENSG00000130234"])
```

---

这些工作流程演示了如何组合多个 gget 模块进行全面的生物信息学分析。使它们适应您的具体研究问题和数据类型。