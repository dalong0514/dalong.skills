<!-- 此文件由机器翻译自 workflows.md -->

# COBRApy 综合工作流程

本文档提供了代谢建模中常见 COBRApy 任务的详细分步工作流程。

## 工作流程 1：通过可视化完成淘汰赛研究

此工作流程演示了如何进行全面的基因敲除研究并可视化结果。

```python
import pandas as pd
import matplotlib.pyplot as plt
from cobra.io import load_model
from cobra.flux_analysis import single_gene_deletion, double_gene_deletion

# Step 1: Load model
model = load_model("ecoli")
print(f"Loaded model: {model.id}")
print(f"Model contains {len(model.reactions)} reactions, {len(model.metabolites)} metabolites, {len(model.genes)} genes")

# Step 2: Get baseline growth rate
baseline = model.slim_optimize()
print(f"Baseline growth rate: {baseline:.3f} /h")

# Step 3: Perform single gene deletions
print("Performing single gene deletions...")
single_results = single_gene_deletion(model)

# Step 4: Classify genes by impact
essential_genes = single_results[single_results["growth"] < 0.01]
severely_impaired = single_results[(single_results["growth"] >= 0.01) &
                                   (single_results["growth"] < 0.5 * baseline)]
moderately_impaired = single_results[(single_results["growth"] >= 0.5 * baseline) &
                                     (single_results["growth"] < 0.9 * baseline)]
neutral_genes = single_results[single_results["growth"] >= 0.9 * baseline]

print(f"\nSingle Deletion Results:")
print(f"  Essential genes: {len(essential_genes)}")
print(f"  Severely impaired: {len(severely_impaired)}")
print(f"  Moderately impaired: {len(moderately_impaired)}")
print(f"  Neutral genes: {len(neutral_genes)}")

# Step 5: Visualize distribution
fig, ax = plt.subplots(figsize=(10, 6))
single_results["growth"].hist(bins=50, ax=ax)
ax.axvline(baseline, color='r', linestyle='--', label='Baseline')
ax.set_xlabel("Growth rate (/h)")
ax.set_ylabel("Number of genes")
ax.set_title("Distribution of Growth Rates After Single Gene Deletions")
ax.legend()
plt.tight_layout()
plt.savefig("single_deletion_distribution.png", dpi=300)

# Step 6: Identify gene pairs for double deletions
# Focus on non-essential genes to find synthetic lethals
target_genes = single_results[single_results["growth"] >= 0.5 * baseline].index.tolist()
target_genes = [list(gene)[0] for gene in target_genes[:50]]  # Limit for performance

print(f"\nPerforming double deletions on {len(target_genes)} genes...")
double_results = double_gene_deletion(
    model,
    gene_list1=target_genes,
    processes=4
)

# Step 7: Find synthetic lethal pairs
synthetic_lethals = double_results[
    (double_results["growth"] < 0.01) &
    (single_results.loc[double_results.index.get_level_values(0)]["growth"].values >= 0.5 * baseline) &
    (single_results.loc[double_results.index.get_level_values(1)]["growth"].values >= 0.5 * baseline)
]

print(f"Found {len(synthetic_lethals)} synthetic lethal gene pairs")
print("\nTop 10 synthetic lethal pairs:")
print(synthetic_lethals.head(10))

# Step 8: Export results
single_results.to_csv("single_gene_deletions.csv")
double_results.to_csv("double_gene_deletions.csv")
synthetic_lethals.to_csv("synthetic_lethals.csv")
```

## 工作流程 2：媒体设计和优化

此工作流程展示了如何系统地设计生长培养基并找到最少的培养基成分。

<<<代码块_1>>>

## 工作流程 3：通过采样进行通量空间探索

此工作流程演示了使用 FVA 和采样进行全面的通量空间分析。

<<<代码块_2>>>

## 工作流程 4：生产菌株设计

此工作流程演示了如何设计目标代谢物的生产菌株。

<<<代码块_3>>>

## 工作流程 5：模型验证和调试

该工作流程展示了验证和调试代谢模型的系统方法。

<<<代码块_4>>>

这些工作流程为常见 COBRApy 任务提供了全面的模板。根据具体研究问题和模型的需要调整它们。