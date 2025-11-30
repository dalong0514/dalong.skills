<!-- 此文件由机器翻译自 SKILL.md -->

---
名称： 科布拉比
描述：“基于约束的代谢模型 (COBRA)。FBA、FVA、基因敲除、通量采样、SBML 模型，用于系统生物学和代谢工程分析。”
---

# COBRApy - 基于约束的重建和分析

## 概述

COBRApy 是一个用于代谢模型基于约束的重建和分析 (COBRA) 的 Python 库，对于系统生物学研究至关重要。使用基因组规模的代谢模型，对细胞代谢进行计算模拟，进行代谢工程分析并预测表型行为。

## 核心能力

COBRApy 提供了分为几个关键领域的综合工具：

### 1.模型管理

从存储库或文件加载现有模型：
```python
from cobra.io import load_model

# Load bundled test models
model = load_model("textbook")  # E. coli core model
model = load_model("ecoli")     # Full E. coli model
model = load_model("salmonella")

# Load from files
from cobra.io import read_sbml_model, load_json_model, load_yaml_model
model = read_sbml_model("path/to/model.xml")
model = load_json_model("path/to/model.json")
model = load_yaml_model("path/to/model.yml")
```

以各种格式保存模型：
<<<代码块_1>>>

### 2. 模型结构和组件

访问和检查模型组件：
<<<代码块_2>>>

### 3. 通量平衡分析 (FBA)

执行标准 FBA 模拟：
<<<代码块_3>>>

简约的 FBA（最小化总通量）：
<<<代码块_4>>>

几何FBA（寻找中心解）：
<<<代码块_5>>>

### 4. 通量变异分析 (FVA)

确定所有反应的通量范围：
<<<代码块_6>>>

### 5. 基因和反应删除研究

进行淘汰分析：
```python
from cobra.flux_analysis import (
    single_gene_deletion,
    single_reaction_deletion,
    double_gene_deletion,
    double_reaction_deletion
)

# Single deletions
gene_results = single_gene_deletion(model)
reaction_results = single_reaction_deletion(model)

# Double deletions (uses multiprocessing)
double_gene_results = double_gene_deletion(
    model,
    processes=4  # Number of CPU cores
)

# Manual knockout using context manager
with model:
    model.genes.get_by_id("b0008").knock_out()
    solution = model.optimize()
    print(f"Growth after knockout: {solution.objective_value}")
# Model automatically reverts after context exit
```

### 6.生长培养基和最小培养基

管理生长介质：
```python
# View current medium
print(model.medium)

# Modify medium (must reassign entire dict)
medium = model.medium
medium["EX_glc__D_e"] = 10.0  # Set glucose uptake
medium["EX_o2_e"] = 0.0       # Anaerobic conditions
model.medium = medium

# Calculate minimal media
from cobra.medium import minimal_medium

# Minimize total import flux
min_medium = minimal_medium(model, minimize_components=False)

# Minimize number of components (uses MILP, slower)
min_medium = minimal_medium(
    model,
    minimize_components=True,
    open_exchanges=True
)
```

### 7. 通量采样

对可行通量空间进行采样：
```python
from cobra.sampling import sample

# Sample using OptGP (default, supports parallel processing)
samples = sample(model, n=1000, method="optgp", processes=4)

# Sample using ACHR
samples = sample(model, n=1000, method="achr")

# Validate samples
from cobra.sampling import OptGPSampler
sampler = OptGPSampler(model, processes=4)
sampler.sample(1000)
validation = sampler.validate(sampler.samples)
print(validation.value_counts())  # Should be all 'v' for valid
```

### 8. 生产信封

计算表型相平面：
```python
from cobra.flux_analysis import production_envelope

# Standard production envelope
envelope = production_envelope(
    model,
    reactions=["EX_glc__D_e", "EX_o2_e"],
    objective="EX_ac_e"  # Acetate production
)

# With carbon yield
envelope = production_envelope(
    model,
    reactions=["EX_glc__D_e", "EX_o2_e"],
    carbon_sources="EX_glc__D_e"
)

# Visualize (use matplotlib or pandas plotting)
import matplotlib.pyplot as plt
envelope.plot(x="EX_glc__D_e", y="EX_o2_e", kind="scatter")
plt.show()
```

### 9. 填补空白

添加反应以使模型可行：
```python
from cobra.flux_analysis import gapfill

# Prepare universal model with candidate reactions
universal = load_model("universal")

# Perform gapfilling
with model:
    # Remove reactions to create gaps for demonstration
    model.remove_reactions([model.reactions.PGI])

    # Find reactions needed
    solution = gapfill(model, universal)
    print(f"Reactions to add: {solution}")
```

### 10. 模型构建

从头开始构建模型：
```python
from cobra import Model, Reaction, Metabolite

# Create model
model = Model("my_model")

# Create metabolites
atp_c = Metabolite("atp_c", formula="C10H12N5O13P3",
                   name="ATP", compartment="c")
adp_c = Metabolite("adp_c", formula="C10H12N5O10P2",
                   name="ADP", compartment="c")
pi_c = Metabolite("pi_c", formula="HO4P",
                  name="Phosphate", compartment="c")

# Create reaction
reaction = Reaction("ATPASE")
reaction.name = "ATP hydrolysis"
reaction.subsystem = "Energy"
reaction.lower_bound = 0.0
reaction.upper_bound = 1000.0

# Add metabolites with stoichiometry
reaction.add_metabolites({
    atp_c: -1.0,
    adp_c: 1.0,
    pi_c: 1.0
})

# Add gene-reaction rule
reaction.gene_reaction_rule = "(gene1 and gene2) or gene3"

# Add to model
model.add_reactions([reaction])

# Add boundary reactions
model.add_boundary(atp_c, type="exchange")
model.add_boundary(adp_c, type="demand")

# Set objective
model.objective = "ATPASE"
```

## 常见工作流程

### 工作流程 1：加载模型并预测增长

```python
from cobra.io import load_model

# Load model
model = load_model("ecoli")

# Run FBA
solution = model.optimize()
print(f"Growth rate: {solution.objective_value:.3f} /h")

# Show active pathways
print(solution.fluxes[solution.fluxes.abs() > 1e-6])
```

### 工作流程 2：基因敲除筛选

```python
from cobra.io import load_model
from cobra.flux_analysis import single_gene_deletion

# Load model
model = load_model("ecoli")

# Perform single gene deletions
results = single_gene_deletion(model)

# Find essential genes (growth < threshold)
essential_genes = results[results["growth"] < 0.01]
print(f"Found {len(essential_genes)} essential genes")

# Find genes with minimal impact
neutral_genes = results[results["growth"] > 0.9 * solution.objective_value]
```

### 工作流程 3：媒体优化

```python
from cobra.io import load_model
from cobra.medium import minimal_medium

# Load model
model = load_model("ecoli")

# Calculate minimal medium for 50% of max growth
target_growth = model.slim_optimize() * 0.5
min_medium = minimal_medium(
    model,
    target_growth,
    minimize_components=True
)

print(f"Minimal medium components: {len(min_medium)}")
print(min_medium)
```

### 工作流程 4：通量不确定性分析

```python
from cobra.io import load_model
from cobra.flux_analysis import flux_variability_analysis
from cobra.sampling import sample

# Load model
model = load_model("ecoli")

# First check flux ranges at optimality
fva = flux_variability_analysis(model, fraction_of_optimum=1.0)

# For reactions with large ranges, sample to understand distribution
samples = sample(model, n=1000)

# Analyze specific reaction
reaction_id = "PFK"
import matplotlib.pyplot as plt
samples[reaction_id].hist(bins=50)
plt.xlabel(f"Flux through {reaction_id}")
plt.ylabel("Frequency")
plt.show()
```

### 工作流程 5：临时更改的上下文管理器

使用上下文管理器进行临时修改：
```python
# Model remains unchanged outside context
with model:
    # Temporarily change objective
    model.objective = "ATPM"

    # Temporarily modify bounds
    model.reactions.EX_glc__D_e.lower_bound = -5.0

    # Temporarily knock out genes
    model.genes.b0008.knock_out()

    # Optimize with changes
    solution = model.optimize()
    print(f"Modified growth: {solution.objective_value}")

# All changes automatically reverted
solution = model.optimize()
print(f"Original growth: {solution.objective_value}")
```

## 关键概念

### DictList 对象
模型使用 `DictList` 对象来表示反应、代谢物和基因 - 行为类似于列表和字典：
```python
# Access by index
first_reaction = model.reactions[0]

# Access by ID
pfk = model.reactions.get_by_id("PFK")

# Query methods
atp_reactions = model.reactions.query("atp")
```

### 通量约束
反应界限定义了可行的通量范围：
- **不可逆**：`lower_bound = 0, upper_bound > 0`
- **可逆**：`lower_bound < 0, upper_bound > 0`
- 使用 `.bounds` 同时设置两个边界以避免不一致

### 基因反应规则 (GPR)
将基因与反应联系起来的布尔逻辑：
```python
# AND logic (both required)
reaction.gene_reaction_rule = "gene1 and gene2"

# OR logic (either sufficient)
reaction.gene_reaction_rule = "gene1 or gene2"

# Complex logic
reaction.gene_reaction_rule = "(gene1 and gene2) or (gene3 and gene4)"
```

### 交流反应
代表代谢物导入/导出的特殊反应：
- 按照惯例以前缀 `EX_` 命名
- 正通量 = 分泌，负通量 = 摄取
- 通过`model.medium`字典进行管理

## 最佳实践

1. **使用上下文管理器**进行临时修改以避免状态管理问题
2. **在分析之前使用`model.slim_optimize()`验证模型**以确保可行性
3.优化后**检查求解状态** - `optimal`表示求解成功
4. **当热力学可行性很重要时，使用无环 FVA**
5. **在FVA中适当设置fraction_of_optimum**以探索次优空间
6. **并行化**计算量大的操作（采样、双删除）
7. **首选 SBML 格式**用于模型交换和长期存储
8. **当性能只需要目标值时使用 slim_optimize()**
9. **验证通量样本**以确保数值稳定性

## 故障排除

**不可行的解决方案**：检查介质约束、反应界限和模型一致性
**缓慢优化**：通过 `model.solver` 尝试不同的求解器（GLPK、CPLEX、Gurobi）
**无界解决方案**：验证交换反应具有适当的上限
**导入错误**：确保正确的文件格式和有效的 SBML 标识符

## 参考文献

有关详细的工作流程和 API 模式，请参阅：
- `references/workflows.md` - 全面的分步工作流程示例
- `references/api_quick_reference.md` - 常见函数签名和模式

官方文档：https://cobrapy.readthedocs.io/en/latest/