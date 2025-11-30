<!-- 此文件由机器翻译自 constraints_mcdm.md -->

# Pymoo 约束和决策参考

pymoo 中约束处理和多标准决策的参考。

## 约束处理

### 定义约束

问题定义中指定了约束：

```python
from pymoo.core.problem import ElementwiseProblem
import numpy as np

class ConstrainedProblem(ElementwiseProblem):
    def __init__(self):
        super().__init__(
            n_var=2,
            n_obj=2,
            n_ieq_constr=2,    # Number of inequality constraints
            n_eq_constr=1,      # Number of equality constraints
            xl=np.array([0, 0]),
            xu=np.array([5, 5])
        )

    def _evaluate(self, x, out, *args, **kwargs):
        # Objectives
        f1 = x[0]**2 + x[1]**2
        f2 = (x[0]-1)**2 + (x[1]-1)**2

        out["F"] = [f1, f2]

        # Inequality constraints (formulated as g(x) <= 0)
        g1 = x[0] + x[1] - 5  # x[0] + x[1] >= 5 → -(x[0] + x[1] - 5) <= 0
        g2 = x[0]**2 + x[1]**2 - 25  # x[0]^2 + x[1]^2 <= 25

        out["G"] = [g1, g2]

        # Equality constraints (formulated as h(x) = 0)
        h1 = x[0] - 2*x[1]

        out["H"] = [h1]
```

**约束制定规则：**
- 不等式：`g(x) <= 0`（负数或零时可行）
- 相等：`h(x) = 0`（零时可行）
- 将 `g(x) >= 0` 转换为 `-g(x) <= 0`

### 约束处理技术

#### 1. 可行性第一（默认）
**机制：** 总是优先选择可行的解决方案而不是不可行的解决方案
**比较：**
1. 都可行 → 按客观值比较
2. 一种可行，一种不可行→可行胜
3. 两者都不可行 → 通过违反约束进行比较

**用途：**
<<<代码块_1>>>

**优点：**
- 适用于任何基于排序的算法
- 简单有效
- 无需参数调整

**缺点：**
- 可能会与小的可行区域发生冲突
- 可以忽略好的不可行的解决方案

#### 二、处罚方法
**机制：** 根据违反约束条件对目标添加惩罚
**公式：** `F_penalized = F + penalty_factor * violation`

**用途：**
<<<代码块_2>>>

**参数：**
- `penalty`：惩罚系数（根据问题规模调整）

**优点：**
- 将受约束问题转换为无约束问题
- 适用于任何优化算法

**缺点：**
- 惩罚参数敏感
- 可能需要针对具体问题进行调整

#### 3. 约束作为目标
**机制：** 将违反约束视为附加目标
**结果：** M+1个目标的多目标问题（M个原始+约束）

**用途：**
<<<代码块_3>>>

**优点：**
- 无需参数调整
- 维护可能有用的不可行的解决方案
- 当可行区域较小时效果很好

**缺点：**
- 增加问题维度
- 更复杂的帕累托前沿分析

#### 4. Epsilon 约束处理
**机制：** 动态可行性阈值
**概念：** 几代人逐渐收紧约束容差

**优点：**
- 平滑过渡到可行区域
- 帮助解决困难的约束景观

**缺点：**
- 特定算法的实现
- 需要参数调整

#### 5. 维修操作员
**机制：**修改不可行的解来满足约束
**应用：** 交叉/变异后，修复后代

**用途：**
<<<代码块_4>>>

**优点：**
- 在整个优化过程中保持可行性
- 可以对领域知识进行编码

**缺点：**
- 需要针对具体问题实施
- 可能会限制搜索

### 约束处理算法

一些算法具有内置的约束处理：

#### SRES（随机排名进化策略）
**目的：** 单目标约束优化
**机制：** 随机排名平衡目标和约束

**用途：**
<<<代码块_5>>>

#### ISRES（改进的 SRES）
**目的：**增强约束优化
**改进：**更好的参数适配

**用途：**
<<<代码块_6>>>

### 约束处理指南

**选择技术基于：**

|问题特征|推荐技术 |
|------------------------------------|------------------------|
|大可行域|可行性第一 |
|小可行域|约束作为目标，修复|
|严重受限| SRES/ISRES、Epsilon 约束 |
|线性约束|修复（投影）|
|非线性约束|可行性第一，惩罚|
|已知可行的解决方案|有偏初始化 |

## 多标准决策 (MCDM)

获得帕累托前沿后，MCDM 帮助选择首选解决方案。

### 决策背景

**帕累托前沿特征：**
- 多个非支配解
- 每个代表不同的权衡
- 没有客观的“最佳”解决方案
- 需要决策者的偏好

### Pymoo 中的 MCDM 方法

#### 1. 伪权重
**概念：** 对每个目标进行加权，选择最小化加权和的解决方案
**公式：** `score = w1*f1 + w2*f2 + ... + wM*fM`

**用途：**
```python
from pymoo.mcdm.pseudo_weights import PseudoWeights

# Define weights (must sum to 1)
weights = np.array([0.3, 0.7])  # 30% weight on f1, 70% on f2

dm = PseudoWeights(weights)
best_idx = dm.do(result.F)
best_solution = result.X[best_idx]
```

**何时使用：**
- 提供明确的偏好表达
- 目标相当
- 可接受的线性权衡

**限制：**
- 需要重量规格
- 线性假设可能无法捕捉偏好
- 对目标缩放敏感

#### 2. 妥协编程
**概念：**选择最接近理想点的解
**度量：** 到理想的距离（例如，欧几里得、切比雪夫）

**用途：**
```python
from pymoo.mcdm.compromise_programming import CompromiseProgramming

dm = CompromiseProgramming()
best_idx = dm.do(result.F, ideal=ideal_point, nadir=nadir_point)
```

**何时使用：**
- 已知或估计的理想目标值
- 平衡考虑所有目标
- 没有明确的体重偏好

#### 3. 交互式决策
**概念：** 迭代偏好细化
**流程：**
1.向决策者展示有代表性的解决方案
2. 收集偏好反馈
3. 集中搜索首选区域
4. 重复直到找到满意的解决方案

**方法：**
- 参考点方法
- 权衡分析
- 渐进式偏好表达

### 决策工作流程

**第 1 步：标准化目标**
```python
# Normalize to [0, 1] for fair comparison
F_norm = (result.F - result.F.min(axis=0)) / (result.F.max(axis=0) - result.F.min(axis=0))
```

**第 2 步：分析权衡**
```python
from pymoo.visualization.scatter import Scatter

plot = Scatter()
plot.add(result.F)
plot.show()

# Identify knee points, extreme solutions
```

**步骤 3：应用 MCDM 方法**
```python
from pymoo.mcdm.pseudo_weights import PseudoWeights

weights = np.array([0.4, 0.6])  # Based on preferences
dm = PseudoWeights(weights)
selected = dm.do(F_norm)
```

**第 4 步：验证选择**
```python
# Visualize selected solution
from pymoo.visualization.petal import Petal

plot = Petal()
plot.add(result.F[selected], label="Selected")
# Add other candidates for comparison
plot.show()
```

### 高级 MCDM 技术

#### 拐点检测
**概念：** 一个目标的微小改进会导致其他目标大幅下降的解决方案

**用途：**
```python
from pymoo.mcdm.knee import KneePoint

km = KneePoint()
knee_idx = km.do(result.F)
knee_solutions = result.X[knee_idx]
```

**何时使用：**
- 没有明确的偏好
- 需要平衡的权衡
- 凸帕累托前沿

#### 超容量贡献
**概念：** 选择对超容量贡献最大的解决方案
**用例：** 维护解决方案的不同子集

**用途：**
```python
from pymoo.indicators.hv import HV

hv = HV(ref_point=reference_point)
hv_contributions = hv.calc_contributions(result.F)

# Select top contributors
top_k = 5
top_indices = np.argsort(hv_contributions)[-top_k:]
selected_solutions = result.X[top_indices]
```

### 决策指南

**当决策者有：**

|偏好信息 |推荐方法 |
|------------------------------------|--------------------|
|明确的目标权重 |伪权重 |
|理想目标值|妥协编程|
|没有优先偏好 |拐点，目视检查|
|相互矛盾的标准 |互动方式|
|需要不同的子集 |超容量贡献 |

**最佳实践：**
1. **在 MCDM 之前将目标标准化**
2. **可视化帕累托前沿**以了解权衡
3. **考虑多种方法**进行稳健选择
4. **与领域专家验证结果**
5. **记录假设**和偏好来源
6. **对权重/参数进行敏感性分析**

### 集成示例

包含约束处理和决策的完整工作流程：

```python
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.mcdm.pseudo_weights import PseudoWeights
import numpy as np

# Define constrained problem
problem = MyConstrainedProblem()

# Setup algorithm with feasibility-first constraint handling
algorithm = NSGA2(
    pop_size=100,
    eliminate_duplicates=True
)

# Optimize
result = minimize(
    problem,
    algorithm,
    ('n_gen', 200),
    seed=1,
    verbose=True
)

# Filter feasible solutions only
feasible_mask = result.CV[:, 0] == 0  # Constraint violation = 0
F_feasible = result.F[feasible_mask]
X_feasible = result.X[feasible_mask]

# Normalize objectives
F_norm = (F_feasible - F_feasible.min(axis=0)) / (F_feasible.max(axis=0) - F_feasible.min(axis=0))

# Apply MCDM
weights = np.array([0.5, 0.5])
dm = PseudoWeights(weights)
best_idx = dm.do(F_norm)

# Get final solution
best_solution = X_feasible[best_idx]
best_objectives = F_feasible[best_idx]

print(f"Selected solution: {best_solution}")
print(f"Objective values: {best_objectives}")
```