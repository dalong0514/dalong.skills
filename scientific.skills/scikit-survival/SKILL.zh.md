<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：scikit-生存
描述：使用 scikit-survival 在 Python 中进行生存分析和事件时间建模的综合工具包。在处理经过审查的生存数据、执行事件时间分析、拟合 Cox 模型、随机生存森林、梯度提升模型或生存 SVM、使用一致性指数或 Brier 评分评估生存预测、处理竞争风险或使用 scikit-survival 库实施任何生存分析工作流程时，可以使用此技能。
---

# scikit-survival：Python 中的生存分析

## 概述

scikit-survival 是一个基于 scikit-learn 构建的用于生存分析的 Python 库。它提供了用于事件发生时间分析的专用工具，可应对某些观察结果仅部分已知的审查数据的独特挑战。

生存分析旨在建立协变量和事件发生时间之间的联系，并考虑审查记录（特别是参与者在观察期间没有经历事件的研究中的右审查数据）。

## 何时使用此技能

在以下情况下使用此技能：
- 执行生存分析或事件时间建模
- 使用审查数据（右审查、左审查或区间审查）
- 拟合Cox比例风险模型（标准或惩罚）
- 建立整体生存模型（随机生存森林、梯度提升）
- 训练生存支持向量机
- 评估生存模型性能（一致性指数、Brier 评分、时间依赖性 AUC）
- 估计 Kaplan-Meier 或 Nelson-Aalen 曲线
- 分析竞争风险
- 预处理生存数据或处理生存数据集中的缺失值
- 使用 scikit-survival 库进行任何分析

## 核心能力

### 1. 模型类型及选择

scikit-survival 提供了多个模型系列，每个模型系列适用于不同的场景：

#### Cox 比例风险模型
**用途**：具有可解释系数的标准生存分析
- `CoxPHSurvivalAnalysis`：基本 Cox 模型
- `CoxnetSurvivalAnalysis`：针对高维数据使用弹性网络惩罚 Cox
- `IPCRidge`：加速失效时间模型的岭回归

**请参阅**：`references/cox-models.md` 有关 Cox 模型、正则化和解释的详细指导

#### 集成方法
**用途**：具有复杂非线性关系的高预测性能
- `RandomSurvivalForest`：稳健的非参数集成方法
- `GradientBoostingSurvivalAnalysis`：基于树的提升以获得最大性能
- `ComponentwiseGradientBoostingSurvivalAnalysis`：具有特征选择的线性增强
- `ExtraSurvivalTrees`：极其随机的树，用于额外的正则化

**请参阅**：`references/ensemble-models.md` 有关集成方法、超参数调整以及何时使用每个模型的综合指导

#### 生存支持向量机
**用途**：具有基于边际学习的中型数据集
- `FastSurvivalSVM`：针对速度进行优化的线性 SVM
- `FastKernelSurvivalSVM`：非线性关系的内核 SVM
- `HingeLossSurvivalSVM`：具有铰链损失的 SVM
- `ClinicalKernelTransform`：临床+分子数据的专用内核

**请参阅**：`references/svm-models.md` 了解详细的 SVM 指南、内核选择和超参数调整

#### 模型选择决策树

```
Start
├─ High-dimensional data (p > n)?
│  ├─ Yes → CoxnetSurvivalAnalysis (elastic net)
│  └─ No → Continue
│
├─ Need interpretable coefficients?
│  ├─ Yes → CoxPHSurvivalAnalysis or ComponentwiseGradientBoostingSurvivalAnalysis
│  └─ No → Continue
│
├─ Complex non-linear relationships expected?
│  ├─ Yes
│  │  ├─ Large dataset (n > 1000) → GradientBoostingSurvivalAnalysis
│  │  ├─ Medium dataset → RandomSurvivalForest or FastKernelSurvivalSVM
│  │  └─ Small dataset → RandomSurvivalForest
│  └─ No → CoxPHSurvivalAnalysis or FastSurvivalSVM
│
└─ For maximum performance → Try multiple models and compare
```

### 2. 数据准备和预处理

在建模之前，正确准备生存数据：

#### 创造生存成果
<<<代码块_1>>>

#### 基本预处理步骤
1. **处理缺失值**：特征的插补策略
2. **对分类变量进行编码**：One-hot编码或标签编码
3. **标准化特征**：对于 SVM 和正则化 Cox 模型至关重要
4. **验证数据质量**：检查负时间、每个功能是否有足够的事件
5. **训练-测试分割**：在分割之间保持相似的审查率

**请参阅**：`references/data-handling.md` 了解完整的预处理工作流程、数据验证和最佳实践

### 3.模型评估

正确的评估对于生存模型至关重要。使用适当的指标来进行审查：

#### 一致性指数（C 指数）
排名/歧视的主要指标：
- **Harrell 的 C 指数**：用于低审查 (<40%)
- **Uno 的 C 指数**：用于中度至高度审查 (>40%) - 更稳健

<<<代码块_2>>>

#### 时间依赖性 AUC
评估特定时间点的歧视：

<<<代码块_3>>>

#### 荆棘分数
评估歧视和校准：

<<<代码块_4>>>
**请参阅**：`references/evaluation-metrics.md` 了解综合评估指南、指标选择以及使用具有交叉验证的评分器

### 4.竞争风险分析

处理具有多个互斥事件类型的情况：

<<<代码块_5>>>

**在以下情况下使用竞争风险：
- 存在多种相互排斥的事件类型（例如，不同原因导致的死亡）
- 一个事件的发生会妨碍其他事件的发生
- 需要对特定事件类型进行概率估计

**参见**：`references/competing-risks.md` 了解详细的竞争风险方法、特定原因的危害模型和解释

### 5.非参数估计

在没有参数假设的情况下估计生存函数：

#### Kaplan-Meier 估计器
<<<代码块_6>>>

#### Nelson-Aalen 估计器
```python
from sksurv.nonparametric import nelson_aalen_estimator

time, cumulative_hazard = nelson_aalen_estimator(y['event'], y['time'])
```

## 典型工作流程

### 工作流程 1：标准生存分析

```python
from sksurv.datasets import load_breast_cancer
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.metrics import concordance_index_ipcw
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

# 1. Load and prepare data
X, y = load_breast_cancer()
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 2. Preprocess
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# 3. Fit model
estimator = CoxPHSurvivalAnalysis()
estimator.fit(X_train_scaled, y_train)

# 4. Predict
risk_scores = estimator.predict(X_test_scaled)

# 5. Evaluate
c_index = concordance_index_ipcw(y_train, y_test, risk_scores)[0]
print(f"C-index: {c_index:.3f}")
```

### 工作流程 2：具有特征选择的高维数据

```python
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sklearn.model_selection import GridSearchCV
from sksurv.metrics import as_concordance_index_ipcw_scorer

# 1. Use penalized Cox for feature selection
estimator = CoxnetSurvivalAnalysis(l1_ratio=0.9)  # Lasso-like

# 2. Tune regularization with cross-validation
param_grid = {'alpha_min_ratio': [0.01, 0.001]}
cv = GridSearchCV(estimator, param_grid,
                  scoring=as_concordance_index_ipcw_scorer(), cv=5)
cv.fit(X, y)

# 3. Identify selected features
best_model = cv.best_estimator_
selected_features = np.where(best_model.coef_ != 0)[0]
```

### 工作流程 3：实现最大性能的集成方法

```python
from sksurv.ensemble import GradientBoostingSurvivalAnalysis
from sklearn.model_selection import GridSearchCV

# 1. Define parameter grid
param_grid = {
    'learning_rate': [0.01, 0.05, 0.1],
    'n_estimators': [100, 200, 300],
    'max_depth': [3, 5, 7]
}

# 2. Grid search
gbs = GradientBoostingSurvivalAnalysis()
cv = GridSearchCV(gbs, param_grid, cv=5,
                  scoring=as_concordance_index_ipcw_scorer(), n_jobs=-1)
cv.fit(X_train, y_train)

# 3. Evaluate best model
best_model = cv.best_estimator_
risk_scores = best_model.predict(X_test)
c_index = concordance_index_ipcw(y_train, y_test, risk_scores)[0]
```

### 工作流程 4：综合模型比较

```python
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.ensemble import RandomSurvivalForest, GradientBoostingSurvivalAnalysis
from sksurv.svm import FastSurvivalSVM
from sksurv.metrics import concordance_index_ipcw, integrated_brier_score

# Define models
models = {
    'Cox': CoxPHSurvivalAnalysis(),
    'RSF': RandomSurvivalForest(n_estimators=100, random_state=42),
    'GBS': GradientBoostingSurvivalAnalysis(random_state=42),
    'SVM': FastSurvivalSVM(random_state=42)
}

# Evaluate each model
results = {}
for name, model in models.items():
    model.fit(X_train_scaled, y_train)
    risk_scores = model.predict(X_test_scaled)
    c_index = concordance_index_ipcw(y_train, y_test, risk_scores)[0]
    results[name] = c_index
    print(f"{name}: C-index = {c_index:.3f}")

# Select best model
best_model_name = max(results, key=results.get)
print(f"\nBest model: {best_model_name}")
```

## 与 scikit-learn 集成

scikit-survival 与 scikit-learn 的生态系统完全集成：

```python
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, GridSearchCV

# Use pipelines
pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('model', CoxPHSurvivalAnalysis())
])

# Use cross-validation
scores = cross_val_score(pipeline, X, y, cv=5,
                         scoring=as_concordance_index_ipcw_scorer())

# Use grid search
param_grid = {'model__alpha': [0.1, 1.0, 10.0]}
cv = GridSearchCV(pipeline, param_grid, cv=5)
cv.fit(X, y)
```

## 最佳实践

1. **始终标准化 SVM 和正则化 Cox 模型的特征**
2. **当审查 > 40% 时，使用 Uno 的 C 指数** 而不是 Harrell 的
3. **报告多个评估指标**（C 指数、综合 Brier 评分、时间相关 AUC）
4. **检查 Cox 模型的比例风险假设**
5. **使用交叉验证**通过适当的评分器进行超参数调整
6. **在建模之前验证数据质量**（检查负时间、每个特征是否有足够的事件）
7. **比较多种模型类型**以找到最佳性能
8. **对随机生存森林使用排列重要性**（不是内置重要性）
9. **当存在多种事件类型时，考虑竞争风险**
10. **文件审查机制**和分析率

## 要避免的常见陷阱

1. **使用 Harrell 的 C 指数和高审查** → 使用 Uno 的 C 指数
2. **不标准化 SVM 的功能** → 始终标准化
3. **忘记将 y_train 传递给 concordance_index_ipcw** → IPCW 计算所需
4. **将竞争事件视为审查** → 使用竞争风险方法
5. **没有检查每个功能是否有足够的事件** → 经验法则：每个功能 10 个以上事件
6. **使用 RSF 的内置特征重要性** → 使用排列重要性
7. **忽略比例风险假设** → 验证或使用替代模型
8. **在交叉验证中未使用适当的评分器** → 使用 as_concordance_index_ipcw_scorer()

## 参考文件

该技能包括特定主题的详细参考文件：

- **`references/cox-models.md`**：Cox 比例风险模型、惩罚 Cox (CoxNet)、IPRidge、正则化策略和解释的完整指南
- **`references/ensemble-models.md`**：随机生存森林、梯度提升、超参数调整、特征重要性和模型选择
- **`references/evaluation-metrics.md`**：一致性指数（Harrell's 与 Uno's）、时间相关 AUC、Brier 分数、综合评估流程
- **`references/data-handling.md`**：数据加载、预处理工作流程、处理缺失数据、特征编码、验证检查
- **`references/svm-models.md`**：生存支持向量机、内核选择、临床内核变换、超参数调整
- **`references/competing-risks.md`**：竞争风险分析、累积发生函数、特定原因危害模型

当特定任务需要详细信息时，加载这些参考文件。

## 其他资源

- **官方文档**：https://scikit-survival.readthedocs.io/
- **GitHub 存储库**：https://github.com/sebp/scikit-survival
- **内置数据集**：使用 `sksurv.datasets` 作为练习数据集（GBSG2、WHAS500、退伍军人肺癌等）
- **API 参考**：类和函数的完整列表位于 https://scikit-survival.readthedocs.io/en/stable/api/index.html

## 快速参考：关键导入

```python
# Models
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis, IPCRidge
from sksurv.ensemble import RandomSurvivalForest, GradientBoostingSurvivalAnalysis
from sksurv.svm import FastSurvivalSVM, FastKernelSurvivalSVM
from sksurv.tree import SurvivalTree

# Evaluation metrics
from sksurv.metrics import (
    concordance_index_censored,
    concordance_index_ipcw,
    cumulative_dynamic_auc,
    brier_score,
    integrated_brier_score,
    as_concordance_index_ipcw_scorer,
    as_integrated_brier_score_scorer
)

# Non-parametric estimation
from sksurv.nonparametric import (
    kaplan_meier_estimator,
    nelson_aalen_estimator,
    cumulative_incidence_competing_risks
)

# Data handling
from sksurv.util import Surv
from sksurv.preprocessing import OneHotEncoder, encode_categorical
from sksurv.datasets import load_gbsg2, load_breast_cancer, load_veterans_lung_cancer

# Kernels
from sksurv.kernels import ClinicalKernelTransform
```