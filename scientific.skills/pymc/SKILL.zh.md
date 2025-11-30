<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：pymc-贝叶斯建模
描述：“使用 PyMC 进行贝叶斯建模。构建分层模型、MCMC (NUTS)、变分推理、LOO/WAIC 比较、后验检查，用于概率编程和推理。”
---

# PyMC 贝叶斯建模

## 概述

PyMC 是一个用于贝叶斯建模和概率编程的 Python 库。使用 PyMC 的现代 API（版本 5.x+）构建、拟合、验证和比较贝叶斯模型，包括分层模型、MCMC 采样 (NUTS)、变分推理和模型比较（LOO、WAIC）。

## 何时使用此技能

该技能应该在以下情况下使用：
- 建立贝叶斯模型（线性/逻辑回归、层次模型、时间序列等）
- 执行 MCMC 采样或变分推理
- 进行事前/事后预测检查
- 诊断采样问题（分歧、收敛、ESS）
- 使用信息标准（LOO、WAIC）比较多个模型
- 通过贝叶斯方法实现不确定性量化
- 使用分层/多级数据结构
- 以原则性的方式处理缺失数据或测量错误

## 标准贝叶斯工作流程

按照以下工作流程构建和验证贝叶斯模型：

### 1. 数据准备

```python
import pymc as pm
import arviz as az
import numpy as np

# Load and prepare data
X = ...  # Predictors
y = ...  # Outcomes

# Standardize predictors for better sampling
X_mean = X.mean(axis=0)
X_std = X.std(axis=0)
X_scaled = (X - X_mean) / X_std
```

**关键做法：**
- 标准化连续预测器（提高采样效率）
- 尽可能以结果为中心
- 显式处理缺失数据（视为参数）
- 为了清楚起见，使用带有 `coords` 的命名维度

### 2. 模型构建

<<<代码块_1>>>

**关键做法：**
- 使用弱信息先验（不是平坦先验）
- 使用 `HalfNormal` 或 `Exponential` 作为比例参数
- 尽可能使用命名维度 (`dims`) 而不是 `shape`
- 使用 `pm.Data()` 作为将更新预测的值

### 3. 事前预测检查

**在拟合之前始终验证先验：**

<<<代码块_2>>>

**检查：**
- 先前的预测是否涵盖合理的值？
- 考虑到领域知识，极值是否合理？
- 如果先验生成不可信的数据，请调整并重新检查

### 4. 拟合模型

<<<代码块_3>>>

**关键参数：**
- `draws=2000`：每个链的样本数
- `tune=1000`：预热样本（丢弃）
- `chains=4`：运行 4 个链进行收敛检查
- `target_accept=0.9`：对于困难的后验来说更高（0.95-0.99）
- 包括 `log_likelihood=True` 用于模型比较

### 5. 检查诊断

**使用诊断脚本：**

<<<代码块_4>>>

**检查：**
- **R-hat < 1.01**：链已收敛
- **ESS > 400**：足够的有效样本
- **无分歧**：NUTS 采样成功
- **跟踪图**：链条应该混合良好（模糊的毛毛虫）

**如果出现问题：**
- 分歧→增加`target_accept=0.95`，使用非中心参数化
- 低 ESS → 采样更多抽奖，重新参数化以减少相关性
- 高 R 帽 → 跑得更远，检查多模态

### 6. 事后预测检查

**验证模型拟合：**

<<<代码块_5>>>

**检查：**
- 后验预测是否捕获了观察到的数据模式？
- 系统偏差是否明显（模型指定错误）？
- 如果合身性较差，请考虑替代型号

### 7. 分析结果

<<<代码块_6>>>

### 8. 做出预测

```python
X_new = ...  # New predictor values
X_new_scaled = (X_new - X_mean) / X_std

with model:
    pm.set_data({'X_scaled': X_new_scaled})
    post_pred = pm.sample_posterior_predictive(
        idata.posterior,
        var_names=['y_obs'],
        random_seed=42
    )

# Extract prediction intervals
y_pred_mean = post_pred.posterior_predictive['y_obs'].mean(dim=['chain', 'draw'])
y_pred_hdi = az.hdi(post_pred.posterior_predictive, var_names=['y_obs'])
```

## 常见模型模式

### 线性回归

对于具有线性关系的连续结果：

```python
with pm.Model() as linear_model:
    alpha = pm.Normal('alpha', mu=0, sigma=10)
    beta = pm.Normal('beta', mu=0, sigma=10, shape=n_predictors)
    sigma = pm.HalfNormal('sigma', sigma=1)

    mu = alpha + pm.math.dot(X, beta)
    y = pm.Normal('y', mu=mu, sigma=sigma, observed=y_obs)
```

**使用模板：** `assets/linear_regression_template.py`

### 逻辑回归

对于二元结果：

```python
with pm.Model() as logistic_model:
    alpha = pm.Normal('alpha', mu=0, sigma=10)
    beta = pm.Normal('beta', mu=0, sigma=10, shape=n_predictors)

    logit_p = alpha + pm.math.dot(X, beta)
    y = pm.Bernoulli('y', logit_p=logit_p, observed=y_obs)
```

### 分层模型

对于分组数据（使用非中心参数化）：

```python
with pm.Model(coords={'groups': group_names}) as hierarchical_model:
    # Hyperpriors
    mu_alpha = pm.Normal('mu_alpha', mu=0, sigma=10)
    sigma_alpha = pm.HalfNormal('sigma_alpha', sigma=1)

    # Group-level (non-centered)
    alpha_offset = pm.Normal('alpha_offset', mu=0, sigma=1, dims='groups')
    alpha = pm.Deterministic('alpha', mu_alpha + sigma_alpha * alpha_offset, dims='groups')

    # Observation-level
    mu = alpha[group_idx]
    sigma = pm.HalfNormal('sigma', sigma=1)
    y = pm.Normal('y', mu=mu, sigma=sigma, observed=y_obs)
```

**使用模板：** `assets/hierarchical_model_template.py`

**关键：** 始终对分层模型使用非中心参数化以避免发散。

### 泊松回归

对于计数数据：

```python
with pm.Model() as poisson_model:
    alpha = pm.Normal('alpha', mu=0, sigma=10)
    beta = pm.Normal('beta', mu=0, sigma=10, shape=n_predictors)

    log_lambda = alpha + pm.math.dot(X, beta)
    y = pm.Poisson('y', mu=pm.math.exp(log_lambda), observed=y_obs)
```

对于过度分散的计数，请改用 `NegativeBinomial`。

### 时间序列

对于自回归过程：

```python
with pm.Model() as ar_model:
    sigma = pm.HalfNormal('sigma', sigma=1)
    rho = pm.Normal('rho', mu=0, sigma=0.5, shape=ar_order)
    init_dist = pm.Normal.dist(mu=0, sigma=sigma)

    y = pm.AR('y', rho=rho, sigma=sigma, init_dist=init_dist, observed=y_obs)
```

## 型号对比

### 比较模型

使用LOO或WAIC进行模型比较：

```python
from scripts.model_comparison import compare_models, check_loo_reliability

# Fit models with log_likelihood
models = {
    'Model1': idata1,
    'Model2': idata2,
    'Model3': idata3
}

# Compare using LOO
comparison = compare_models(models, ic='loo')

# Check reliability
check_loo_reliability(models)
```

**释义：**
- **Δloo < 2**：模型相似，选择更简单的模型
- **2 < Δloo < 4**：更好模型的证据较弱
- **4 < Δloo < 10**：中等证据
- **Δloo > 10**：更好模型的有力证据

**检查 Pareto-k 值：**
- k < 0.7：LOO 可靠
- k > 0.7：考虑 WAIC 或 k 倍 CV

### 模型平均

当模型相似时，平均预测：

```python
from scripts.model_comparison import model_averaging

averaged_pred, weights = model_averaging(models, var_name='y_obs')
```
## 发行版选择指南

### 对于先修者

**尺度参数** (σ, τ):
- `pm.HalfNormal('sigma', sigma=1)` - 默认选择
- `pm.Exponential('sigma', lam=1)` - 替代方案
- `pm.Gamma('sigma', alpha=2, beta=1)` - 更多信息

**无界参数**：
- `pm.Normal('theta', mu=0, sigma=1)` - 对于标准化数据
- `pm.StudentT('theta', nu=3, mu=0, sigma=1)` - 对异常值具有鲁棒性

**正参数**：
- `pm.LogNormal('theta', mu=0, sigma=1)`
- `pm.Gamma('theta', alpha=2, beta=1)`

**概率**：
- `pm.Beta('p', alpha=2, beta=2)` - 信息量较弱
- `pm.Uniform('p', lower=0, upper=1)` - 非信息性（谨慎使用）

**相关矩阵**：
- `pm.LKJCorr('corr', n=n_vars, eta=2)` - eta=1 统一，eta>1 更喜欢同一性

### 对于可能性

**持续成果**：
- `pm.Normal('y', mu=mu, sigma=sigma)` - 连续数据的默认值
- `pm.StudentT('y', nu=nu, mu=mu, sigma=sigma)` - 对异常值具有鲁棒性

**计数数据**：
- `pm.Poisson('y', mu=lambda)` - 等分散计数
- `pm.NegativeBinomial('y', mu=mu, alpha=alpha)` - 计数过度分散
- `pm.ZeroInflatedPoisson('y', psi=psi, mu=mu)` - 多余的零

**二元结果**：
- `pm.Bernoulli('y', p=p)` 或 `pm.Bernoulli('y', logit_p=logit_p)`

**分类结果**：
- `pm.Categorical('y', p=probs)`

**请参阅：** `references/distributions.md` 获取全面的分发参考

## 采样和推理

### MCMC 与坚果

大多数型号的默认值和推荐值：

```python
idata = pm.sample(
    draws=2000,
    tune=1000,
    chains=4,
    target_accept=0.9,
    random_seed=42
)
```

**需要时调整：**
- 分歧 → `target_accept=0.95` 或更高
- 慢速采样→使用ADVI进行初始化
- 离散参数 → 对离散变量使用 `pm.Metropolis()`

### 变分推理

探索或初始化的快速近似：

```python
with model:
    approx = pm.fit(n=20000, method='advi')

    # Use for initialization
    start = approx.sample(return_inferencedata=False)[0]
    idata = pm.sample(start=start)
```

**权衡：**
- 比MCMC快得多
- 近似值（可能低估不确定性）
- 适合大型模型或快速探索

**请参阅：** `references/sampling_inference.md` 了解详细的采样指南

## 诊断脚本

### 综合诊断

```python
from scripts.model_diagnostics import create_diagnostic_report

create_diagnostic_report(
    idata,
    var_names=['alpha', 'beta', 'sigma'],
    output_dir='diagnostics/'
)
```

创建：
- 跟踪图
- 排名图（混合检查）
- 自相关图
- 能量图
- ESS的演变
- 汇总统计数据 CSV

### 快速诊断检查

```python
from scripts.model_diagnostics import check_diagnostics

results = check_diagnostics(idata)
```

检查 R-hat、ESS、分歧和树深度。

## 常见问题及解决方案

### 分歧

**症状：** `idata.sample_stats.diverging.sum() > 0`

**解决方案：**
1.增加`target_accept=0.95`或`0.99`
2.使用非中心参数化（分层模型）
3.添加更强的先验来约束参数
4. 检查型号规格是否错误

### 有效样本量低

**症状：** `ESS < 400`

**解决方案：**
1. 更多抽奖示例：`draws=5000`
2. 重新参数化以减少后验相关性
3. 使用 QR 分解进行相关预测变量的回归

### 高 R 帽

**症状：** `R-hat > 1.01`

**解决方案：**
1. 运行更长的链：`tune=2000, draws=5000`
2. 检查多模态
3. 使用 ADVI 改进初始化

### 慢速采样

**解决方案：**
1.使用ADVI初始化
2. 降低模型复杂度
3. 增加并行度：`cores=8, chains=8`
4. 如果合适的话使用变分推理

## 最佳实践

### 模型构建

1. **始终标准化预测变量**以获得更好的采样
2. **使用信息量较弱的先验**（不平坦）
3. **为了清晰起见，使用命名维度** (`dims`)
4. **分层模型的非中心参数化**
5. **在拟合之前检查先前的预测**

### 采样

1. **运行多个链**（至少 4 个）以实现收敛
2. **使用 `target_accept=0.9`** 作为基线（如果需要，可以更高）
3. **包含 `log_likelihood=True`** 用于模型比较
4. **设置随机种子**以实现可重复性

### 验证

1. **在解释之前检查诊断**（R-hat、ESS、分歧）
2. **后预测检查** 用于模型验证
3. **在适当的时候比较多个模型**
4. **报告不确定性**（HDI区间，不仅仅是点估计）

### 工作流程

1.从简单开始，逐渐增加复杂性
2. 事前预测检查→拟合→诊断→事后预测检查
3. 基于检查迭代模型规范
4. 记录假设和先前的选择

## 资源

该技能包括：

### 参考文献 (`references/`)

- **`distributions.md`**：按类别（连续、离散、多元、混合、时间序列）组织的 PyMC 分布的综合目录。在选择先验或可能性时使用。

- **`sampling_inference.md`**：采样算法（NUTS、Metropolis、SMC）、变分推理（ADVI、SVGD）和处理采样问题的详细指南。当遇到收敛问题或选择推理方法时使用。
- **`workflows.md`**：常见模型类型、数据准备、事先选择和模型验证的完整工作流程示例和代码模式。用作标准贝叶斯分析的食谱。

### 脚本 (`scripts/`)

- **`model_diagnostics.py`**：自动诊断检查和报告生成。功能：`check_diagnostics()`用于快速检查，`create_diagnostic_report()`用于通过绘图进行综合分析。

- **`model_comparison.py`**：使用 LOO/WAIC 的模型比较实用程序。函数：`compare_models()`、`check_loo_reliability()`、`model_averaging()`。

### 模板 (`assets/`)

- **`linear_regression_template.py`**：具有完整工作流程的贝叶斯线性回归的完整模板（数据准备、事先检查、拟合、诊断、预测）。

- **`hierarchical_model_template.py`**：具有非中心参数化和组级分析的分层/多级模型的完整模板。

## 快速参考

### 模型构建
```python
with pm.Model(coords={'var': names}) as model:
    # Priors
    param = pm.Normal('param', mu=0, sigma=1, dims='var')
    # Likelihood
    y = pm.Normal('y', mu=..., sigma=..., observed=data)
```

### 采样
```python
idata = pm.sample(draws=2000, tune=1000, chains=4, target_accept=0.9)
```

### 诊断
```python
from scripts.model_diagnostics import check_diagnostics
check_diagnostics(idata)
```

### 模型比较
```python
from scripts.model_comparison import compare_models
compare_models({'m1': idata1, 'm2': idata2}, ic='loo')
```

### 预测
```python
with model:
    pm.set_data({'X': X_new})
    pred = pm.sample_posterior_predictive(idata.posterior)
```

## 附加说明

- PyMC 与 ArviZ 集成以实现可视化和诊断
- 使用`pm.model_to_graphviz(model)`可视化模型结构
- 使用 `idata.to_netcdf('results.nc')` 保存结果
- 使用 `az.from_netcdf('results.nc')` 加载
- 对于非常大的模型，请考虑小批量 ADVI 或数据子采样