<!-- 此文件由机器翻译自 workflows.md -->

# PyMC 工作流程和常见模式

该参考提供了在 PyMC 中构建、验证和分析贝叶斯模型的标准工作流程和模式。

## 标准贝叶斯工作流程

### 完整的工作流程模板

```python
import pymc as pm
import arviz as az
import numpy as np
import matplotlib.pyplot as plt

# 1. PREPARE DATA
# ===============
X = ...  # Predictor variables
y = ...  # Observed outcomes

# Standardize predictors for better sampling
X_scaled = (X - X.mean(axis=0)) / X.std(axis=0)

# 2. BUILD MODEL
# ==============
with pm.Model() as model:
    # Define coordinates for named dimensions
    coords = {
        'predictors': ['var1', 'var2', 'var3'],
        'obs_id': np.arange(len(y))
    }

    # Priors
    alpha = pm.Normal('alpha', mu=0, sigma=1)
    beta = pm.Normal('beta', mu=0, sigma=1, dims='predictors')
    sigma = pm.HalfNormal('sigma', sigma=1)

    # Linear predictor
    mu = alpha + pm.math.dot(X_scaled, beta)

    # Likelihood
    y_obs = pm.Normal('y_obs', mu=mu, sigma=sigma, observed=y, dims='obs_id')

# 3. PRIOR PREDICTIVE CHECK
# ==========================
with model:
    prior_pred = pm.sample_prior_predictive(samples=1000, random_seed=42)

# Visualize prior predictions
az.plot_ppc(prior_pred, group='prior', num_pp_samples=100)
plt.title('Prior Predictive Check')
plt.show()

# 4. FIT MODEL
# ============
with model:
    # Quick VI exploration (optional)
    approx = pm.fit(n=20000, random_seed=42)

    # Full MCMC inference
    idata = pm.sample(
        draws=2000,
        tune=1000,
        chains=4,
        target_accept=0.9,
        random_seed=42,
        idata_kwargs={'log_likelihood': True}  # For model comparison
    )

# 5. CHECK DIAGNOSTICS
# ====================
# Summary statistics
print(az.summary(idata, var_names=['alpha', 'beta', 'sigma']))

# R-hat and ESS
summary = az.summary(idata)
if (summary['r_hat'] > 1.01).any():
    print("WARNING: Some R-hat values > 1.01, chains may not have converged")

if (summary['ess_bulk'] < 400).any():
    print("WARNING: Some ESS values < 400, consider more samples")

# Check divergences
divergences = idata.sample_stats.diverging.sum().item()
print(f"Number of divergences: {divergences}")

# Trace plots
az.plot_trace(idata, var_names=['alpha', 'beta', 'sigma'])
plt.tight_layout()
plt.show()

# 6. POSTERIOR PREDICTIVE CHECK
# ==============================
with model:
    pm.sample_posterior_predictive(idata, extend_inferencedata=True, random_seed=42)

# Visualize fit
az.plot_ppc(idata, num_pp_samples=100)
plt.title('Posterior Predictive Check')
plt.show()

# 7. ANALYZE RESULTS
# ==================
# Posterior distributions
az.plot_posterior(idata, var_names=['alpha', 'beta', 'sigma'])
plt.tight_layout()
plt.show()

# Forest plot for coefficients
az.plot_forest(idata, var_names=['beta'], combined=True)
plt.title('Coefficient Estimates')
plt.show()

# 8. PREDICTIONS FOR NEW DATA
# ============================
X_new = ...  # New predictor values
X_new_scaled = (X_new - X.mean(axis=0)) / X.std(axis=0)

with model:
    # Update data
    pm.set_data({'X': X_new_scaled})

    # Sample predictions
    post_pred = pm.sample_posterior_predictive(
        idata.posterior,
        var_names=['y_obs'],
        random_seed=42
    )

# Prediction intervals
y_pred_mean = post_pred.posterior_predictive['y_obs'].mean(dim=['chain', 'draw'])
y_pred_hdi = az.hdi(post_pred.posterior_predictive, var_names=['y_obs'])

# 9. SAVE RESULTS
# ===============
idata.to_netcdf('model_results.nc')  # Save for later
```

## 模型构建模式

### 线性回归

<<<代码块_1>>>

### 逻辑回归

<<<代码块_2>>>

### 分层/多级模型

<<<代码块_3>>>

### 泊松回归（计数数据）

<<<代码块_4>>>

### 时间序列（自回归）

<<<代码块_5>>>

### 混合模型

<<<代码块_6>>>

## 数据准备最佳实践

### 标准化

标准化连续预测变量以实现更好的采样：

```python
# Standardize
X_mean = X.mean(axis=0)
X_std = X.std(axis=0)
X_scaled = (X - X_mean) / X_std

# Model with scaled data
with pm.Model() as model:
    beta_scaled = pm.Normal('beta_scaled', 0, 1)
    # ... rest of model ...

# Transform back to original scale
beta_original = beta_scaled / X_std
alpha_original = alpha - (beta_scaled * X_mean / X_std).sum()
```

### 处理缺失数据

将缺失值视为参数：

```python
# Identify missing values
missing_idx = np.isnan(X)
X_observed = np.where(missing_idx, 0, X)  # Placeholder

with pm.Model() as model:
    # Prior for missing values
    X_missing = pm.Normal('X_missing', mu=0, sigma=1, shape=missing_idx.sum())

    # Combine observed and imputed
    X_complete = pm.math.switch(missing_idx.flatten(), X_missing, X_observed.flatten())

    # ... rest of model using X_complete ...
```

### 居中和缩放

对于回归模型、中心预测变量和结果：

```python
# Center
X_centered = X - X.mean(axis=0)
y_centered = y - y.mean()

with pm.Model() as model:
    # Simpler prior on intercept
    alpha = pm.Normal('alpha', mu=0, sigma=1)  # Intercept near 0 when centered
    beta = pm.Normal('beta', mu=0, sigma=1, shape=n_predictors)

    mu = alpha + pm.math.dot(X_centered, beta)
    sigma = pm.HalfNormal('sigma', sigma=1)

    y_obs = pm.Normal('y_obs', mu=mu, sigma=sigma, observed=y_centered)
```

## 事先选择指南

### 信息较弱的先验

当您的先验知识有限时使用：

```python
# For standardized predictors
beta = pm.Normal('beta', mu=0, sigma=1)

# For scale parameters
sigma = pm.HalfNormal('sigma', sigma=1)

# For probabilities
p = pm.Beta('p', alpha=2, beta=2)  # Slight preference for middle values
```

### 信息先验

使用领域知识：

```python
# Effect size from literature: Cohen's d ≈ 0.3
beta = pm.Normal('beta', mu=0.3, sigma=0.1)

# Physical constraint: probability between 0.7-0.9
p = pm.Beta('p', alpha=8, beta=2)  # Check with prior predictive!
```

### 事先预测检查

始终验证先验：

```python
with model:
    prior_pred = pm.sample_prior_predictive(samples=1000)

# Check if predictions are reasonable
print(f"Prior predictive range: {prior_pred.prior_predictive['y'].min():.2f} to {prior_pred.prior_predictive['y'].max():.2f}")
print(f"Observed range: {y_obs.min():.2f} to {y_obs.max():.2f}")

# Visualize
az.plot_ppc(prior_pred, group='prior')
```

## 模型比较工作流程

### 比较多个模型

```python
import arviz as az

# Fit multiple models
models = {}
idatas = {}

# Model 1: Simple linear
with pm.Model() as models['linear']:
    # ... define model ...
    idatas['linear'] = pm.sample(idata_kwargs={'log_likelihood': True})

# Model 2: With interaction
with pm.Model() as models['interaction']:
    # ... define model ...
    idatas['interaction'] = pm.sample(idata_kwargs={'log_likelihood': True})

# Model 3: Hierarchical
with pm.Model() as models['hierarchical']:
    # ... define model ...
    idatas['hierarchical'] = pm.sample(idata_kwargs={'log_likelihood': True})

# Compare using LOO
comparison = az.compare(idatas, ic='loo')
print(comparison)

# Visualize comparison
az.plot_compare(comparison)
plt.show()

# Check LOO reliability
for name, idata in idatas.items():
    loo = az.loo(idata, pointwise=True)
    high_pareto_k = (loo.pareto_k > 0.7).sum().item()
    if high_pareto_k > 0:
        print(f"Warning: {name} has {high_pareto_k} observations with high Pareto-k")
```

### 模型重量

```python
# Get model weights (pseudo-BMA)
weights = comparison['weight'].values

print("Model probabilities:")
for name, weight in zip(comparison.index, weights):
    print(f"  {name}: {weight:.2%}")

# Model averaging (weighted predictions)
def weighted_predictions(idatas, weights):
    preds = []
    for (name, idata), weight in zip(idatas.items(), weights):
        pred = idata.posterior_predictive['y_obs'].mean(dim=['chain', 'draw'])
        preds.append(weight * pred)
    return sum(preds)

averaged_pred = weighted_predictions(idatas, weights)
```

## 诊断和故障排除

### 诊断采样问题

```python
def diagnose_sampling(idata, var_names=None):
    """Comprehensive sampling diagnostics"""

    # Check convergence
    summary = az.summary(idata, var_names=var_names)

    print("=== Convergence Diagnostics ===")
    bad_rhat = summary[summary['r_hat'] > 1.01]
    if len(bad_rhat) > 0:
        print(f"⚠️  {len(bad_rhat)} variables with R-hat > 1.01")
        print(bad_rhat[['r_hat']])
    else:
        print("✓ All R-hat values < 1.01")

    # Check effective sample size
    print("\n=== Effective Sample Size ===")
    low_ess = summary[summary['ess_bulk'] < 400]
    if len(low_ess) > 0:
        print(f"⚠️  {len(low_ess)} variables with ESS < 400")
        print(low_ess[['ess_bulk', 'ess_tail']])
    else:
        print("✓ All ESS values > 400")

    # Check divergences
    print("\n=== Divergences ===")
    divergences = idata.sample_stats.diverging.sum().item()
    if divergences > 0:
        print(f"⚠️  {divergences} divergent transitions")
        print("   Consider: increase target_accept, reparameterize, or stronger priors")
    else:
        print("✓ No divergences")

    # Check tree depth
    print("\n=== NUTS Statistics ===")
    max_treedepth = idata.sample_stats.tree_depth.max().item()
    hits_max = (idata.sample_stats.tree_depth == max_treedepth).sum().item()
    if hits_max > 0:
        print(f"⚠️  Hit max treedepth {hits_max} times")
        print("   Consider: reparameterize or increase max_treedepth")
    else:
        print(f"✓ No max treedepth issues (max: {max_treedepth})")

    return summary

# Usage
diagnose_sampling(idata, var_names=['alpha', 'beta', 'sigma'])
```

### 常见修复

|问题 |解决方案 |
|---------|----------|
|分歧|增加`target_accept=0.95`，使用非居中参数化 |
|低 ESS |对更多抽奖进行采样，重新参数化以减少相关性 |
|高 R 帽 |运行更长的链，检查多模态，改进初始化 |
|慢采样|使用ADVI初始化，重新参数化，降低模型复杂度 |
|偏后 |检查先前的预测，确保可能性是正确的 |

## 使用命名尺寸（dims）

### 暗淡的好处

- 更具可读性的代码
- 更轻松的子集设置和分析
- 更好的 xarray 集成

```python
# Define coordinates
coords = {
    'predictors': ['age', 'income', 'education'],
    'groups': ['A', 'B', 'C'],
    'time': pd.date_range('2020-01-01', periods=100, freq='D')
}

with pm.Model(coords=coords) as model:
    # Use dims instead of shape
    beta = pm.Normal('beta', mu=0, sigma=1, dims='predictors')
    alpha = pm.Normal('alpha', mu=0, sigma=1, dims='groups')
    y = pm.Normal('y', mu=0, sigma=1, dims=['groups', 'time'], observed=data)

# After sampling, dimensions are preserved
idata = pm.sample()

# Easy subsetting
beta_age = idata.posterior['beta'].sel(predictors='age')
group_A = idata.posterior['alpha'].sel(groups='A')
```

## 保存和加载结果

```python
# Save InferenceData
idata.to_netcdf('results.nc')

# Load InferenceData
loaded_idata = az.from_netcdf('results.nc')

# Save model for later predictions
import pickle

with open('model.pkl', 'wb') as f:
    pickle.dump({'model': model, 'idata': idata}, f)

# Load model
with open('model.pkl', 'rb') as f:
    saved = pickle.load(f)
    model = saved['model']
    idata = saved['idata']
```