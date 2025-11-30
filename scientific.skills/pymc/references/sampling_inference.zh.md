<!-- 此文件由机器翻译自 sampling_inference.md -->

# PyMC 采样和推理方法

本参考文献涵盖了 PyMC 中可用于后验推理的采样算法和推理方法。

## MCMC 采样方法

### 初级采样函数

**`pm.sample(draws=1000, tune=1000, chains=4, **kwargs)`**

PyMC 中 MCMC 采样的主界面。

**关键参数：**
- `draws`：每个链绘制的样本数（默认值：1000）
- `tune`：调整/预热样本的数量（默认值：1000，丢弃）
- `chains`：并行链数量（默认值：4）
- `cores`：要使用的 CPU 核心数（默认值：全部可用）
- `target_accept`：步长调整的目标接受率（默认值：0.8，对于困难的后验增加到 0.9-0.95）
- `random_seed`：用于再现性的随机种子
- `return_inferencedata`：返回 ArviZ InferenceData 对象（默认值：True）
- `idata_kwargs`：用于创建 InferenceData 的附加 kwargs（例如，`{"log_likelihood": True}` 用于模型比较）

**返回：** InferenceData 对象包含后验样本、采样统计数据和诊断信息

**示例：**
```python
with pm.Model() as model:
    # ... define model ...
    idata = pm.sample(draws=2000, tune=1000, chains=4, target_accept=0.9)
```

### 采样算法

PyMC 根据模型结构自动选择合适的采样器，但您可以手动指定算法。

#### NUTS（禁止掉头采样器）

**连续参数的默认算法**。高效哈密顿蒙特卡罗变体。

- 自动调整步长和质量矩阵
- 自适应：在调整过程中探索后验几何
- 最适合平滑、连续的后部
- 可能会与高相关性或多模态作斗争

**说明书规格：**
<<<代码块_1>>>

**何时调整：**
- 如果发现差异，则增加 `target_accept` (0.9-0.99)
- 使用 `init='adapt_diag'` 进行更快的初始化（默认）
- 使用 `init='jitter+adapt_diag'` 进行困难的初始化

#### 大都会

通用 Metropolis-Hastings 采样器。

- 适用于连续变量和离散变量
- 平滑连续后部的效率低于 NUTS
- 对于离散参数或不可微分模型有用
- 需要手动调整

**示例：**
<<<代码块_2>>>

#### 切片采样器

单变量分布的切片采样。

- 无需调整
- 适用于困难的单变量后验
- 对于高维度可能会很慢

**示例：**
<<<代码块_3>>>

#### 复合步骤

结合不同的采样器以获得不同的参数。

**示例：**
<<<代码块_4>>>

### 采样诊断

PyMC 自动计算诊断。在信任结果之前检查这些：

#### 有效样本量 (ESS)

测量相关样本中的独立信息。

- **经验法则**：ESS > 每个链 400（4 个链总共 1600）
- 低 ESS 表示高自相关
- 通过以下方式访问：`az.ess(idata)`

#### R-hat（Gelman-Rubin 统计）

衡量跨链的收敛性。

- **经验法则**：所有参数的 R-hat < 1.01
- R-hat > 1.01 表示不收敛
- 通过以下方式访问：`az.rhat(idata)`

#### 分歧

指出 NUTS 陷入困境的地区。

- **经验法则**：0 分歧（或很少）
- 差异表明样本存在偏差
- **修复**：增加 `target_accept`、重新参数化或使用更强的先验
- 通过以下方式访问：`idata.sample_stats.diverging.sum()`

#### 能量图

可视化哈密顿蒙特卡罗能量转变。

<<<代码块_5>>>

能量分布之间良好的分离表明健康的采样。

### 处理采样问题

#### 分歧

<<<代码块_6>>>

#### 慢速采样

```python
# Use fewer tuning steps if model is simple
idata = pm.sample(tune=500)

# Increase cores for parallelization
idata = pm.sample(cores=8, chains=8)

# Use variational inference for initialization
with model:
    approx = pm.fit()  # Run ADVI
    idata = pm.sample(start=approx.sample(return_inferencedata=False)[0])
```

#### 高自相关

```python
# Increase draws
idata = pm.sample(draws=5000)

# Reparameterize to reduce correlation
# Consider using QR decomposition for regression models
```

## 变分推理

更快地进行大型模型的近似推理或快速探索。

### ADVI（自动微分变分推理）

**`pm.fit(n=10000, method='advi', **kwargs)`>**

使用更简单的分布（通常是平均场高斯分布）来近似后验。

**关键参数：**
- `n`：迭代次数（默认值：10000）
- `method`：VI算法（'advi'，'fullrank_advi'，'svgd'）
- `random_seed`：随机种子

**返回：** 用于采样和分析的近似对象

**示例：**
```python
with model:
    approx = pm.fit(n=50000)
    # Draw samples from approximation
    idata = approx.sample(1000)
    # Or sample for MCMC initialization
    start = approx.sample(return_inferencedata=False)[0]
```

**权衡：**
- **优点**：比 MCMC 快得多，可扩展到大数据
- **缺点**：近似值，可能会错过后部结构，低估不确定性

### 全等级 ADVI

捕获参数之间的相关性。

```python
with model:
    approx = pm.fit(method='fullrank_advi')
```

比平均场更准确，但速度更慢。

### SVGD（斯坦因变分梯度下降）

非参数变分推理。

```python
with model:
    approx = pm.fit(method='svgd', n=20000)
```
更好地捕获多模态，但计算成本更高。

## 事前和事后预测采样

### 事先预测采样

来自先验分布的样本（在查看数据之前）。

**`pm.sample_prior_predictive(samples=500, **kwargs)`>**

**目的：**
- 验证先验是否合理
- 在拟合之前检查隐含的预测
- 确保模型生成合理的数据

**示例：**
```python
with model:
    prior_pred = pm.sample_prior_predictive(samples=1000)

# Visualize prior predictions
az.plot_ppc(prior_pred, group='prior')
```

### 后验预测采样

来自后验预测分布的样本（拟合后）。

**`pm.sample_posterior_predictive(trace, **kwargs)`>**

**目的：**
- 通过后验预测检查进行模型验证
- 生成新数据的预测
- 评估拟合优度

**示例：**
```python
with model:
    # After sampling
    idata = pm.sample()

    # Add posterior predictive samples
    pm.sample_posterior_predictive(idata, extend_inferencedata=True)

# Posterior predictive check
az.plot_ppc(idata)
```

### 新数据的预测

更新数据和样本预测分布：

```python
with model:
    # Original model fit
    idata = pm.sample()

    # Update with new predictor values
    pm.set_data({'X': X_new})

    # Sample predictions
    post_pred_new = pm.sample_posterior_predictive(
        idata.posterior,
        var_names=['y_pred']
    )
```

## 最大后验 (MAP) 估计

求后验模式（点估计）。

**`pm.find_MAP(start=None, method='L-BFGS-B', **kwargs)`>**

**何时使用：**
- 快速点估计
- MCMC的初始化
- 当不需要完整的后路时

**示例：**
```python
with model:
    map_estimate = pm.find_MAP()
    print(map_estimate)
```

**限制：**
- 不量化不确定性
- 可以在多模态后验中找到局部最优
- 对先前的规格敏感

## 推理建议

### 标准工作流程

1. **从 ADVI 开始**进行快速探索：
   ```python
   approx = pm.fit(n=20000)
   ```

2. **运行 MCMC** 进行完整推理：
   ```python
   idata = pm.sample(draws=2000, tune=1000)
   ```

3. **检查诊断**：
   ```python
   az.summary(idata, var_names=['~mu_log__'])  # Exclude transformed vars
   ```

4. **后验预测示例**：
   ```python
   pm.sample_posterior_predictive(idata, extend_inferencedata=True)
   ```

### 选择推理方法

|场景|推荐方法 |
|----------|--------------------|
|中小型模型，需要充分的不确定性| MCMC 与坚果 |
|大模型，初步探索|安维维|
|离散参数|都市化还是边缘化？
|具有分歧的层次模型|非中心参数化 + NUTS |
|数据量非常大|小批量 ADVI |
|快速点估计| MAP 或 ADVI |

### 重新参数化技巧

**分层模型的非中心参数化**：

```python
# Centered (can cause divergences):
mu = pm.Normal('mu', 0, 10)
sigma = pm.HalfNormal('sigma', 1)
theta = pm.Normal('theta', mu, sigma, shape=n_groups)

# Non-centered (better sampling):
mu = pm.Normal('mu', 0, 10)
sigma = pm.HalfNormal('sigma', 1)
theta_offset = pm.Normal('theta_offset', 0, 1, shape=n_groups)
theta = pm.Deterministic('theta', mu + sigma * theta_offset)
```

**相关预测变量的 QR 分解**：

```python
import numpy as np

# QR decomposition
Q, R = np.linalg.qr(X)

with pm.Model():
    # Uncorrelated coefficients
    beta_tilde = pm.Normal('beta_tilde', 0, 1, shape=p)

    # Transform back to original scale
    beta = pm.Deterministic('beta', pm.math.solve(R, beta_tilde))

    mu = pm.math.dot(Q, beta_tilde)
    sigma = pm.HalfNormal('sigma', 1)
    y = pm.Normal('y', mu, sigma, observed=y_obs)
```

## 高级采样

### 顺序蒙特卡罗 (SMC)

对于复杂的后验或模型证据估计：

```python
with model:
    idata = pm.sample_smc(draws=2000, chains=4)
```

适合多模态后路或当 NUTS 遇到困难时。

### 自定义初始化

提供起始值：

```python
start = {'mu': 0, 'sigma': 1}
with model:
    idata = pm.sample(start=start)
```

或者使用 MAP 估计：

```python
with model:
    start = pm.find_MAP()
    idata = pm.sample(start=start)
```