<!-- 此文件由机器翻译自 bayesian_statistics.md -->

# 贝叶斯统计分析

本文件提供了进行和解释贝叶斯统计分析的指导，贝叶斯统计分析为频率主义（经典）统计提供了替代框架。

## 贝叶斯哲学与频率主义哲学

### 根本差异

|方面|常客 |贝叶斯 |
|--------|-------------|----------|
| **概率解释** |事件的长期频率|置信度/不确定度 |
| **参数** |已修复但未知 |具有分布的随机变量 |
| **推论** |基于抽样分布|基于后验分布 |
| **主要输出** | p 值、置信区间 |后验概率、可信区间 |
| **先前信息** |未正式注册|通过先验明确纳入 |
| **假设检验** |拒绝/未能拒绝 null |给定数据的假设概率 |
| **样本大小** |通常需要最低限度|可以处理任何样本大小 |
| **解读** |间接（给定 H₀ 的数据概率）|直接（给定数据的假设概率）|

### 关键问题差异

**频率主义者**：“如果零假设为真，那么观察到如此极端或更极端的数据的概率是多少？”

**贝叶斯**：“根据观察到的数据，假设成立的概率是多少？”

贝叶斯问题更直观，直接解决研究人员想知道的问题。

---

## 贝叶斯定理

**公式**：
```
P(θ|D) = P(D|θ) × P(θ) / P(D)
```

**用言语来说**：
<<<代码块_1>>>

其中：
- **θ (theta)**：感兴趣的参数（例如，平均差、相关性）
- **D**：观察到的数据
- **P(θ|D)**：后验分布（看到数据后对 θ 的信念）
- **P(D|θ)**：可能性（给定 θ 的数据的概率）
- **P(θ)**：先验分布（在看到数据之前对 θ 的信念）
- **P(D)**：边际可能性/证据（标准化常数）

---

## 先验分布

### 先验的类型

#### 1. 信息先验

**何时使用**：当您拥有以下方面的丰富先验知识时：
- 以前的研究
- 专业知识
- 理论
- 试点数据

**示例**：荟萃分析显示效应大小 d ≈ 0.40，SD = 0.15
- 先前：正常（0.40，0.15）

**优点**：
- 整合现有知识
- 更高效（需要更小的样品）
- 可以用小数据稳定估计

**缺点**：
- 主观（但主观性可以是力量）
- 必须合理且透明
- 如果先验与数据存在强烈冲突，可能会引起争议

---

#### 2. 信息较弱的先验

**何时使用**：大多数应用程序的默认选择

**特点**：
- 规范化估计（防止极端值）
- 对中等数据的后验影响最小
- 防止计算问题

**先验示例**：
- 效果大小：正常(0, 1) 或柯西(0, 0.707)
- 方差：半柯西(0, 1)
- 相关性：Uniform(-1, 1) 或 Beta(2, 2)

**优点**：
- 平衡客观性和规范化
- 计算稳定
- 广泛接受

---

#### 3. 非信息性（平坦/统一）先验

**何时使用**：当试图“客观”时

**示例**：对于任何值，Uniform(-Infinity, Infinity)

**⚠️注意**：
- 可能导致后路不正确
- 可能会产生不合理的结果
- 并非真正的“非信息性”（仍然做出假设）
- 在现代贝叶斯实践中通常不推荐

**更好的选择**：使用信息量较弱的先验

---

### 先验敏感性分析

**始终进行**：测试结果如何随不同先验变化

**流程**：
1. 将模型与默认/事先计划的进行拟合
2. 拟合具有更分散先验的模型
3. 拟合更集中的先验模型
4. 比较后验分布

**报告**：
- 如果结果相似：证据是强有力的
- 如果结果差异很大：数据不足以压倒先前的结果

**Python 示例**：
<<<代码块_2>>>

---

## 贝叶斯假设检验

### 贝叶斯因子 (BF)

**它是什么**：两个相互竞争的假设的证据比率

**公式**：
<<<代码块_3>>>

**解释**：

| BF₁₀ |证据|
|------|----------|
| >100 |对 H₁ 具有决定性作用 |
| 30-100 | H₁ 非常强 |
| 10-30 | 10-30 H₁ 强 |
| 3-10 | 3-10 H₁ 中等 |
| 1-3 | 1-3 H₁ 轶事 |
| 1 |没有证据 |
| 1/3-1 | 1/3-1 H₀ 的轶事 |
| 1/10-1/3 | H₀ 中等 |
| 1/30-1/10 |强 H₀ |
| 1/100-1/30 | H₀ 非常强 |
| <1/100 |对 H₀ 具有决定性作用 |

**相对于 p 值的优势**：
1. 可以为原假设提供证据
2.不依赖于采样意图（不存在“偷看”问题）
3. 直接量化证据
4.可以更新更多数据

**Python计算**：
<<<代码块_4>>>

---

### 实用等效区域 (ROPE)

**目的**：定义可忽略效应大小的范围

**流程**：
1. 定义 ROPE（例如，d ∈ [-0.1, 0.1] 对于可忽略的影响）
2.计算ROPE内部后部的%
3. 做出决定：
   - >95% 的 ROPE：接受实际等效性
   - >95% 超出绳索：拒绝等效
   - 否则：不确定

**优点**：直接测试实际意义

**Python 示例**：
<<<代码块_5>>>

---

## 贝叶斯估计

### 可信区间

**它是什么**：包含 X% 概率参数的区间

**95% 可信区间解释**：
> “真实参数有 95% 的概率位于该区间内。”

**这就是人们认为置信区间的含义**（但在频率论框架中并非如此）

**类型**：

#### 等尾区间 (ETI)
- 第 2.5 至 97.5 个百分位
- 计算简单
- 可能不包括倾斜分布的模式

#### 最高密度区间 (HDI)
- 包含 95% 分布的最窄区间
- 始终包含模式
- 更适合偏态分布

**Python计算**：
<<<代码块_6>>>

---

### 后验分布

**解释后验分布**：

1. **集中趋势**：
   - 平均值：平均后验值
   - 中位数：第 50 个百分位
   - 众数：最可能值（MAP - 最大后验概率）

2. **不确定性**：
   - SD：后路扩散
   - 可信区间：量化不确定性

3. **形状**：
   - 对称：与正常相似
   - 倾斜：不对称的不确定性
   - 多模式：多个合理值

**可视化**：
```python
import matplotlib.pyplot as plt
import arviz as az

# Posterior plot with HDI
az.plot_posterior(trace, hdi_prob=0.95)

# Trace plot (check convergence)
az.plot_trace(trace)

# Forest plot (multiple parameters)
az.plot_forest(trace)
```

---

## 常见贝叶斯分析

### 贝叶斯 T 检验

**目的**：比较两组（贝叶斯替代 t 检验）

**输出**：
1.均差后验分布
2. 95%可信区间
3.贝叶斯因子（BF₁₀）
4. 方向假设的概率（例如，P(μ₁ > μ2)）

**Python 实现**：
```python
import pymc as pm
import arviz as az

# Bayesian independent samples t-test
with pm.Model() as model:
    # Priors for group means
    mu1 = pm.Normal('mu1', mu=0, sigma=10)
    mu2 = pm.Normal('mu2', mu=0, sigma=10)

    # Prior for pooled standard deviation
    sigma = pm.HalfNormal('sigma', sigma=10)

    # Likelihood
    y1 = pm.Normal('y1', mu=mu1, sigma=sigma, observed=group1)
    y2 = pm.Normal('y2', mu=mu2, sigma=sigma, observed=group2)

    # Derived quantity: mean difference
    diff = pm.Deterministic('diff', mu1 - mu2)

    # Sample posterior
    trace = pm.sample(2000, tune=1000, return_inferencedata=True)

# Analyze results
print(az.summary(trace, var_names=['mu1', 'mu2', 'diff']))

# Probability that group1 > group2
prob_greater = np.mean(trace.posterior['diff'].values > 0)
print(f"P(μ₁ > μ₂) = {prob_greater:.3f}")

# Plot posterior
az.plot_posterior(trace, var_names=['diff'], ref_val=0)
```

---

### 贝叶斯方差分析

**目的**：比较三个或更多组

**型号**：
```python
import pymc as pm

with pm.Model() as anova_model:
    # Hyperpriors
    mu_global = pm.Normal('mu_global', mu=0, sigma=10)
    sigma_between = pm.HalfNormal('sigma_between', sigma=5)
    sigma_within = pm.HalfNormal('sigma_within', sigma=5)

    # Group means (hierarchical)
    group_means = pm.Normal('group_means',
                            mu=mu_global,
                            sigma=sigma_between,
                            shape=n_groups)

    # Likelihood
    y = pm.Normal('y',
                  mu=group_means[group_idx],
                  sigma=sigma_within,
                  observed=data)

    trace = pm.sample(2000, tune=1000, return_inferencedata=True)

# Posterior contrasts
contrast_1_2 = trace.posterior['group_means'][:,:,0] - trace.posterior['group_means'][:,:,1]
```

---

### 贝叶斯相关性

**目的**：估计两个变量之间的相关性

**优点**：提供相关值的分布

**Python 实现**：
```python
import pymc as pm

with pm.Model() as corr_model:
    # Prior on correlation
    rho = pm.Uniform('rho', lower=-1, upper=1)

    # Convert to covariance matrix
    cov_matrix = pm.math.stack([[1, rho],
                                [rho, 1]])

    # Likelihood (bivariate normal)
    obs = pm.MvNormal('obs',
                     mu=[0, 0],
                     cov=cov_matrix,
                     observed=np.column_stack([x, y]))

    trace = pm.sample(2000, tune=1000, return_inferencedata=True)

# Summarize correlation
print(az.summary(trace, var_names=['rho']))

# Probability that correlation is positive
prob_positive = np.mean(trace.posterior['rho'].values > 0)
```

---

### 贝叶斯线性回归

**目的**：预测变量和结果之间的模型关系

**优点**：
- 所有参数的不确定性
- 自然正则化（通过先验）
- 可以结合先验知识
- 预测的可信区间

**Python 实现**：
```python
import pymc as pm

with pm.Model() as regression_model:
    # Priors for coefficients
    alpha = pm.Normal('alpha', mu=0, sigma=10)  # Intercept
    beta = pm.Normal('beta', mu=0, sigma=10, shape=n_predictors)
    sigma = pm.HalfNormal('sigma', sigma=10)

    # Expected value
    mu = alpha + pm.math.dot(X, beta)

    # Likelihood
    y_obs = pm.Normal('y_obs', mu=mu, sigma=sigma, observed=y)

    trace = pm.sample(2000, tune=1000, return_inferencedata=True)

# Posterior predictive checks
with regression_model:
    ppc = pm.sample_posterior_predictive(trace)

az.plot_ppc(ppc)

# Predictions with uncertainty
with regression_model:
    pm.set_data({'X': X_new})
    posterior_pred = pm.sample_posterior_predictive(trace)
```

---

## 分层（多级）模型

**何时使用**：
- 嵌套/集群数据（学校内的学生）
- 重复措施
- 荟萃分析
- 不同群体的效果不同

**关键概念**：部分池化
- 完全池化：忽略组（有偏见）
- 无合并：单独分析组（高方差）
- 部分池化：跨组借用力量（贝叶斯）

**示例：不同的截距**：
```python
with pm.Model() as hierarchical_model:
    # Hyperpriors
    mu_global = pm.Normal('mu_global', mu=0, sigma=10)
    sigma_between = pm.HalfNormal('sigma_between', sigma=5)
    sigma_within = pm.HalfNormal('sigma_within', sigma=5)

    # Group-level intercepts
    alpha = pm.Normal('alpha',
                     mu=mu_global,
                     sigma=sigma_between,
                     shape=n_groups)

    # Likelihood
    y_obs = pm.Normal('y_obs',
                     mu=alpha[group_idx],
                     sigma=sigma_within,
                     observed=y)

    trace = pm.sample()
```

---

## 型号对比

### 方法

#### 1. 贝叶斯因子
- 直接比较模型证据
- 对先前的规格敏感
- 可能是计算密集型的

#### 2. 信息标准

**WAIC（广泛适用的信息标准）**：
- AIC 的贝叶斯模拟
- 越低越好
- 考虑参数的有效数量

**LOO（留一交叉验证）**：
- 估计样本外预测误差
- 越低越好
- 比WAIC更稳健

**Python计算**：
```python
import arviz as az

# Calculate WAIC and LOO
waic = az.waic(trace)
loo = az.loo(trace)

print(f"WAIC: {waic.elpd_waic:.2f}")
print(f"LOO: {loo.elpd_loo:.2f}")

# Compare multiple models
comparison = az.compare({
    'model1': trace1,
    'model2': trace2,
    'model3': trace3
})
print(comparison)
```

---

## 检查贝叶斯模型

### 1. 收敛诊断

**R-hat（Gelman-Rubin 统计）**：
- 比较链内和链间方差
- 接近 1.0 的值表示收敛
- R 帽 < 1.01：好
- R-hat > 1.05：收敛性差

**有效样本量 (ESS)**：
- 独立样本数
- 越高越好
- 建议每个链 ESS > 400

**跟踪图**：
- 看起来应该像“毛毛虫”
- 没有趋势，没有卡住的链条

**Python 检查**：
```python
# Automatic summary with diagnostics
print(az.summary(trace, var_names=['parameter']))

# Visual diagnostics
az.plot_trace(trace)
az.plot_rank(trace)  # Rank plots
```

---

### 2. 事后预测检查

**目的**：模型生成的数据是否与观察到的数据相似？

**流程**：
1. 从后验生成预测
2、与实际数据对比
3.寻找系统性差异

**Python 实现**：
```python
with model:
    ppc = pm.sample_posterior_predictive(trace)

# Visual check
az.plot_ppc(ppc, num_pp_samples=100)

# Quantitative checks
obs_mean = np.mean(observed_data)
pred_means = [np.mean(sample) for sample in ppc.posterior_predictive['y_obs']]
p_value = np.mean(pred_means >= obs_mean)  # Bayesian p-value
```

---

## 报告贝叶斯结果

### T 检验报告示例
> “进行贝叶斯独立样本 t 检验来比较 A 组和 B 组。使用弱信息先验：均值差使用正态 (0, 1)，合并标准差使用 Half-Cauchy(0, 1)。均值差的后验分布平均值为 5.2 (95% CI [2.3, 8.1])，表明 A 组得分高于 B 组。贝叶斯因子 BF₁₀ = 23.5 为各组之间的差异提供了强有力的证据，A 组的平均值超过 B 组平均值的概率为 99.7%。”

### 回归报告示例

> “贝叶斯线性回归拟合了信息量较弱的先验（系数为正态 (0, 10)，残差 SD 为 Half-Cauchy(0, 5)）。该模型解释了显着方差（R² = 0.47, 95% CI [0.38, 0.55]）。学习时间（β = 0.52, 95% CI [0.38, 0.66]）和先前的 GPA (β = 0.31，95% CI [0.17, 0.45]）是可靠的预测因子（95% CI 排除零）。

---

## 优点和局限性

### 优势

1. **直观解释**：关于参数的直接概率陈述
2. **结合先验知识**：使用所有可用信息
3. **灵活**：轻松处理复杂模型
4. **无 p-hacking**：可以在数据到达时查看数据
5. **量化不确定性**：完全后验分布
6. **小样本**：适用于任何样本大小

### 限制

1. **计算**：需要MCMC采样（可能很慢）
2. **先前的规范**：需要思考和论证
3. **复杂性**：学习曲线更陡
4. **软件**：比频率论方法更少的工具
5. **沟通**：可能需要教育审稿人/读者

---

## 关键的 Python 包

- **PyMC**：完整的贝叶斯建模框架
- **ArviZ**：可视化和诊断
- **Bambi**：回归模型的高级接口
- **PyStan**：Stan 的 Python 接口
- **TensorFlow Probability**：使用 TensorFlow 进行贝叶斯推理

---

## 何时使用贝叶斯方法

**当**时使用贝叶斯：
- 您有要合并的先前信息
- 你想要直接的概率陈述
- 样本量较小
- 模型复杂（层次结构、缺失数据等）
- 您想要在数据到达时更新分析

**在以下情况下，Frequentist 可能就足够了：
- 大样本标准分析
- 没有先验信息
- 计算资源有限
- 审稿人不熟悉贝叶斯方法