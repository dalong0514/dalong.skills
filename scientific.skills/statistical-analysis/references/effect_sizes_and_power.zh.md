<!-- 此文件由机器翻译自 effect_sizes_and_power.md -->

# 效果大小和功效分析

本文件提供了计算、解释和报告效应大小以及为研究计划进行功效分析的指导。

## 为什么效应量很重要

1. **统计意义≠实际意义**：p值仅表明效应是否存在，而不是效应有多大
2. **依赖于样本量**：样本量较大时，微不足道的影响会变得“显着”
3. **解释**：效应大小提供了幅度和实际重要性
4. **荟萃分析**：效应大小能够合并跨研究的结果
5. **功效分析**：确定样本量所需

**黄金法则**：始终与 p 值一起报告效应大小。

---

## 按分析类型划分的效应量

### T 检验和均值差

#### Cohen's d（标准化均值差）

**公式**：
- 独立组：d = (M₁ - M2) / SD_pooled
- 配对组：d = M_diff / SD_diff

**解释**（科恩，1988）：
- 小：|d| = 0.20
- 中：|d| = 0.50
- 大：|d| = 0.80

**依赖于上下文的解释**：
- 在教育领域：d = 0.40 是成功干预的典型值
- 在心理学中：d = 0.40 被认为是有意义的
- 在医学上：小效应可能具有临床重要性

**Python计算**：
```python
import pingouin as pg
import numpy as np

# Independent t-test with effect size
result = pg.ttest(group1, group2, correction=False)
cohens_d = result['cohen-d'].values[0]

# Manual calculation
mean_diff = np.mean(group1) - np.mean(group2)
pooled_std = np.sqrt((np.var(group1, ddof=1) + np.var(group2, ddof=1)) / 2)
cohens_d = mean_diff / pooled_std

# Paired t-test
result = pg.ttest(pre, post, paired=True)
cohens_d = result['cohen-d'].values[0]
```

**d 的置信区间**：
<<<代码块_1>>>

---

#### Hedges' g（偏差校正 d）

**为什么使用它**：Cohen 的 d 在小样本 (n < 20) 下有轻微向上偏差

**公式**：g = d × Correction_Factor，其中 Correction_Factor = 1 - 3/(4df - 1)

**Python计算**：
<<<代码块_2>>>

**在以下情况下使用 Hedges' g：
- 样本量较小（每组 n < 20）
- 进行荟萃分析（荟萃分析标准）

---

#### 玻璃的 Δ (Delta)

**何时使用**：当一组是具有已知变异性的对照时

**公式**：Δ = (M₁ - M2) / SD_control

**用例**：
- 临床试验（使用对照组SD）
- 当治疗影响变异性时

---

### 方差分析

#### Eta 平方 (η²)

**它衡量什么**：由因子解释的总方差的比例

**公式**： η² = SS_effect / SS_total

**解释**：
- 小：η² = 0.01（1% 方差）
- 中等：η² = 0.06（6% 方差）
- 大：η² = 0.14（14% 方差）

**限制**：存在多个因素的偏差（总和 > 1.0）

**Python计算**：
<<<代码块_3>>>

---

#### 部分 Eta 平方 (η²_p)

**它衡量什么**：由因素解释的方差比例，不包括其他因素

**公式**： η²_p = SS_effect / (SS_effect + SS_error)

**解释**：与 η² 相同的基准

**何时使用**：多因子方差分析（因子设计中的标准）

**Python计算**：
<<<代码块_4>>>

---

#### 欧米伽平方 (ω²)

**它衡量什么**：解释总体方差的较少偏差估计

**为什么使用它**： η² 高估了效应大小； ω² 提供更好的总体估计

**公式**： ω² = (SS_effect - df_effect × MS_error) / (SS_total + MS_error)

**解释**：与 η² 相同的基准，但通常值较小

**Python计算**：
<<<代码块_5>>>

---

#### 科恩的 f

**测量的内容**：方差分析的效应大小（类似于 Cohen 的 d）

**公式**：f = √(η² / (1 - η²))

**解释**：
- 小：f = 0.10
- 中：f = 0.25
- 大：f = 0.40

**Python计算**：
<<<代码块_6>>>

**用于功效分析**：ANOVA 功效计算所需

---

### 相关性

#### Pearson 的 r / Spearman 的 ρ

**解释**：
- 小：|r| = 0.10
- 中：|r| = 0.30
- 大：|r| = 0.50

**重要说明**：
- r² = 决定系数（解释的方差比例）
- r = 0.30 表示 9% 共享方差 (0.30² = 0.09)
- 考虑方向（正/负）和背景

**Python计算**：
```python
import pingouin as pg

# Pearson correlation with CI
result = pg.corr(x, y, method='pearson')
r = result['r'].values[0]
ci = [result['CI95%'][0][0], result['CI95%'][0][1]]

# Spearman correlation
result = pg.corr(x, y, method='spearman')
rho = result['r'].values[0]
```

---

### 回归

#### R²（决定系数）

**它测量什么**：模型解释的 Y 方差比例

**解释**：
- 小：R² = 0.02
- 中：R² = 0.13
- 大：R² = 0.26

**取决于上下文**：
- 物理科学：R² > 0.90 预期
- 社会科学：R² > 0.30 视为良好
- 行为预测：R² > 0.10 可能有意义

**Python计算**：
```python
from sklearn.metrics import r2_score
from statsmodels.api import OLS

# Using statsmodels
model = OLS(y, X).fit()
r_squared = model.rsquared
adjusted_r_squared = model.rsquared_adj

# Manual
r_squared = 1 - (SS_residual / SS_total)
```

---

#### 调整后的 R²

**为什么使用它**：添加预测变量时，R² 会人为地增加；调整后的 R² 会惩罚模型复杂性

**公式**：R²_adj = 1 - (1 - R²) × (n - 1) / (n - k - 1)

**何时使用**：始终与 R² 一起报告多元回归

---

#### 标准化回归系数 (β)
**它衡量什么**：预测变量的一标准差变化对结果的影响（以标准差单位表示）

**解释**：类似于科恩的d
- 小：|β| = 0.10
- 中：|β| = 0.30
- 大：|β| = 0.50

**Python计算**：
```python
from scipy import stats

# Standardize variables first
X_std = (X - X.mean()) / X.std()
y_std = (y - y.mean()) / y.std()

model = OLS(y_std, X_std).fit()
beta = model.params
```

---

#### f²（科恩回归的 f 平方）

**它衡量什么**：单个预测变量或模型比较的效应大小

**公式**：f² = R²_AB - R²_A / (1 - R²_AB)

其中：
- R²_AB = R² 对于带有预测器的完整模型
- R²_A = R² 对于没有预测器的简化模型

**解释**：
- 小：f² = 0.02
- 中：f² = 0.15
- 大：f² = 0.35

**Python计算**：
```python
# Compare two nested models
model_full = OLS(y, X_full).fit()
model_reduced = OLS(y, X_reduced).fit()

r2_full = model_full.rsquared
r2_reduced = model_reduced.rsquared

f_squared = (r2_full - r2_reduced) / (1 - r2_full)
```

---

### 分类数据分析

#### 克莱默的 V

**它测量什么**：χ2 测试的关联强度（适用于任何表大小）

**公式**：V = √(χ² / (n × (k - 1)))

其中 k = min(行、列)

**解释**（对于 k > 2）：
- 小：V = 0.07
- 中：V = 0.21
- 大：V = 0.35

**对于 2×2 表**：使用 phi 系数 (φ)

**Python计算**：
```python
from scipy.stats.contingency import association

# Cramér's V
cramers_v = association(contingency_table, method='cramer')

# Phi coefficient (for 2x2)
phi = association(contingency_table, method='pearson')
```

---

#### 优势比 (OR) 和风险比 (RR)

**对于 2×2 列联表**：

|           |结果 + |结果- |
|------------|------------|------------|
|暴露|一个 |乙|
|未曝光 | c | d |

**优势比**：OR = (a/b) / (c/d) = ad / bc

**解释**：
- OR = 1：无关联
- OR > 1：正相关（几率增加）
- OR < 1：负关联（几率降低）
- OR = 2：两倍的几率
- OR = 0.5：一半的可能性

**风险比率**：RR = (a/(a+b)) / (c/(c+d))

**何时使用**：
- 队列研究：使用 RR（更容易解释）
- 病例对照研究：使用 OR（RR 不可用）
- 逻辑回归：OR是自然输出

**Python计算**：
```python
import statsmodels.api as sm

# From contingency table
odds_ratio = (a * d) / (b * c)

# Confidence interval
table = np.array([[a, b], [c, d]])
oddsratio, pvalue = stats.fisher_exact(table)

# From logistic regression
model = sm.Logit(y, X).fit()
odds_ratios = np.exp(model.params)  # Exponentiate coefficients
ci = np.exp(model.conf_int())  # Exponentiate CIs
```

---

### 贝叶斯效应大小

#### 贝叶斯因子 (BF)

**它衡量什么**：替代假设与零假设的证据比率

**解释**：
- BF₁₀ = 1：H₁ 和 H₀ 的证据相同
- BF₁₀ = 3：H₁ 的可能性是 H₀ 的 3 倍（中等证据）
- BF₁₀ = 10：H₁ 的可能性是 H₀ 的 10 倍（强有力的证据）
- BF₁₀ = 100：H₁ 的可能性比 H₀ 高 100 倍（决定性证据）
- BF₁₀ = 0.33：H₀ 的可能性是 H₁ 的 3 倍
- BF₁₀ = 0.10：H₀ 的可能性是 H₁ 的 10 倍

**分类**（Jeffreys，1961）：
- 1-3：轶事证据
- 3-10：中等证据
- 10-30：强有力的证据
- 30-100：非常有力的证据
- >100：决定性证据

**Python计算**：
```python
import pingouin as pg

# Bayesian t-test
result = pg.ttest(group1, group2, correction=False)
# Note: pingouin doesn't include BF; use other packages

# Using JASP or BayesFactor (R) via rpy2
# Or implement using numerical integration
```

---

## 功率分析

### 概念

**统计功效**：检测到效果（如果存在）的概率 (1 - β)

**常规标准**：
- 功效 = 0.80（80% 的机会检测到效果）
- α = 0.05（5% I 类错误率）

**四个相互关联的参数**（给定 3 个，可以求解第 4 个）：
1. 样本量（n）
2. 效应大小（d、f等）
3. 显着性水平（α）
4. 功率（1 - β）

---

### 先验功效分析（规划）

**目的**：研究前确定所需的样本量

**步骤**：
1. 指定预期效果大小（来自文献、试点数据或最小有意义的效果）
2. 设置α水平（通常为0.05）
3. 设置所需功率（通常为 0.80）
4. 计算所需n

**Python 实现**：
```python
from statsmodels.stats.power import (
    tt_ind_solve_power,
    zt_ind_solve_power,
    FTestAnovaPower,
    NormalIndPower
)

# T-test power analysis
n_required = tt_ind_solve_power(
    effect_size=0.5,  # Cohen's d
    alpha=0.05,
    power=0.80,
    ratio=1.0,  # Equal group sizes
    alternative='two-sided'
)

# ANOVA power analysis
anova_power = FTestAnovaPower()
n_per_group = anova_power.solve_power(
    effect_size=0.25,  # Cohen's f
    ngroups=3,
    alpha=0.05,
    power=0.80
)

# Correlation power analysis
from pingouin import power_corr
n_required = power_corr(r=0.30, power=0.80, alpha=0.05)
```

---

### 事后功效分析（研究后）

**⚠️ 注意**：事后权力存在争议，通常不推荐

**为什么会出现问题**：
- 观察功效是 p 值的直接函数
- 如果 p > 0.05，功率始终较低
- 不提供除 p 值之外的其他信息
- 可能会产生误导

**何时可以接受**：
- 未来研究的学习计划
- 使用多项研究（不仅仅是您自己的）的效应量
- 明确的目标是复制的样本量

**更好的选择**：
- 报告效应大小的置信区间
- 进行敏感性分析
- 报告最小可检测效应大小

---

### 敏感性分析

**目的**：确定给定研究参数的最小可检测效应大小

**何时使用**：学习完成后，了解学习的能力

**Python 实现**：
```python
# What effect size could we detect with n=50 per group?
detectable_effect = tt_ind_solve_power(
    effect_size=None,  # Solve for this
    nobs1=50,
    alpha=0.05,
    power=0.80,
    ratio=1.0,
    alternative='two-sided'
)

print(f"With n=50 per group, we could detect d ≥ {detectable_effect:.2f}")
```

---

## 报告效果大小

### APA 风格指南

**T 检验示例**：
> “A 组 (M = 75.2，SD = 8.5) 得分显着高于 B 组 (M = 68.3，SD = 9.2)，t(98) = 3.82，p < .001，d = 0.77，95% CI [0.36, 1.18]。”

**方差分析示例**：
> “治疗条件对测试分数有显着的主效应，F(2, 87) = 8.45，p < .001，η²p = .16。使用 Tukey 的 HSD 进行事后比较显示......”

**相关性示例**：
> “学习时间和考试成绩之间存在中度正相关，r(148) = .42，p < .001，95% CI [.27, .55]。”

**回归示例**：
> “回归模型显着预测了考试成绩，F(3, 146) = 45.2，p < .001，R² = .48。学习时间 (β = .52，p < .001) 和之前的 GPA (β = .31，p < .001) 是显着的预测因子。”

**贝叶斯示例**：
> “贝叶斯独立样本 t 检验为组间差异提供了强有力的证据，BF₁₀ = 23.5，表明 H₁ 下的数据可能性比 H₀ 下的可能性高 23.5 倍。”

---

## 效应大小陷阱

1. **不要只依赖基准**：背景很重要；微小的影响也可能有意义
2. **报告置信区间**：置信区间显示效应大小估计的精确度
3. **区分统计意义与实际意义**：大n可以使微不足道的影响变得“显着”
4. **考虑成本效益**：如果干预成本低廉，即使很小的影响也可能是有价值的
5. **多种结果**：效果大小因结果而异；报告全部
6. **不要挑选**：报告所有计划分析的效果
7. **发表偏差**：发表的效果常常被高估

---

## 快速参考表

|分析|效应大小|小|中等|大|
|----------|-------------|--------|--------|--------|
| T 检验 |科恩的 d | 0.20 | 0.20 0.50 | 0.50 0.80 |
|方差分析 | η², ω² | 0.01 | 0.01 0.06 | 0.06 0.14 | 0.14
|方差分析 |科恩的 f | 0.10 | 0.10 0.25 | 0.25 0.40 | 0.40
|相关性| r, ρ | 0.10 | 0.10 0.30 | 0.30 0.50 | 0.50
|回归 | R²| 0.02 | 0.02 0.13 | 0.13 0.26 | 0.26
|回归 | f² | 0.02 | 0.02 0.15 | 0.15 0.35 | 0.35
|卡方|克拉梅尔的 V | 0.07 | 0.07 0.21 | 0.21 0.35 | 0.35
|卡方 (2×2) | φ | 0.10 | 0.10 0.30 | 0.30 0.50 | 0.50

---

## 资源

- 科恩，J. (1988)。 *行为科学的统计功效分析*（第二版）
- 拉肯斯，D. (2013)。计算和报告效应大小
- 埃利斯，P.D. (2010)。 *效果大小的基本指南*