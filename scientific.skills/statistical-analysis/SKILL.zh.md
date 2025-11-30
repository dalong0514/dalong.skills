<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：统计分析
描述：“统计分析工具包。假设检验（t 检验、方差分析、卡方）、回归、相关性、贝叶斯统计、功效分析、假设检查、APA 报告，用于学术研究。”
---

# 统计分析

## 概述

统计分析是检验假设和量化关系的系统过程。通过假设检查和 APA 报告进行假设检验（t 检验、方差分析、卡方）、回归、相关性和贝叶斯分析。将这项技能应用于学术研究。

## 何时使用此技能

该技能应该在以下情况下使用：
- 进行统计假设检验（t 检验、方差分析、卡方）
- 执行回归或相关分析
- 运行贝叶斯统计分析
- 检查统计假设和诊断
- 计算效应大小并进行功效分析
- 以APA格式报告统计结果
- 分析实验或观察数据以进行研究

---

## 核心能力

### 1. 测试选择和计划
- 根据研究问题和数据特征选择适当的统计检验
- 进行先验功效分析以确定所需的样本量
- 计划分析策略，包括多重比较校正

### 2. 假设检查
- 在运行测试之前自动验证所有相关假设
- 提供诊断可视化（Q-Q 图、残差图、箱线图）
- 当假设被违反时建议采取补救措施

### 3. 统计测试
- 假设检验：t 检验、方差分析、卡方、非参数替代
- 回归：线性、多重、逻辑、带诊断
- 相关性：Pearson、Spearman，具有置信区间
- 贝叶斯替代方案：贝叶斯 t 检验、方差分析、贝叶斯因子回归

### 4. 效应大小和解释
- 计算并解释所有分析的适当效应量
- 提供效果估计的置信区间
- 区分统计意义和实际意义

### 5. 专业报告
- 生成APA风格的统计报告
- 创建可供出版的图表
- 提供完整的解释以及所有必需的统计数据

---

## 工作流决策树

使用此决策树来确定您的分析路径：

```
START
│
├─ Need to SELECT a statistical test?
│  └─ YES → See "Test Selection Guide"
│  └─ NO → Continue
│
├─ Ready to check ASSUMPTIONS?
│  └─ YES → See "Assumption Checking"
│  └─ NO → Continue
│
├─ Ready to run ANALYSIS?
│  └─ YES → See "Running Statistical Tests"
│  └─ NO → Continue
│
└─ Need to REPORT results?
   └─ YES → See "Reporting Results"
```

---

## 测试选择指南

### 快速参考：选择正确的测试

使用 `references/test_selection_guide.md` 进行全面指导。快速参考：

**比较两组：**
- 独立、连续、正态 → 独立 t 检验
- 独立、连续、非正态→Mann-Whitney U 检验
- 配对、连续、正态 → 配对 t 检验
- 配对、连续、非正态 → Wilcoxon 符号秩检验
- 二元结果→卡方或费舍尔精确检验

**比较 3 个以上组：**
- 独立、连续、正态→单向方差分析
- 独立、连续、非正态 → Kruskal-Wallis 检验
- 配对、连续、正态 → 重复测量方差分析
- 配对、连续、非正态→弗里德曼检验

**关系：**
- 两个连续变量 → Pearson（正态）或 Spearman 相关（非正态）
- 具有预测变量的连续结果 → 线性回归
- 具有预测变量的二元结果 → 逻辑回归

**贝叶斯替代方案：**
所有测试都有贝叶斯版本，提供：
- 关于假设的直接概率陈述
- 贝叶斯因素量化证据
- 支持原假设的能力
- 请参阅`references/bayesian_statistics.md`

---

## 假设检验

### 系统假设验证

**在解释测试结果之前始终检查假设。**

使用提供的 `scripts/assumption_checks.py` 模块进行自动检查：

<<<代码块_1>>>

这执行：
1. **异常值检测**（IQR 和 z 分数方法）
2. **正态性检验**（Shapiro-Wilk 检验 + Q-Q 图）
3. **方差齐性**（Levene 检验 + 箱线图）
4. **解释和建议**

### 个人假设检查

对于有针对性的检查，请使用单独的函数：

<<<代码块_2>>>

### 违反假设时该怎么办

**违反常态：**
- 轻度违规 + 每组 n > 30 → 继续进行参数测试（稳健）
- 中度违规 → 使用非参数替代方案
- 严重违规→转换数据或使用非参数检验

**违反方差齐性：**
- 对于 t 检验 → 使用 Welch 的 t 检验
- 对于方差分析 → 使用 Welch 方差分析或 Brown-Forsythe 方差分析
- 对于回归→使用稳健的标准误差或加权最小二乘法
**违反线性（回归）：**
- 添加多项式项
- 变换变量
- 使用非线性模型或GAM

请参阅 `references/assumptions_and_diagnostics.md` 以获取全面指导。

---

## 运行统计测试

### Python 库

用于统计分析的主要库：
- **scipy.stats**：核心统计测试
- **统计模型**：高级回归和诊断
- **pingouin**：用户友好的统计测试与效果大小
- **pymc**：贝叶斯统计建模
- **arviz**：贝叶斯可视化和诊断

### 示例分析

#### 具有完整报告的 T 测试

<<<代码块_3>>>

#### 方差分析与事后检验

<<<代码块_4>>>

#### 带诊断的线性回归

<<<代码块_5>>>

#### 贝叶斯 T 检验

<<<代码块_6>>>

---

## 效果大小

### 始终计算效果大小

**效应大小量化了幅度，而 p 值仅表明效应的存在。**

请参阅 `references/effect_sizes_and_power.md` 以获取全面指导。

### 快速参考：常见效果大小

|测试|效应大小|小|中等|大|
|------|-------------|--------|--------|--------|
| T 检验 |科恩的 d | 0.20 | 0.20 0.50 | 0.50 0.80 |
|方差分析 | η²_p | 0.01 | 0.01 0.06 | 0.06 0.14 | 0.14
|相关性| r | 0.10 | 0.10 0.30 | 0.30 0.50 | 0.50
|回归 | R² | 0.02 | 0.02 0.13 | 0.13 0.26 | 0.26
|卡方|克拉梅尔的 V | 0.07 | 0.07 0.21 | 0.21 0.35 | 0.35

**重要**：基准是指导方针。背景很重要！

### 计算效果大小

大多数效果大小由 pinouin 自动计算：

```python
# T-test returns Cohen's d
result = pg.ttest(x, y)
d = result['cohen-d'].values[0]

# ANOVA returns partial eta-squared
aov = pg.anova(dv='score', between='group', data=df)
eta_p2 = aov['np2'].values[0]

# Correlation: r is already an effect size
corr = pg.corr(x, y)
r = corr['r'].values[0]
```

### 效应大小的置信区间

始终报告 CI 以显示精度：

```python
from pingouin import compute_effsize_from_t

# For t-test
d, ci = compute_effsize_from_t(
    t_statistic,
    nx=len(group1),
    ny=len(group2),
    eftype='cohen'
)
print(f"d = {d:.2f}, 95% CI [{ci[0]:.2f}, {ci[1]:.2f}]")
```

---

## 功率分析

### 先验功效分析（研究计划）

在数据收集之前确定所需的样本量：

```python
from statsmodels.stats.power import (
    tt_ind_solve_power,
    FTestAnovaPower
)

# T-test: What n is needed to detect d = 0.5?
n_required = tt_ind_solve_power(
    effect_size=0.5,
    alpha=0.05,
    power=0.80,
    ratio=1.0,
    alternative='two-sided'
)
print(f"Required n per group: {n_required:.0f}")

# ANOVA: What n is needed to detect f = 0.25?
anova_power = FTestAnovaPower()
n_per_group = anova_power.solve_power(
    effect_size=0.25,
    ngroups=3,
    alpha=0.05,
    power=0.80
)
print(f"Required n per group: {n_per_group:.0f}")
```

### 敏感性分析（研究后）

确定您可以检测到的效应大小：

```python
# With n=50 per group, what effect could we detect?
detectable_d = tt_ind_solve_power(
    effect_size=None,  # Solve for this
    nobs1=50,
    alpha=0.05,
    power=0.80,
    ratio=1.0,
    alternative='two-sided'
)
print(f"Study could detect d ≥ {detectable_d:.2f}")
```

**注意**：一般不建议进行事后功效分析（研究后计算功效）。请改用敏感性分析。

请参阅 `references/effect_sizes_and_power.md` 了解详细指导。

---

## 报告结果

### APA 风格的统计报告

请遵循 `references/reporting_standards.md` 中的准则。

### 基本报告要素

1. **描述性统计**：所有组/变量的 M、SD、n
2. **检验统计数据**：检验名称、统计量、df、精确 p 值
3. **效果大小**：具有置信区间
4. **假设检查**：进行了哪些测试、结果、采取的行动
5. **所有计划分析**：包括非重大发现

### 报告模板示例

#### 独立 T 检验

```
Group A (n = 48, M = 75.2, SD = 8.5) scored significantly higher than
Group B (n = 52, M = 68.3, SD = 9.2), t(98) = 3.82, p < .001, d = 0.77,
95% CI [0.36, 1.18], two-tailed. Assumptions of normality (Shapiro-Wilk:
Group A W = 0.97, p = .18; Group B W = 0.96, p = .12) and homogeneity
of variance (Levene's F(1, 98) = 1.23, p = .27) were satisfied.
```

#### 单因素方差分析

```
A one-way ANOVA revealed a significant main effect of treatment condition
on test scores, F(2, 147) = 8.45, p < .001, η²_p = .10. Post hoc
comparisons using Tukey's HSD indicated that Condition A (M = 78.2,
SD = 7.3) scored significantly higher than Condition B (M = 71.5,
SD = 8.1, p = .002, d = 0.87) and Condition C (M = 70.1, SD = 7.9,
p < .001, d = 1.07). Conditions B and C did not differ significantly
(p = .52, d = 0.18).
```

#### 多元回归

```
Multiple linear regression was conducted to predict exam scores from
study hours, prior GPA, and attendance. The overall model was significant,
F(3, 146) = 45.2, p < .001, R² = .48, adjusted R² = .47. Study hours
(B = 1.80, SE = 0.31, β = .35, t = 5.78, p < .001, 95% CI [1.18, 2.42])
and prior GPA (B = 8.52, SE = 1.95, β = .28, t = 4.37, p < .001,
95% CI [4.66, 12.38]) were significant predictors, while attendance was
not (B = 0.15, SE = 0.12, β = .08, t = 1.25, p = .21, 95% CI [-0.09, 0.39]).
Multicollinearity was not a concern (all VIF < 1.5).
```

#### 贝叶斯分析

```
A Bayesian independent samples t-test was conducted using weakly
informative priors (Normal(0, 1) for mean difference). The posterior
distribution indicated that Group A scored higher than Group B
(M_diff = 6.8, 95% credible interval [3.2, 10.4]). The Bayes Factor
BF₁₀ = 45.3 provided very strong evidence for a difference between
groups, with a 99.8% posterior probability that Group A's mean exceeded
Group B's mean. Convergence diagnostics were satisfactory (all R̂ < 1.01,
ESS > 1000).
```

---

## 贝叶斯统计

### 何时使用贝叶斯方法

在以下情况下考虑贝叶斯方法：
- 您有要合并的先前信息
- 您想要关于假设的直接概率陈述
- 样本量较小或计划连续数据收集
- 您需要量化原假设的证据
- 模型复杂（分层、缺失数据）

请参阅 `references/bayesian_statistics.md` 以获取以下方面的综合指南：
- 贝叶斯定理及其解释
- 先前的规范（信息性、弱信息性、非信息性）
- 使用贝叶斯因子进行贝叶斯假设检验
- 可信区间与置信区间
- 贝叶斯 t 检验、方差分析、回归和分层模型
- 模型收敛性检查和后验预测检查

### 主要优势

1. **直观解释**：“根据数据，参数有 95% 的概率位于该区间内”
2. **无效的证据**：可以量化无效的支持
3. **灵活**：无 p-hacking 问题；可以在数据到达时对其进行分析
4. **不确定性量化**：完全后验分布

---

## 资源

该技能包括全面的参考资料：

### 参考目录

- **test_selection_guide.md**：用于选择适当统计测试的决策树
- **assections_and_diagnostics.md**：检查和处理假设违规的详细指南
- **effect_sizes_and_power.md**：计算、解释和报告效应大小；进行功率分析
- **bayesian_statistics.md**：贝叶斯分析方法完整指南
- **reporting_standards.md**：带有示例的 APA 风格报告指南

### 脚本目录

- **assump_checks.py**：通过可视化自动进行假设检查
  - `comprehensive_assumption_check()`：完整的工作流程
- `check_normality()`：使用 Q-Q 图进行正态性测试
  - `check_homogeneity_of_variance()`：带有箱线图的 Levene 检验
  - `check_linearity()`：回归线性检查
  - `detect_outliers()`：IQR 和 z 分数异常值检测

---

## 最佳实践

1. **预先注册分析** 尽可能区分验证性和探索性
2. **在解释结果之前始终检查假设**
3. **报告效应大小**以及置信区间
4. **报告所有计划的分析**，包括不显着的结果
5. **区分统计意义和实际意义**
6. **分析前后的数据可视化**
7. **检查回归/方差分析的诊断**（残差图、VIF 等）
8. **进行敏感性分析**以评估稳健性
9. **共享数据和代码**以实现可重复性
10. **对违规、转变和决策保持透明**

---

## 要避免的常见陷阱

1. **P-hacking**：在事情变得重要之前不要测试多种方法
2. **倾听**：不要将探索性发现作为证实性证据
3. **忽略假设**：检查并报告违规情况
4. **混淆显着性和重要性**：p < .05 ≠ 有意义的效果
5. **不报告效应大小**：对于解释至关重要
6. **挑选结果**：报告所有计划的分析
7. **误解 p 值**：它们不是假设成立的概率
8. **多重比较**：在适当的时候纠正家族错误
9. **忽略缺失数据**：了解机制（MCAR、MAR、MNAR）
10. **过度解释非显着性结果**：缺乏证据≠不存在证据

---

## 入门清单

开始统计分析时：

- [ ] 定义研究问题和假设
- [ ] 确定适当的统计测试（使用 test_selection_guide.md）
- [ ] 进行功效分析以确定样本量
- [ ] 加载和检查数据
- [ ] 检查缺失数据和异常值
- [ ] 使用假设_checks.py 验证假设
- [ ] 运行初步分析
- [ ] 使用置信区间计算效应大小
- [ ] 如果需要的话进行事后测试（带修正）
- [ ] 创建可视化
- [ ] 按照reporting_standards.md写入结果
- [ ] 进行敏感性分析
- [ ] 分享数据和代码

---

## 支持和进一步阅读

对于以下问题：
- **测试选择**：参见references/test_selection_guide.md
- **假设**：请参阅references/assurations_and_diagnostics.md
- **效果大小**：参见references/effect_sizes_and_power.md
- **贝叶斯方法**：参见references/bayesian_statistics.md
- **报告**：参见references/reporting_standards.md

**主要教材**：
- 科恩，J. (1988)。 *行为科学的统计功效分析*
- 菲尔德，A.（2013）。 *使用 IBM SPSS Statistics 发现统计数据*
- Gelman, A. 和 Hill, J. (2006)。 *使用回归和多级/分层模型进行数据分析*
- 克鲁施克，J.K. (2014)。 *进行贝叶斯数据分析*

**在线资源**：
- APA 风格指南：https://apastyle.apa.org/
- 统计咨询：交叉验证 (stats.stackexchange.com)