<!-- 此文件由机器翻译自 cox-models.md -->

# Cox 比例风险模型

## 概述

Cox 比例风险模型是将协变量与事件时间相关联的半参数模型。个体 *i* 的危险函数表示为：

**h_i(t) = h_0(t) × exp(β^T x_i)**

其中：
- h_0(t) 是基线危险函数（未指定）
- β 是系数向量
- x_i 是个体 *i* 的协变量向量

关键假设是两个人之间的风险比随着时间的推移保持不变（比例风险）。

## CoxPHS 生存分析

用于生存分析的基本 Cox 比例风险模型。

### 何时使用
- 使用审查数据进行标准生存分析
- 需要可解释的系数（对数风险比）
- 比例风险假设成立
- 数据集的特征相对较少

### 关键参数
- `alpha`：正则化参数（默认：0，无正则化）
- `ties`：处理绑定事件时间的方法（'breslow' 或 'efron'）
- `n_iter`：优化的最大迭代次数

### 用法示例
```python
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.datasets import load_gbsg2

# Load data
X, y = load_gbsg2()

# Fit Cox model
estimator = CoxPHSurvivalAnalysis()
estimator.fit(X, y)

# Get coefficients (log hazard ratios)
coefficients = estimator.coef_

# Predict risk scores
risk_scores = estimator.predict(X)
```

## Coxnet 生存分析

具有用于特征选择和正则化的弹性净惩罚的 Cox 模型。

### 何时使用
- 高维数据（很多特征）
- 需要自动特征选择
- 想要处理多重共线性
- 需要稀疏模型

### 处罚类型
- **山脊 (L2)**：alpha_min_ratio=1.0，l1_ratio=0
  - 缩小所有系数
  - 当所有功能都相关时很好

- **套索 (L1)**：l1_ratio=1.0
  - 执行特征选择（将系数设置为零）
  - 适用于稀疏模型

- **弹性网络**：0 < l1_ratio < 1
  - L1 和 L2 的组合
  - 平衡特征选择和分组

### 关键参数
- `l1_ratio`：L1 和 L2 惩罚之间的平衡（0=Ridge，1=Lasso）
- `alpha_min_ratio`：正则化路径中最小与最大惩罚的比率
- `n_alphas`：正则化路径上的 alpha 数量
- `fit_baseline_model`：是否拟合未惩罚基线模型

### 用法示例
<<<代码块_1>>>

### Alpha 选择的交叉验证
<<<代码块_2>>>

## IPCRridge

加速失效时间模型的审查加权岭回归的逆概率。

### 何时使用
- 优先考虑加速失效时间 (AFT) 框架而不是比例风险
- 需要对特征如何加速/减慢生存时间进行建模
- 高审查率
- 希望通过 Ridge 惩罚进行正则化

### 与 Cox 模型的主要区别
AFT 模型假设特征将生存时间乘以常数因子，而不是乘以危险率。该模型直接预测对数生存时间。

### 用法示例
<<<代码块_3>>>

## 模型比较与选择

### 在模型之间进行选择

**在以下情况下使用 CoxPHSurvivalAnalysis：**
- 少量到中等数量的功能
- 想要可解释的风险比
- 标准生存分析设置

**在以下情况下使用 CoxnetSurvivalAnalysis：**
- 高维数据 (p >> n)
- 需要特征选择
- 想要确定重要的预测因素
- 存在多重共线性

**在以下情况下使用 IPCRidge：**
- AFT框架更合适
- 高审查率
- 想要直接模拟时间而不是危险

### 检查比例风险假设

应使用以下方法验证比例风险假设：
- 舍恩菲尔德残差
- 对数生存图
- 统计测试（可在其他软件包中使用，例如生命线）

如果违反，请考虑：
- 通过违反协变量进行分层
- 时变系数
- 替代模型（AFT、参数模型）

## 解释

### Cox 模型系数
- 正系数：危险增加（生存期缩短）
- 负系数：危险降低（生存时间更长）
- 协变量增加一单位的风险比 = exp(β)
- 示例：β=0.693 → HR=2.0（风险加倍）

### 风险评分
- 较高的风险评分=较高的事件风险=较短的预期生存期
- 风险评分是相对的；使用生存函数进行绝对预测