<!-- 此文件由机器翻译自 ensemble-models.md -->

# 用于生存分析的集成模型

## 随机生存森林

### 概述

随机生存森林将随机森林算法扩展到使用审查数据的生存分析。他们在引导样本和聚合预测的基础上构建多个决策树。

### 他们是如何工作的

1. **引导采样**：每棵树都建立在训练数据的不同引导样本上
2. **特征随机性**：在每个节点，仅考虑特征的随机子集进行分割
3. **生存函数估计**：在终端节点，Kaplan-Meier 和 Nelson-Aalen 估计器计算生存函数
4. **集成聚合**：最终预测所有树的平均生存函数

### 何时使用

- 特征和生存之间复杂的非线性关系
- 无需对功能形式进行假设
- 希望通过最少的调整获得可靠的预测
- 需要估计特征重要性
- 有足够的样本量（通常 n > 100）

### 关键参数

- `n_estimators`：树的数量（默认值：100）
  - 更多的树=更稳定的预测，但速度更慢
  - 典型范围：100-1000

- `max_depth`：树的最大深度
  - 控制树的复杂性
  - None = 节点扩展至 pure 或 min_samples_split

- `min_samples_split`：分割节点的最小样本（默认值：6）
  - 更大的值=更多的正则化

- `min_samples_leaf`：叶节点的最小样本（默认值：3）
  - 防止过度拟合小群体

- `max_features`：每次分割时要考虑的特征数量
  - 'sqrt': sqrt(n_features) - 良好的默认值
  - 'log2': log2(n_features)
  - 无：所有功能

- `n_jobs`：并行作业数（-1 使用所有处理器）

### 用法示例

```python
from sksurv.ensemble import RandomSurvivalForest
from sksurv.datasets import load_breast_cancer

# Load data
X, y = load_breast_cancer()

# Fit Random Survival Forest
rsf = RandomSurvivalForest(n_estimators=1000,
                           min_samples_split=10,
                           min_samples_leaf=15,
                           max_features="sqrt",
                           n_jobs=-1,
                           random_state=42)
rsf.fit(X, y)

# Predict risk scores
risk_scores = rsf.predict(X)

# Predict survival functions
surv_funcs = rsf.predict_survival_function(X)

# Predict cumulative hazard functions
chf_funcs = rsf.predict_cumulative_hazard_function(X)
```

### 功能重要性

**重要**：基于分裂杂质的内置特征重要性对于生存数据来说并不可靠。请改用基于排列的特征重要性。

<<<代码块_1>>>

## 梯度提升生存分析

### 概述

梯度提升通过顺序添加弱学习器来纠正先前学习器的错误来构建集成。模型为： **f(x) = Σ β_m g(x; θ_m)**

### 模型类型

#### GradientBoosting 生存分析

使用回归树作为基础学习器。可以捕捉复杂的非线性关系。

**何时使用：**
- 需要对复杂的非线性关系进行建模
- 想要高预测性能
- 有足够的数据以避免过度拟合
- 可以仔细调整超参数

#### ComponentwiseGradientBoosting 生存分析

使用分量最小二乘作为基础学习器。通过自动特征选择生成线性模型。

**何时使用：**
- 想要可解释的线性模型
- 需要自动特征选择（如Lasso）
- 拥有高维数据
- 更喜欢稀疏模型

### 损失函数

#### Cox 的部分似然（默认）

保持比例风险框架，但用加法系综模型替换线性模型。

**适合：**
- 标准生存分析设置
- 当比例风险合理时
- 大多数用例

#### 加速故障时间 (AFT)

假设特征以恒定因子加速或减慢生存时间。损失函数：**(1/n) Σ ω_i (log y_i - f(x_i))²**

**适合：**
- AFT 框架优先于比例风险
- 想要直接模拟时间
- 需要解释对生存时间的影响

### 正则化策略

防止过度拟合的三种主要技术：

1. **学习率** (`learning_rate < 1`)
   - 减少每个基础学习者的贡献
   - 较小的值需要更多的迭代，但需要更好的泛化
   - 典型范围：0.01 - 0.1

2. **退出** (`dropout_rate > 0`)
   - 在训练期间随机删除以前的学习者
   - 迫使学习者变得更加坚强
   - 典型范围：0.01 - 0.2

3. **二次采样** (`subsample < 1`)
   - 每次迭代使用随机数据子集
   - 增加随机性并减少过度拟合
   - 典型范围：0.5 - 0.9

**建议**：将较小的学习率与提前停止相结合以获得最佳性能。

### 关键参数

- `loss`：损失函数（'coxph' 或 'ipcwls'）
- `learning_rate`：缩小每棵树的贡献（默认值：0.1）
- `n_estimators`：提升迭代次数（默认值：100）
- `subsample`：每次迭代的样本分数（默认值：1.0）
- `dropout_rate`：学习者的辍学率（默认值：0.0）
- `max_depth`：树的最大深度（默认值：3）
- `min_samples_split`：分割节点的最小样本（默认值：2）
- `min_samples_leaf`：叶的最小样本（默认值：1）
- `max_features`：每次分割时要考虑的功能

### 用法示例

<<<代码块_2>>>

### 提前停止

使用验证集来防止过度拟合：

<<<代码块_3>>>

### 超参数调整

<<<代码块_4>>>

## ComponentwiseGradientBoosting 生存分析

### 概述

使用分量最小二乘法，生成稀疏线性模型，并具有类似于 Lasso 的自动特征选择。

### 何时使用

- 想要可解释的线性模型
- 需要自动特征选择
- 拥有具有许多不相关特征的高维数据
- 更喜欢基于系数的解释

### 用法示例

<<<代码块_5>>>

## 额外生存树

极其随机的生存树 - 类似于随机生存森林，但在分割选择中具有额外的随机性。

### 何时使用

- 想要比随机生存森林更多的正则化
- 数据有限
- 需要更快的训练

### 主要区别

它不是为选定的特征寻找最佳分割，而是随机选择分割点，从而为整体增加更多多样性。

<<<代码块_6>>>

## 型号对比

|型号|复杂性 |可解释性|性能|速度|
|--------|---------|------------------|--------------|--------|
|随机生存森林 |中等|低|高|中等|
|梯度提升生存分析 |高|低|最高|慢|
| ComponentwiseGradientBoosting 生存分析 |低|高|中等|快|
|额外生存树 |中等|低|中高|快|

**一般建议：**
- **最佳整体性能**：带有调整的 GradientBoostingSurvivalAnalysis
- **最佳平衡**：RandomSurvivalForest
- **最佳可解释性**：ComponentwiseGradientBoostingSurvivalAnalysis
- **最快的训练**：ExtraSurvivalTrees