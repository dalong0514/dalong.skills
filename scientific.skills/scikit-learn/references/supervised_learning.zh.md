<!-- 此文件由机器翻译自 supervised_learning.md -->

# 监督学习参考

## 概述

监督学习算法从标记的训练数据中学习，以对新数据进行预测。 Scikit-learn 为分类和回归任务提供了全面的实现。

## 线性模型

### 回归

**线性回归 (`sklearn.linear_model.LinearRegression`)**
- 普通最小二乘回归
- 快速、可解释、无超参数
- 在以下情况下使用：线性关系、可解释性很重要
- 示例：
```python
from sklearn.linear_model import LinearRegression

model = LinearRegression()
model.fit(X_train, y_train)
predictions = model.predict(X_test)
```

**岭回归 (`sklearn.linear_model.Ridge`)**
- L2正则化防止过度拟合
- 关键参数：`alpha`（正则化强度，默认=1.0）
- 使用时：存在多重共线性，需要正则化
- 示例：
<<<代码块_1>>>

**套索 (`sklearn.linear_model.Lasso`)**
- 带特征选择的 L1 正则化
- 关键参数：`alpha`（正则化强度）
- 使用时：想要稀疏模型、特征选择
- 可以将一些系数减少到恰好为零
- 示例：
<<<代码块_2>>>

**ElasticNet (`sklearn.linear_model.ElasticNet`)**
- 结合 L1 和 L2 正则化
- 关键参数：`alpha`、`l1_ratio`（0=山脊，1=套索）
- 使用时：需要特征选择和正则化
- 示例：
<<<代码块_3>>>

### 分类

**逻辑回归 (`sklearn.linear_model.LogisticRegression`)**
- 二元和多类分类
- 关键参数：`C`（逆正则化）、`penalty`（'l1'、'l2'、'elasticnet'）
- 返回概率估计
- 使用时：需要概率预测、可解释性
- 示例：
<<<代码块_4>>>

**随机梯度下降 (SGD)**
- `SGDClassifier`、`SGDRegressor`
- 高效的大规模学习
- 关键参数：`loss`、`penalty`、`alpha`、`learning_rate`
- 使用时机：非常大的数据集（>10^4 样本）
- 示例：
<<<代码块_5>>>

## 支持向量机

**SVC (`sklearn.svm.SVC`)**
- 使用核方法进行分类
- 关键参数：`C`、`kernel`（'线性'、'rbf'、'聚'）、`gamma`
- 使用场合：中小型数据集、复杂的决策边界
- 注意：不能很好地扩展到大型数据集
- 示例：
<<<代码块_6>>>

**SVR (`sklearn.svm.SVR`)**
- 使用核方法进行回归
- 与SVC类似的参数
- 附加参数：`epsilon`（管宽度）
- 示例：
```python
from sklearn.svm import SVR

model = SVR(kernel='rbf', C=1.0, epsilon=0.1)
model.fit(X_train, y_train)
```

## 决策树

**决策树分类器/决策树回归器**
- 非参数模型学习决策规则
- 关键参数：
  - `max_depth`：最大树深度（防止过度拟合）
  - `min_samples_split`：分割节点的最小样本
  - `min_samples_leaf`：叶子中的最小样本
  - `criterion`：用于分类的“基尼”、“熵”； 'squared_error', 'absolute_error' 用于回归
- 何时使用：需要可解释的模型、非线性关系、混合特征类型
- 容易过度拟合 - 使用集成或修剪
- 示例：
```python
from sklearn.tree import DecisionTreeClassifier

model = DecisionTreeClassifier(
    max_depth=5,
    min_samples_split=20,
    min_samples_leaf=10,
    criterion='gini'
)
model.fit(X_train, y_train)

# Visualize the tree
from sklearn.tree import plot_tree
plot_tree(model, feature_names=feature_names, class_names=class_names)
```

## 集成方法

### 随机森林

**随机森林分类器/随机森林回归器**
- 决策树与装袋的集成
- 关键参数：
  - `n_estimators`：树的数量（默认=100）
  - `max_depth`：最大树深度
  - `max_features`：分割时要考虑的功能（'sqrt'、'log2' 或 int）
  - `min_samples_split`、`min_samples_leaf`：控制树的生长
- 使用场合：需要高精度，可以承受计算
- 提供功能重要性
- 示例：
```python
from sklearn.ensemble import RandomForestClassifier

model = RandomForestClassifier(
    n_estimators=100,
    max_depth=10,
    max_features='sqrt',
    n_jobs=-1  # Use all CPU cores
)
model.fit(X_train, y_train)

# Feature importance
importances = model.feature_importances_
```

### 梯度提升

**GradientBoostingClassifier / GradientBoostingRegressor**
- 基于残差的顺序集成构建树
- 关键参数：
  - `n_estimators`：升压级数
  - `learning_rate`：缩小每棵树的贡献
  - `max_depth`：单个树的深度（通常为 3-5）
  - `subsample`：用于训练每棵树的样本比例
- 使用时机：需要高精度，可以承受训练时间
- 经常取得最佳表现
- 示例：
```python
from sklearn.ensemble import GradientBoostingClassifier

model = GradientBoostingClassifier(
    n_estimators=100,
    learning_rate=0.1,
    max_depth=3,
    subsample=0.8
)
model.fit(X_train, y_train)
```

**HistGradientBoostingClassifier / HistGradientBoostingRegressor**
- 使用基于直方图的算法实现更快的梯度提升
- 对缺失值和分类特征的本机支持
- 关键参数：与GradientBoosting类似
- 使用时机：大型数据集，需要更快的训练
- 示例：
```python
from sklearn.ensemble import HistGradientBoostingClassifier

model = HistGradientBoostingClassifier(
    max_iter=100,
    learning_rate=0.1,
    max_depth=None,  # No limit by default
    categorical_features='from_dtype'  # Auto-detect categorical
)
model.fit(X_train, y_train)
```

### 其他集成方法

**阿达助推**
- 自适应增强专注于错误分类的样本
- 关键参数：`n_estimators`、`learning_rate`、`estimator`（基本估计器）
- 使用时间：需要简单的增强方法
- 示例：
```python
from sklearn.ensemble import AdaBoostClassifier

model = AdaBoostClassifier(n_estimators=50, learning_rate=1.0)
model.fit(X_train, y_train)
```

**投票分类器/回归器**
- 结合多个模型的预测
- 类型：“硬”（多数票）或“软”（平均概率）
- 使用时间：想要集成不同的模型类型
- 示例：
```python
from sklearn.ensemble import VotingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC

model = VotingClassifier(
    estimators=[
        ('lr', LogisticRegression()),
        ('dt', DecisionTreeClassifier()),
        ('svc', SVC(probability=True))
    ],
    voting='soft'
)
model.fit(X_train, y_train)
```

**堆叠分类器/回归器**
- 根据基本模型的预测训练元模型
- 比投票更复杂
- 关键参数：`final_estimator`（元学习器）
- 示例：
```python
from sklearn.ensemble import StackingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC

model = StackingClassifier(
    estimators=[
        ('dt', DecisionTreeClassifier()),
        ('svc', SVC())
    ],
    final_estimator=LogisticRegression()
)
model.fit(X_train, y_train)
```

## K-最近邻

**KNeighbors分类器 / KNeighbors回归器**
- 基于距离的非参数方法
- 关键参数：
  - `n_neighbors`：邻居数量（默认=5）
  - `weights`：“均匀”或“距离”
  - `metric`：距离度量（'euclidean'、'manhattan' 等）
- 使用时机：数据集较小，需要简单的基线
- 对大型数据集的预测速度较慢
- 示例：
```python
from sklearn.neighbors import KNeighborsClassifier

model = KNeighborsClassifier(n_neighbors=5, weights='distance')
model.fit(X_train, y_train)
```

## 朴素贝叶斯

**高斯NB、多项式NB、伯努利NB**
- 基于贝叶斯定理的概率分类器
- 快速训练和预测
- GaussianNB：连续特征（假设高斯分布）
- MultinomialNB：计数特征（文本分类）
- BernoulliNB：二元特征
- 使用场合：文本分类、快速基线、概率预测
- 示例：
```python
from sklearn.naive_bayes import GaussianNB, MultinomialNB

# For continuous features
model_gaussian = GaussianNB()

# For text/count data
model_multinomial = MultinomialNB(alpha=1.0)  # alpha is smoothing parameter
model_multinomial.fit(X_train, y_train)
```

## 神经网络

**MLP分类器/MLP回归器**
- 多层感知器（前馈神经网络）
- 关键参数：
  - `hidden_layer_sizes`：隐藏层大小的元组，例如 (100, 50)
  - `activation`: 'relu', 'tanh', '逻辑'
  - `solver`: 'adam', 'sgd', 'lbfgs'
  - `alpha`：L2正则化参数
  - `learning_rate`：“常量”、“自适应”
- 使用场合：复杂的非线性模式、大型数据集
- 需要特征缩放
- 示例：
```python
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler

# Scale features first
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)

model = MLPClassifier(
    hidden_layer_sizes=(100, 50),
    activation='relu',
    solver='adam',
    alpha=0.0001,
    max_iter=1000
)
model.fit(X_train_scaled, y_train)
```

## 算法选择指南

### 选择基于：

**数据集大小：**
- 小（<1k 样本）：KNN、SVM、决策树
- 中（1k-100k）：随机森林、梯度提升、线性模型
- 大型 (>100k)：SGD、线性模型、HistGradientBoosting

**可解释性：**
- 高：线性模型、决策树
- 中：随机森林（特征重要性）
- 低：带有 RBF 内核的 SVM、神经网络

**准确度与速度：**
- 快速训练：朴素贝叶斯、线性模型、KNN
- 高精度：梯度提升、随机森林、堆叠
- 快速预测：线性模型、朴素贝叶斯
- 预测速度慢：KNN（在大型数据集上）、SVM

**功能类型：**
- 连续：大多数算法运行良好
- 类别：树、HistGradientBoosting（原生支持）
- 混合：树、梯度提升
- 文本：朴素贝叶斯、具有 TF-IDF 的线性模型

**共同出发点：**
1.Logistic Regression（分类）/Linear Regression（回归）——快速基线
2. 随机森林 - 不错的默认选择
3. 梯度提升 - 优化以获得最佳精度