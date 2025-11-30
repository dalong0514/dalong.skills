<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：scikit学习
描述：使用 scikit-learn 使用 Python 进行机器学习。在处理监督学习（分类、回归）、无监督学习（聚类、降维）、模型评估、超参数调整、预处理或构建 ML 管道时使用。提供有关算法、预处理技术、管道和最佳实践的全面参考文档。
---

# Scikit 学习

## 概述

该技能为使用 scikit-learn（用于经典机器学习的行业标准 Python 库）执行机器学习任务提供了全面的指导。使用此技能进行分类、回归、聚类、降维、预处理、模型评估和构建可用于生产的 ML 管道。

## 安装

```bash
# Install scikit-learn using uv
uv uv pip install scikit-learn

# Optional: Install visualization dependencies
uv uv pip install matplotlib seaborn

# Commonly used with
uv uv pip install pandas numpy
```

## 何时使用此技能

在以下情况下使用 scikit-learn 技能：

- 建立分类或回归模型
- 执行聚类或降维
- 机器学习数据的预处理和转换
- 通过交叉验证评估模型性能
- 通过网格或随机搜索调整超参数
- 为生产工作流程创建机器学习管道
- 比较任务的不同算法
- 使用结构化（表格）和文本数据
- 需要可解释的经典机器学习方法

## 快速入门

### 分类示例

<<<代码块_1>>>

### 具有混合数据的完整管道

<<<代码块_2>>>

## 核心能力

### 1.监督学习

用于分类和回归任务的综合算法。

**关键算法：**
- **线性模型**：逻辑回归、线性回归、Ridge、Lasso、ElasticNet
- **基于树**：决策树、随机森林、梯度提升
- **支持向量机**：具有各种内核的 SVC、SVR
- **集成方法**：AdaBoost、投票、堆叠
- **神经网络**：MLPClassifier、MLPRegressor
- **其他**：朴素贝叶斯、K 最近邻

**何时使用：**
- 分类：预测离散类别（垃圾邮件检测、图像分类、欺诈检测）
- 回归：预测连续值（价格预测、需求预测）

**请参阅：** `references/supervised_learning.md` 了解详细的算法文档、参数和使用示例。

### 2.无监督学习

通过聚类和降维发现未标记数据中的模式。

**聚类算法：**
- **基于分区**：K-Means、MiniBatchKMeans
- **基于密度**：DBSCAN、HDBSCAN、OPTICS
- **分层**：凝聚聚类
- **概率**：高斯混合模型
- **其他**：MeanShift、SpectralClustering、BIRCH

**降维：**
- **线性**：PCA、TruncatedSVD、NMF
- **流形学习**：t-SNE、UMAP、Isomap、LLE
- **特征提取**：FastICA、LatentDirichletAllocation

**何时使用：**
- 客户细分、异常检测、数据可视化
- 降低特征维度，探索性数据分析
- 主题建模、图像压缩

**请参阅：** `references/unsupervised_learning.md` 了解详细文档。

### 3.模型评估与选择

用于稳健模型评估、交叉验证和超参数调整的工具。

**交叉验证策略：**
- KFold、StratifiedKFold（分类）
- TimeSeriesSplit（时间数据）
- GroupKFold（分组样本）

**超参数调整：**
- GridSearchCV（穷举搜索）
- RandomizedSearchCV（随机采样）
- HalvingGridSearchCV（连续减半）

**指标：**
- **分类**：准确度、精确度、召回率、F1 分数、ROC AUC、混淆矩阵
- **回归**：MSE、RMSE、MAE、R²、MAPE
- **聚类**：轮廓得分、Calinski-Harabasz、Davies-Bouldin

**何时使用：**
- 客观比较模型性能
- 寻找最佳超参数
- 通过交叉验证防止过度拟合
- 通过学习曲线了解模型行为

**请参阅：** `references/model_evaluation.md` 了解全面的指标和调整策略。

### 4.数据预处理

将原始数据转换为适合机器学习的格式。

**缩放和标准化：**
- StandardScaler（零均值，单位方差）
- MinMaxScaler（有界范围）
- RobustScaler（对异常值具有鲁棒性）
- 标准化器（样本标准化）

**编码分类变量：**
- OneHotEncoder（名义类别）
- OrdinalEncoder（有序类别）
- LabelEncoder（目标编码）

**处理缺失值：**
- SimpleImputer（平均值、中值、最频繁）
- KNNImputer（k 最近邻）
- IterativeImputer（多元插补）

**特征工程：**
- 多项式特征（交互项）
- KBinsDiscretizer（分箱）
- 特征选择（RFE、SelectKBest、SelectFromModel）

**何时使用：**
- 在训练任何需要缩放特征的算法（SVM、KNN、神经网络）之前
- 将分类变量转换为数字格式
- 系统地处理缺失数据
- 为线性模型创建非线性特征

**请参阅：** `references/preprocessing.md` 了解详细的预处理技术。

### 5. 管道和组合

构建可重复、可用于生产的 ML 工作流程。

**关键部件：**
- **管道**：按顺序链接变压器和估计器
- **ColumnTransformer**：对不同的列应用不同的预处理
- **FeatureUnion**：并行组合多个变压器
- **TransformedTargetRegressor**：变换目标变量

**好处：**
- 防止交叉验证中的数据泄漏
- 简化代码并提高可维护性
- 启用联合超参数调整
- 确保训练和预测之间的一致性

**何时使用：**
- 始终使用 Pipelines 进行生产工作流程
- 混合数字和分类特征时（使用 ColumnTransformer）
- 使用预处理步骤执行交叉验证时
- 当超参数调整包括预处理参数时

**请参阅：** `references/pipelines_and_composition.md` 了解全面的管道模式。

## 示例脚本

### 分类管道

运行完整的分类工作流程，包括预处理、模型比较、超参数调整和评估：

<<<代码块_3>>>

该脚本演示了：
- 处理混合数据类型（数字和分类）
- 使用交叉验证进行模型比较
- 使用 GridSearchCV 进行超参数调整
- 多指标综合评价
- 特征重要性分析

### 聚类分析

通过算法比较和可视化进行聚类分析：

<<<代码块_4>>>

该脚本演示了：
- 寻找最佳簇数（肘部法、轮廓分析）
- 比较多种聚类算法（K-Means、DBSCAN、Agglomerative、Gaussian Mixture）
- 在没有基本事实的情况下评估聚类质量
- 使用 PCA 投影可视化结果

## 参考文档

该技能包括用于深入研究特定主题的综合参考文件：

### 快速参考
**文件：** `references/quick_reference.md`
- 常见导入模式和安装说明
- 常见任务的快速工作流程模板
- 算法选择备忘单
- 常见模式和陷阱
- 性能优化技巧

### 监督学习
**文件：** `references/supervised_learning.md`
- 线性模型（回归和分类）
- 支持向量机
- 决策树和集成方法
- K 最近邻、朴素贝叶斯、神经网络
- 算法选择指南

### 无监督学习
**文件：** `references/unsupervised_learning.md`
- 所有带有参数和用例的聚类算法
- 降维技术
- 异常值和新颖性检测
- 高斯混合模型
- 方法选择指南

### 模型评估
**文件：** `references/model_evaluation.md`
- 交叉验证策略
- 超参数调整方法
- 分类、回归和聚类指标
- 学习和验证曲线
- 模型选择的最佳实践

### 预处理
**文件：** `references/preprocessing.md`
- 特征缩放和标准化
- 编码分类变量
- 缺失值插补
- 特征工程技术
- 定制变压器

### 管道和组合
**文件：** `references/pipelines_and_composition.md`
- 管道建设及使用
- 用于混合数据类型的ColumnTransformer
- 用于并行转换的FeatureUnion
- 完整的端到端示例
- 最佳实践

## 常见工作流程

### 构建分类模型

1. **加载和探索数据**
   <<<代码块_5>>>

2. **通过分层分割数据**
   <<<代码块_6>>>

3. **创建预处理管道**
   ```python
   from sklearn.pipeline import Pipeline
   from sklearn.preprocessing import StandardScaler
   from sklearn.compose import ColumnTransformer

   # Handle numeric and categorical features separately
   preprocessor = ColumnTransformer([
       ('num', StandardScaler(), numeric_features),
       ('cat', OneHotEncoder(), categorical_features)
   ])
   ```

4. **建立完整的管道**
   ```python
   model = Pipeline([
       ('preprocessor', preprocessor),
       ('classifier', RandomForestClassifier(random_state=42))
   ])
   ```

5. **调整超参数**
   ```python
   from sklearn.model_selection import GridSearchCV

   param_grid = {
       'classifier__n_estimators': [100, 200],
       'classifier__max_depth': [10, 20, None]
   }

   grid_search = GridSearchCV(model, param_grid, cv=5)
   grid_search.fit(X_train, y_train)
   ```

6. **在测试集上进行评估**
   ```python
   from sklearn.metrics import classification_report

   best_model = grid_search.best_estimator_
   y_pred = best_model.predict(X_test)
   print(classification_report(y_test, y_pred))
   ```

### 执行聚类分析

1. **预处理数据**
   ```python
   from sklearn.preprocessing import StandardScaler

   scaler = StandardScaler()
   X_scaled = scaler.fit_transform(X)
   ```

2. **找到最佳簇数**
   ```python
   from sklearn.cluster import KMeans
   from sklearn.metrics import silhouette_score

   scores = []
   for k in range(2, 11):
       kmeans = KMeans(n_clusters=k, random_state=42)
       labels = kmeans.fit_predict(X_scaled)
       scores.append(silhouette_score(X_scaled, labels))

   optimal_k = range(2, 11)[np.argmax(scores)]
   ```

3. **应用聚类**
   ```python
   model = KMeans(n_clusters=optimal_k, random_state=42)
   labels = model.fit_predict(X_scaled)
   ```

4. **通过降维进行可视化**
   ```python
   from sklearn.decomposition import PCA

   pca = PCA(n_components=2)
   X_2d = pca.fit_transform(X_scaled)

   plt.scatter(X_2d[:, 0], X_2d[:, 1], c=labels, cmap='viridis')
   ```

## 最佳实践

### 始终使用管道
管道防止数据泄露并确保一致性：
```python
# Good: Preprocessing in pipeline
pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('model', LogisticRegression())
])

# Bad: Preprocessing outside (can leak information)
X_scaled = StandardScaler().fit_transform(X)
```

### 仅适合训练数据
永远不要拟合测试数据：
```python
# Good
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)  # Only transform

# Bad
scaler = StandardScaler()
X_all_scaled = scaler.fit_transform(np.vstack([X_train, X_test]))
```

### 使用分层分割进行分类
保留类别分布：
```python
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=42
)
```

### 设置随机状态以实现可重复性
```python
model = RandomForestClassifier(n_estimators=100, random_state=42)
```

### 选择适当的指标
- 平衡数据：准确性、F1 分数
- 不平衡数据：精确率、召回率、ROC AUC、平衡准确率
- 成本敏感：定义自定义记分器

### 在需要时扩展功能
需要特征缩放的算法：
- SVM、KNN、神经网络
- PCA，具有正则化的线性/逻辑回归
- K-Means聚类

不需要缩放的算法：
- 基于树的模型（决策树、随机森林、梯度提升）
- 朴素贝叶斯

## 常见问题故障排除

### 收敛警告
**问题：**模型未收敛
**解决方案：**增加`max_iter`或缩放特征
```python
model = LogisticRegression(max_iter=1000)
```

### 测试集表现不佳
**问题：** 过度拟合
**解决方案：** 使用正则化、交叉验证或更简单的模型
```python
# Add regularization
model = Ridge(alpha=1.0)

# Use cross-validation
scores = cross_val_score(model, X, y, cv=5)
```

### 大型数据集的内存错误
**解决方案：** 使用专为大数据设计的算法
```python
# Use SGD for large datasets
from sklearn.linear_model import SGDClassifier
model = SGDClassifier()

# Or MiniBatchKMeans for clustering
from sklearn.cluster import MiniBatchKMeans
model = MiniBatchKMeans(n_clusters=8, batch_size=100)
```

## 其他资源

- 官方文档：https://scikit-learn.org/stable/
- 用户指南：https://scikit-learn.org/stable/user_guide.html
- API 参考：https://scikit-learn.org/stable/api/index.html
- 示例库：https://scikit-learn.org/stable/auto_examples/index.html