<!-- 此文件由机器翻译自 explainers.md -->

# SHAP 解释器参考

本文档提供有关所有 SHAP 解释器类、其参数、方法以及何时使用每种类型的全面信息。

## 概述

SHAP 为不同模型类型提供专门的解释器，每种模型都针对特定架构进行了优化。一般`shap.Explainer`类会根据模型类型自动选择合适的算法。

## 核心讲解器类

### shap.Explainer（自动选择器）

**目的**：通过选择最合适的解释器算法，自动使用 Shapley 值来解释任何机器学习模型或 Python 函数。

**构造函数参数**：
- `model`：要解释的模型（函数或模型对象）
- `masker`：用于特征操作的背景数据或掩码对象
- `algorithm`：可选覆盖以强制特定的解释器类型
- `output_names`：模型输出的名称
- `feature_names`：输入特征的名称

**何时使用**：不确定使用哪个解释器时的默认选择；根据模型类型自动选择最佳算法。

### 树解释器

**目的**：使用 Tree SHAP 算法快速准确地计算基于树的集成模型的 SHAP 值。

**构造函数参数**：
- `model`：基于树的模型（XGBoost、LightGBM、CatBoost、PySpark 或 scikit-learn 树）
- `data`：用于特征集成的后台数据集（与tree_path_dependent可选）
- `feature_perturbation`：如何处理相关功能
  - `"interventional"`：需要背景数据；遵循因果推理规则
  - `"tree_path_dependent"`：不需要背景数据；每片叶子使用训练样本
  - `"auto"`：如果提供数据，则默认为干预，否则为tree_path_dependent
- `model_output`：解释什么模型输出
  - `"raw"`：标准模型输出（默认）
  - `"probability"`：概率转换输出
  - `"log_loss"`：损失函数的自然对数
  - 自定义方法名称，例如 `"predict_proba"`
- `feature_names`：可选的功能命名

**支持的型号**：
- XGBoost（xgboost.XGBClassifier，xgboost.XGBRegressor，xgboost.Booster）
- LightGBM（lightgbm.LGBMClassifier、lightgbm.LGBMRegressor、lightgbm.Booster）
- CatBoost（catboost.CatBoostClassifier，catboost.CatBoostRegressor）
- PySpark MLlib 树模型
- scikit-learn（DecisionTreeClassifier、DecisionTreeRegressor、RandomForestClassifier、RandomForestRegressor、ExtraTreesClassifier、ExtraTreesRegressor、GradientBoostingClassifier、GradientBoostingRegressor）

**关键方法**：
- `shap_values(X)`：计算样本的SHAP值；返回数组，其中每行代表特征属性
- `shap_interaction_values(X)`：估计特征对之间的交互效果；提供具有主效应和成对相互作用的矩阵
- `explain_row(row)`：用详细的属性信息解释各个行

**何时使用**：
- 所有基于树的模型的主要选择
- 当需要精确的 SHAP 值（而不是近似值）时
- 当计算速度对于大型数据集很重要时
- 适用于随机森林、梯度增强或 XGBoost 等模型

**示例**：
```python
import shap
import xgboost

# Train model
model = xgboost.XGBClassifier().fit(X_train, y_train)

# Create explainer
explainer = shap.TreeExplainer(model)

# Compute SHAP values
shap_values = explainer.shap_values(X_test)

# Compute interaction values
shap_interaction = explainer.shap_interaction_values(X_test)
```

### 深度解释器

**用途**：使用 DeepLIFT 算法的增强版本近似深度学习模型的 SHAP 值。

**构造函数参数**：
- `model`：依赖于框架的规范
  - **TensorFlow**：（input_tensor，output_tensor）的元组，其中输出是单维的
  - **PyTorch**：`nn.Module` 对象或 `(model, layer)` 元组用于特定于层的解释
- `data`：用于功能集成的后台数据集
  - **TensorFlow**：numpy 数组或 pandas DataFrame
  - **PyTorch**：火炬张量
  - **推荐大小**：100-1000 个样本（不是完整的训练集）以平衡准确性和计算成本
- `session`（仅限 TensorFlow）：可选会话对象；如果没有则自动检测
- `learning_phase_flags`：自定义学习阶段张量，用于在推理过程中处理批量归一化/丢失

**支持的框架**：
- **TensorFlow**：全面支持，包括 Keras 模型
- **PyTorch**：与 nn.Module 架构完全集成

**关键方法**：
- `shap_values(X)`：返回应用于数据 X 的模型的近似 SHA 值
- `explain_row(row)`：解释具有属性值和预期输出的单行
- `save(file)` / `load(file)`：解释器对象的序列化支持
- `supports_model_with_masker(model, masker)`：模型类型的兼容性检查器

**何时使用**：
- 适用于 TensorFlow 或 PyTorch 中的深度神经网络
- 使用卷积神经网络 (CNN) 时
- 对于循环神经网络 (RNN) 和 Transformer
- 当深度学习架构需要特定于模型的解释时

**主要设计特点**：
期望估计的方差大约为 1/√N，其中 N 是背景样本的数量，从而实现准确性与效率的权衡。

**示例**：
<<<代码块_1>>>

### 内核解释器

**目的**：使用带有加权线性回归的内核 SHAP 方法进行与模型无关的 SHAP 值计算。

**构造函数参数**：
- `model`：采用样本矩阵并返回模型输出的函数或模型对象
- `data`：用于模拟缺失特征的背景数据集（numpy 数组、pandas DataFrame 或稀疏矩阵）
- `feature_names`：可选的功能名称列表；如果可用，自动从 DataFrame 列名称派生
- `link`：特征重要性与模型输出之间的连接函数
  - `"identity"`：直接关系（默认）
  - `"logit"`：对于概率输出

**关键方法**：
- `shap_values(X, **kwargs)`：计算样本预测的 SHA 值
  - `nsamples`：每个预测的评估计数（“自动”或整数）；较高的值会减少方差
  - `l1_reg`：特征选择正则化（“num_features(int)”、“aic”、“bic”或 float）
  - 返回数组，其中每行的总和等于模型输出与期望值之间的差值
- `explain_row(row)`：用归因值和期望值解释各个预测
- `save(file)` / `load(file)`：保留并恢复解释器对象

**何时使用**：
- 对于没有专门解释器的黑盒模型
- 使用自定义预测函数时
- 适用于任何模型类型（神经网络、SVM、集成方法等）
- 当需要与模型无关的解释时
- **注**：比专业讲解员慢；仅当不存在专门选项时使用

**示例**：
<<<代码块_2>>>

### 线性解释器

**目的**：用于解释特征相关性的线性模型的专门解释器。

**构造函数参数**：
- `model`：线性模型或（系数、截距）元组
- `masker`：特征相关性的背景数据
- `feature_perturbation`：如何处理特征相关性
  - `"interventional"`：假设功能独立
  - `"correlation_dependent"`：考虑特征相关性

**支持的型号**：
- scikit-learn 线性模型（LinearRegression、LogisticRegression、Ridge、Lasso、ElasticNet）
- 具有系数和截距的自定义线性模型

**何时使用**：
- 对于线性回归和逻辑回归模型
- 当特征相关性对于解释准确性很重要时
- 当需要极快的解释时
- 适用于 GLM 和其他线性模型类型

**示例**：
<<<代码块_3>>>

### 梯度解释器

**目的**：使用预期梯度来近似神经网络的 SHAP 值。

**构造函数参数**：
- `model`：深度学习模型（TensorFlow 或 PyTorch）
- `data`：用于集成的背景示例
- `batch_size`：梯度计算的批量大小
- `local_smoothing`：为平滑而添加的噪声量（默认为 0）

**何时使用**：
- 作为神经网络 DeepExplainer 的替代品
- 当首选基于梯度的解释时
- 对于可获得梯度信息的可微模型

**示例**：
<<<代码块_4>>>

### 排列解释器

**目的**：通过迭代输入的排列来近似 Shapley 值。

**构造函数参数**：
- `model`：预测函数
- `masker`：背景数据或屏蔽对象
- `max_evals`：每个样本的最大模型评估数

**何时使用**：
- 当需要精确的 Shapley 值但没有专门的解释器时
- 对于排列易于处理的小功能集
- 作为 KernelExplainer 更准确的替代品（但速度较慢）

**示例**：
<<<代码块_5>>>

## 解释器选择指南

**选择解释者的决策树**：
1. **你的模型是基于树的吗？**（XGBoost、LightGBM、CatBoost、随机森林等）
   - 是→使用`TreeExplainer`（快速且准确）
   - 否 → 继续步骤 2

2. **您的模型是深度神经网络吗？**（TensorFlow、PyTorch、Keras）
   - 是 → 使用 `DeepExplainer` 或 `GradientExplainer`
   - 否 → 继续步骤 3

3. **您的模型是线性的吗？**（线性/逻辑回归，GLM）
   - 是→使用`LinearExplainer`（非常快）
   - 否 → 继续步骤 4

4. **您需要与模型无关的解释吗？**
   - 是→使用`KernelExplainer`（速度较慢，但适用于任何模型）
   - 如果计算预算允许并且需要高精度 → 使用 `PermutationExplainer`

5. **不确定或想要自动选择？**
   - 使用`shap.Explainer`（自动选择最佳算法）

## 解释器中的通用参数

**背景数据/掩码**：
- 目的：代表建立基线期望的“典型”输入
- 尺寸建议：50-1000 个样本（复杂模型更多）
- 选择：从训练数据或 kmeans 选择的代表中随机抽取样本

**功能名称**：
- 自动从pandas DataFrames中提取
- 可以为numpy数组手动指定
- 对于情节的可解释性很重要

**型号输出规格**：
- 原始模型输出与转换输出（概率、对数赔率）
- 对于正确解释 SHAP 值至关重要
- 示例：对于 XGBoost 分类器，SHAP 解释了逻辑转换之前的边际输出（对数赔率）

## 性能考虑因素

**速度排名**（最快到最慢）：
1. `LinearExplainer` - 几乎瞬时
2. `TreeExplainer` - 非常快，可扩展性良好
3. `DeepExplainer` - 神经网络速度快
4. `GradientExplainer` - 神经网络速度快
5. `KernelExplainer` - 慢，仅在必要时使用
6. `PermutationExplainer` - 非常慢，但对于小特征集来说最准确

**内存注意事项**：
- `TreeExplainer`：低内存开销
- `DeepExplainer`：内存与背景样本大小成正比
- `KernelExplainer`：对于大型背景数据集可能会占用大量内存
- 对于大型数据集：使用批处理或样本子集

## 解释器输出：解释对象

所有解释器都会返回 `shap.Explanation` 对象，其中包含：
- `values`：SHAP 值（numpy 数组）
- `base_values`：预期模型输出（基线）
- `data`：原始特征值
- `feature_names`：功能名称

解释对象支持：
- 切片：第一个样本的 `explanation[0]`
- 数组运算：与numpy运算兼容
- 直接绘图：可以传递给绘图函数