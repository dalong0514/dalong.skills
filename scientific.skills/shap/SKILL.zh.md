<!-- 此文件由机器翻译自 SKILL.md -->

---
名称： 形状
描述：使用 SHAP（SHapley Additive exPlanations）进行模型可解释性和可解释性。在解释机器学习模型预测、计算特征重要性、生成 SHAP 图（瀑布图、蜂群图、条形图、散点图、力图、热图）、调试模型、分析模型偏差或公平性、比较模型或实现可解释的 AI 时，请使用此技能。适用于基于树的模型（XGBoost、LightGBM、随机森林）、深度学习（TensorFlow、PyTorch）、线性模型和任何黑盒模型。
---

# SHAP（SHapley 添加剂解释）

## 概述

SHAP 是一种使用合作博弈论中的 Shapley 值来解释机器学习模型输出的统一方法。该技能为以下方面提供全面指导：

- 计算任何模型类型的 SHAP 值
- 创建可视化以了解功能重要性
- 调试和验证模型行为
- 分析公平性和偏见
- 在生产中实施可解释的人工智能

SHAP 适用于所有模型类型：基于树的模型（XGBoost、LightGBM、CatBoost、随机森林）、深度学习模型（TensorFlow、PyTorch、Keras）、线性模型和黑盒模型。

## 何时使用此技能

**当用户询问时触发此技能**：
- “解释哪些特征在我的模型中最重要”
- “生成 SHAP 图”（瀑布图、蜂群图、条形图、散点图、力图、热图等）
- “为什么我的模型会做出这样的预测？”
- “计算我的模型的 SHAP 值”
- “使用 SHAP 可视化特征重要性”
- “调试我的模型的行为”或“验证我的模型”
-“检查我的模型是否存在偏见”或“分析公平性”
- “比较不同模型的特征重要性”
-“实施可解释的人工智能”或“为我的模型添加解释”
- “了解功能交互”
- “创建模型解释仪表板”

## 快速入门指南

### 第 1 步：选择正确的讲解员

**决策树**：

1. **基于树的模型？**（XGBoost、LightGBM、CatBoost、随机森林、梯度提升）
   - 使用`shap.TreeExplainer`（快速、准确）

2. **深度神经网络？**（TensorFlow、PyTorch、Keras、CNN、RNN、Transformers）
   - 使用`shap.DeepExplainer`或`shap.GradientExplainer`

3. **线性模型？**（线性/逻辑回归，GLM）
   - 使用`shap.LinearExplainer`（非常快）

4. **还有其他模型吗？**（SVM、自定义函数、黑盒模型）
   - 使用`shap.KernelExplainer`（与模型无关但速度较慢）

5. **不确定？**
   - 使用`shap.Explainer`（自动选择最佳算法）

**有关所有解释器类型的详细信息，请参阅 `references/explainers.md`。**

### 步骤 2：计算 SHAP 值

```python
import shap

# Example with tree-based model (XGBoost)
import xgboost as xgb

# Train model
model = xgb.XGBClassifier().fit(X_train, y_train)

# Create explainer
explainer = shap.TreeExplainer(model)

# Compute SHAP values
shap_values = explainer(X_test)

# The shap_values object contains:
# - values: SHAP values (feature attributions)
# - base_values: Expected model output (baseline)
# - data: Original feature values
```

### 步骤 3：可视化结果

**为了全球理解**（整个数据集）：
<<<代码块_1>>>

**对于个人预测**：
<<<代码块_2>>>

**对于特征关系**：
<<<代码块_3>>>

**有关所有绘图类型的综合指南，请参阅 `references/plots.md`。**

## 核心工作流程

此技能支持多种常见工作流程。选择与当前任务匹配的工作流程。

### 工作流程1：基本模型解释

**目标**：了解驱动模型预测的因素

**步骤**：
1. 训练模型并创建适当的解释器
2. 计算测试集的SHAP值
3. 生成全局重要性图（蜂群或条形图）
4. 检查最重要的特征关系（散点图）
5.解释具体预测（瀑布图）

**示例**：
<<<代码块_4>>>

### 工作流程2：模型调试

**目标**：识别并修复模型问题

**步骤**：
1. 计算SHA值
2. 识别预测错误
3. 解释错误分类的样本
4. 检查意外的特征重要性（数据泄漏）
5. 验证特征关系是否有意义
6. 检查功能交互

**请参阅 `references/workflows.md` 了解详细的调试工作流程。**

### 工作流程 3：特征工程

**目标**：使用 SHAP 见解来改进功能

**步骤**：
1. 计算基线模型的SHA值
2. 识别非线性关系（转换候选对象）
3. 识别特征交互（交互术语的候选）
4. 设计新功能
5. 重新训练并比较 SHAP 值
6. 验证改进

**请参阅 `references/workflows.md` 了解详细的特征工程工作流程。**

### 工作流程 4：模型比较

**目标**：比较多个模型以选择最佳可解释选项

**步骤**：
1. 训练多个模型
2. 计算每个的 SHA 值
3. 比较全局特征重要性
4. 检查特征排名的一致性
5. 分析跨模型的具体预测
6. 根据准确性、可解释性和一致性进行选择

**请参阅 `references/workflows.md` 了解详细的模型比较工作流程。**

### 工作流程 5：公平性和偏见分析

**目标**：检测并分析跨人口群体的模型偏差

**步骤**：
1. 识别受保护的属性（性别、种族、年龄等）
2. 计算 SHAP 值
3. 比较各组特征的重要性
4. 检查受保护属性SHAP重要性
5. 识别代理特征
6. 如果发现偏见，实施缓解策略

**请参阅 `references/workflows.md` 了解详细的公平性分析工作流程。**

### 工作流程 6：生产部署

**目标**：将 SHAP 解释集成到生产系统中

**步骤**：
1. 训练并保存模型
2. 创建并保存解释器
3. 搭建讲解服务
4. 创建带有解释的预测 API 端点
5. 实现缓存和优化
6. 监控解释质量

**有关详细的生产部署工作流程，请参阅 `references/workflows.md`。**

## 关键概念

### 形状值

**定义**：SHAP 值量化每个特征对预测的贡献，以与预期模型输出（基线）的偏差来衡量。

**属性**：
- **相加性**：SHAP 值总和为预测与基线之间的差异
- **公平**：基于博弈论中的 Shapley 值
- **一致性**：如果某个功能变得更加重要，它的 SHAP 值就会增加

**解释**：
- 正 SHA 值 → 特征将预测推高
- 负 SHAP 值 → 特征降低预测值
- 大小 → 特征影响的强度
- SHAP 值的总和 → 相对于基线的总预测变化

**示例**：
<<<代码块_5>>>

### 背景数据/基线

**目的**：代表建立基线期望的“典型”输入

**选择**：
- 从训练数据中随机抽取样本（50-1000 个样本）
- 或者使用kmeans选择代表性样本
- 对于 DeepExplainer/KernelExplainer：100-1000 个样本平衡准确性和速度

**影响**：基线影响 SHAP 值大小，但不影响相对重要性

### 模型输出类型

**关键考虑因素**：了解模型的输出内容

- **原始输出**：用于回归或树边距
- **概率**：用于分类概率
- **对数赔率**：用于逻辑回归（在 sigmoid 之前）

**示例**：XGBoost 分类器默认解释边际输出（对数赔率）。要解释概率，请在 TreeExplainer 中使用 `model_output="probability"`。

## 常见模式

### 模式1：完整的模型分析

<<<代码块_6>>>

### 模式 2：群组比较

```python
# Define cohorts
cohort1_mask = X_test['Group'] == 'A'
cohort2_mask = X_test['Group'] == 'B'

# Compare feature importance
shap.plots.bar({
    "Group A": shap_values[cohort1_mask],
    "Group B": shap_values[cohort2_mask]
})
```

### 模式 3：调试错误

```python
# Find errors
errors = model.predict(X_test) != y_test
error_indices = np.where(errors)[0]

# Explain errors
for idx in error_indices[:5]:
    print(f"Sample {idx}:")
    shap.plots.waterfall(shap_values[idx])

    # Investigate key features
    shap.plots.scatter(shap_values[:, "Suspicious_Feature"])
```

## 性能优化

### 速度考虑因素

**解释速度**（最快到最慢）：
1. `LinearExplainer` - 几乎瞬时
2. `TreeExplainer` - 非常快
3. `DeepExplainer` - 神经网络快速
4. `GradientExplainer` - 神经网络速度快
5. `KernelExplainer` - 慢（仅在必要时使用）
6. `PermutationExplainer` - 非常慢但准确

### 优化策略

**对于大型数据集**：
```python
# Compute SHAP for subset
shap_values = explainer(X_test[:1000])

# Or use batching
batch_size = 100
all_shap_values = []
for i in range(0, len(X_test), batch_size):
    batch_shap = explainer(X_test[i:i+batch_size])
    all_shap_values.append(batch_shap)
```

**对于可视化**：
```python
# Sample subset for plots
shap.plots.beeswarm(shap_values[:1000])

# Adjust transparency for dense plots
shap.plots.scatter(shap_values[:, "Feature"], alpha=0.3)
```

**对于生产**：
```python
# Cache explainer
import joblib
joblib.dump(explainer, 'explainer.pkl')
explainer = joblib.load('explainer.pkl')

# Pre-compute for batch predictions
# Only compute top N features for API responses
```

## 故障排除

### 问题：解释者选择错误
**问题**：将 KernelExplainer 用于树模型（缓慢且不必要）
**解决方案**：对于基于树的模型始终使用 TreeExplainer

### 问题：背景数据不足
**问题**：DeepExplainer/KernelExplainer 的背景样本太少
**解决方案**：使用100-1000个代表性样本

### 问题：令人困惑的单位
**问题**：将对数赔率解释为概率
**解决办法**：检查模型输出类型；了解值是概率、对数赔率还是原始输出

### 问题：绘图不显示
**问题**：Matplotlib 后端问题
**解决办法**：确保后端设置正确；如果需要，请使用 `plt.show()`

### 问题：太多的功能使绘图变得混乱
**问题**：默认 max_display=10 可能太多或太少
**解决方案**：调整`max_display`参数或使用特征聚类

### 问题：计算速度慢
**问题**：计算非常大的数据集的 SHAP
**解决方案**：样本子集，使用批处理，或确保使用专门的解释器（而不是 KernelExplainer）

## 与其他工具集成

### Jupyter 笔记本
- 交互式力图无缝工作
- 内联绘图显示 `show=True` （默认）
- 结合Markdown进行叙述性解释

### MLflow / 实验跟踪
```python
import mlflow

with mlflow.start_run():
    # Train model
    model = train_model(X_train, y_train)

    # Compute SHAP
    explainer = shap.TreeExplainer(model)
    shap_values = explainer(X_test)

    # Log plots
    shap.plots.beeswarm(shap_values, show=False)
    mlflow.log_figure(plt.gcf(), "shap_beeswarm.png")
    plt.close()

    # Log feature importance metrics
    mean_abs_shap = np.abs(shap_values.values).mean(axis=0)
    for feature, importance in zip(X_test.columns, mean_abs_shap):
        mlflow.log_metric(f"shap_{feature}", importance)
```

### 生产 API
```python
class ExplanationService:
    def __init__(self, model_path, explainer_path):
        self.model = joblib.load(model_path)
        self.explainer = joblib.load(explainer_path)

    def predict_with_explanation(self, X):
        prediction = self.model.predict(X)
        shap_values = self.explainer(X)

        return {
            'prediction': prediction[0],
            'base_value': shap_values.base_values[0],
            'feature_contributions': dict(zip(X.columns, shap_values.values[0]))
        }
```

## 参考文档

该技能包括按主题组织的综合参考文档：

### 参考文献/explainers.md
所有解释器类的完整指南：
- `TreeExplainer` - 对基于树的模型的快速、准确的解释
- `DeepExplainer` - 深度学习模型（TensorFlow、PyTorch）
- `KernelExplainer` - 与模型无关（适用于任何模型）
- `LinearExplainer` - 线性模型的快速解释
- `GradientExplainer` - 基于梯度的神经网络
- `PermutationExplainer` - 精确但对于任何模型来说都很慢

包括：构造函数参数、方法、支持的模型、何时使用、示例、性能注意事项。

### 参考文献/plots.md
综合可视化指南：
- **瀑布图** - 个人预测细分
- **蜂群图** - 价值分布的全球重要性
- **条形图** - 清晰的特征重要性摘要
- **散点图** - 特征预测关系和交互
- **力图** - 交互式附加力可视化
- **热图** - 多样本比较网格
- **小提琴图** - 以分布为中心的替代方案
- **决策图** - 多类预测路径

包括：参数、用例、示例、最佳实践、绘图选择指南。

### 参考文献/workflows.md
详细的工作流程和最佳实践：
- 基本模型解释工作流程
- 模型调试和验证
- 特征工程指导
- 模型比较与选择
- 公平性和偏见分析
- 深度学习模型讲解
- 生产部署
- 时间序列模型解释
- 常见的陷阱和解决方案
- 先进技术
- MLOps 集成

包括：分步说明、代码示例、决策标准、故障排除。

### 参考文献/理论.md
理论基础：
- 博弈论中的沙普利值
- 数学公式和属性
- 与其他解释方法的连接（LIME、DeepLIFT等）
- SHAP计算算法（Tree SHAP、Kernel SHAP等）
- 条件期望和基线选择
- 解释SHAP值
- 互动价值
- 理论限制和考虑因素

包括：数学基础、证明、比较、高级主题。

## 使用指南

**何时加载参考文件**：
- 当用户需要有关特定解释器类型或参数的详细信息时加载`explainers.md`
- 当用户需要详细的可视化指导或探索绘图选项时加载`plots.md`
- 当用户执行复杂的多步骤任务（调试、公平性分析、生产部署）时加载`workflows.md`
- 当用户询问理论基础、Shapley 值或数学细节时加载 `theory.md`

**默认方法**（不加载参考）：
- 使用此 SKILL.md 进行基本解释和快速入门
- 提供标准工作流程和通用模式
- 如果需要更多详细信息，可提供参考文件

**加载参考**：
```python
# To load reference files, use the Read tool with appropriate file path:
# /path/to/shap/references/explainers.md
# /path/to/shap/references/plots.md
# /path/to/shap/references/workflows.md
# /path/to/shap/references/theory.md
```

## 最佳实践总结

1. **选择正确的解释器**：尽可能使用专门的解释器（TreeExplainer、DeepExplainer、LinearExplainer）；除非必要，否则避免 KernelExplainer

2. **从全局开始，然后走向局部**：从蜂群图/条形图开始了解整体情况，然后深入研究瀑布图/散点图了解细节

3. **使用多种可视化**：不同的图揭示不同的见解；结合全局（蜂群）+局部（瀑布）+关系（分散）视图

4. **选择合适的背景数据**：使用训练数据中的50-1000个代表性样本

5. **了解模型输出单位**：了解是否解释概率、对数赔率或原始输出

6. **使用领域知识进行验证**：SHA 显示模型行为；使用领域专业知识来解释和验证

7. **优化性能**：用于可视化的示例子集、大型数据集的批处理、生产中的缓存解释器

8. **检查数据泄漏**：特征重要性出乎意料的高可能表明数据质量问题

9. **考虑特征相关性**：使用 TreeExplainer 的相关性感知选项或特征聚类来获取冗余特征

10. **记住SHAP显示关联，而不是因果关系**：使用领域知识进行因果解释

## 安装

```bash
# Basic installation
uv pip install shap

# With visualization dependencies
uv pip install shap matplotlib

# Latest version
uv pip install -U shap
```

**依赖项**：numpy、pandas、scikit-learn、matplotlib、scipy
**可选**：xgboost、lightgbm、tensorflow、torch（取决于模型类型）

## 其他资源

- **官方文档**：https://shap.readthedocs.io/
- **GitHub 存储库**：https://github.com/slundberg/shap
- **原始论文**：Lundberg & Lee (2017) - “解释模型预测的统一方法”
- **《自然 MI 论文》**：Lundberg 等人。 (2020) - “通过可解释的树木人工智能从本地解释到全球理解”

该技能提供了 SHAP 的全面覆盖，以实现所有用例和模型类型的模型可解释性。