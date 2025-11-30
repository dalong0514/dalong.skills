<!-- 此文件由机器翻译自 examples.md -->

# 德纳里奥示例

## 完整的端到端研究示例

此示例展示了从数据到发表的完整研究流程。

### 设置

```python
from denario import Denario, Journal
import os

# Create project directory
os.makedirs("climate_research", exist_ok=True)
den = Denario(project_dir="./climate_research")
```

### 定义研究背景

<<<代码块_1>>>

### 执行完整管道

<<<代码块_2>>>

### 查看输出

<<<代码块_3>>>

## 增强输入描述

改进数据描述以更好地产生想法。

### 基本描述

<<<代码块_4>>>

### 增强描述

<<<代码块_5>>>

## 文献检索整合

将现有研究纳入您的工作流程。

### 示例：查找相关工作

<<<代码块_6>>>

## 从数据中产生研究想法

专注于创意生成，无需完整的流程执行。

### 示例：集思广益研究问题

```python
den = Denario(project_dir="./idea_generation")

# Provide comprehensive data description
den.set_data_description("""
Available datasets:
1. Social media sentiment data (1M tweets, 2020-2023)
2. Stock market prices (S&P 500, daily, 2020-2023)
3. Economic indicators (GDP, unemployment, inflation)

Tools: pandas, sklearn, statsmodels, Prophet, VADER sentiment analysis

Domain: Computational social science and finance
Research interests: Market prediction, sentiment analysis, causal inference
""")

# Generate multiple ideas (conceptual - depends on denario API)
den.get_idea()

# Review the generated idea in idea.md
# Decide whether to proceed or regenerate
```

## 根据现有结果撰写论文

当分析已经完成时，使用第纳尔生成论文。

### 示例：格式化现有研究

```python
den = Denario(project_dir="./paper_generation")

# Provide all components manually
den.set_data_description("""
Completed analysis of traffic pattern data from urban sensors
Dataset: 6 months of traffic flow measurements from 100 intersections
Analysis completed using R and Python
""")

den.set_idea("""
Research question: Optimize traffic light timing using reinforcement learning
to reduce congestion and improve traffic flow efficiency
""")

den.set_method("""
# Methodology

## Data Collection
Traffic flow data collected from 100 intersections in downtown area from
January-June 2023. Measurements include vehicle counts, wait times, and
queue lengths at 1-minute intervals.

## Model Development
Developed a Deep Q-Network (DQN) reinforcement learning agent to optimize
traffic light timing. State space includes current queue lengths and
historical flow patterns. Actions correspond to light timing adjustments.

## Training
Trained the agent using historical data with a reward function based on
total wait time reduction. Used experience replay and target networks for
stable learning.

## Validation
Validated using held-out test data and compared against:
- Current fixed-timing system
- Actuated control system
- Alternative RL algorithms (A3C, PPO)

## Metrics
- Average wait time reduction
- Total throughput improvement
- Queue length distribution
- Computational efficiency
""")

den.set_results("""
# Results

## Training Performance
The DQN agent converged after 500,000 training episodes. Training time: 12 hours
on NVIDIA V100 GPU.

## Wait Time Reduction
- Current system: Average wait time 45.2 seconds
- DQN system: Average wait time 32.8 seconds
- Improvement: 27.4% reduction (p < 0.001)

## Throughput Analysis
- Vehicles processed per hour increased from 2,850 to 3,420 (+20%)
- Peak hour congestion reduced by 35%

## Comparison with Baselines
- Actuated control: 38.1 seconds average wait (DQN still 14% better)
- A3C: 34.5 seconds (DQN slightly better, 5%)
- PPO: 33.2 seconds (DQN marginally better, 1%)

## Queue Length Analysis
Maximum queue length reduced from 42 vehicles to 28 vehicles during peak hours.

## Figures
- Figure 1: Training curve showing convergence
- Figure 2: Wait time distribution comparison
- Figure 3: Throughput over time of day
- Figure 4: Heatmap of queue lengths across intersections
""")

# Generate publication-ready paper
den.get_paper(journal=Journal.APS)
```

## Gemini 的快速模式

使用 Google 的 Gemini 模型可以加快执行速度。

### 示例：快速原型设计

```python
# Configure for fast mode (conceptual - check denario documentation)
# This would involve setting appropriate LLM backend

den = Denario(project_dir="./fast_research")

# Same workflow, optimized for speed
den.set_data_description("""
Quick analysis needed: Monthly sales data (2 years)
Goal: Identify seasonal patterns and forecast next quarter
Tools: pandas, Prophet
""")

# Fast execution
den.get_idea()
den.get_method()
den.get_results()
den.get_paper()

# Trade-off: Faster execution, potentially less detailed analysis
```

## 混合工作流程：自定义想法 + 自动化方法

结合手动和自动方法。

### 示例：定向研究

```python
den = Denario(project_dir="./hybrid_workflow")

# Describe data
den.set_data_description("""
Medical imaging dataset: 10,000 chest X-rays
Labels: Normal, pneumonia, COVID-19
Format: 224x224 grayscale PNG files
Tools: TensorFlow, Keras, scikit-learn, OpenCV
""")

# Provide specific research direction
den.set_idea("""
Develop a transfer learning approach using pre-trained ResNet50 for multi-class
classification of chest X-rays, with focus on interpretability using Grad-CAM
to identify diagnostic regions
""")

# Let denario develop the methodology
den.get_method()

# Review methodology, then execute
den.get_results()

# Generate paper
den.get_paper(journal=Journal.APS)
```

## 时间序列分析示例

时态数据的专门示例。

### 示例：经济预测

```python
den = Denario(project_dir="./time_series_analysis")

den.set_data_description("""
Dataset: Monthly unemployment rates (US, 1950-2023)
Additional features: GDP growth, inflation, interest rates
Format: Multivariate time-series DataFrame
Tools: statsmodels, Prophet, pmdarima, sklearn

Analysis goals:
- Model unemployment trends
- Forecast next 12 months
- Identify leading indicators
- Assess forecast uncertainty

Data characteristics:
- Seasonal patterns (annual cycles)
- Structural breaks (recessions)
- Autocorrelation present
- Non-stationary (unit root)
""")

den.get_idea()
# Might generate: "Develop a SARIMAX model incorporating economic indicators
# as exogenous variables to forecast unemployment with confidence intervals"

den.get_method()
den.get_results()
den.get_paper(journal=Journal.APS)
```

## 机器学习管道示例

通过验证完成 ML 工作流程。

### 示例：预测建模

```python
den = Denario(project_dir="./ml_pipeline")

den.set_data_description("""
Dataset: Customer churn prediction
- 50,000 customers, 30 features (demographics, usage patterns, service history)
- Binary target: churned (1) or retained (0)
- Imbalanced: 20% churn rate
- Features: Numerical and categorical mixed

Available tools:
- pandas for preprocessing
- sklearn for modeling (RF, XGBoost, logistic regression)
- imblearn for handling imbalance
- SHAP for feature importance

Goals:
- Build predictive model for churn
- Identify key churn factors
- Provide actionable insights
- Achieve >85% AUC-ROC
""")

den.get_idea()
# Might generate: "Develop an ensemble model combining XGBoost and Random Forest
# with SMOTE oversampling, and use SHAP values to identify interpretable
# churn risk factors"

den.get_method()
# Will include: train/test split, cross-validation, hyperparameter tuning,
# performance metrics, feature importance analysis

den.get_results()
# Executes full ML pipeline, generates:
# - Model performance metrics
# - ROC curves
# - Feature importance plots
# - Confusion matrices

den.get_paper(journal=Journal.APS)
```

## 有效使用技巧

### 提供丰富的上下文

更多背景→更好的想法和方法：

```python
# Include:
# - Data characteristics (size, format, quality issues)
# - Available tools and libraries
# - Domain-specific knowledge
# - Research objectives and constraints
# - Known challenges or considerations
```

### 迭代中间输出

在每个阶段进行审查和完善：

```python
# Generate
den.get_idea()

# Review idea.md
# If needed, refine:
den.set_idea("Refined version of the idea")

# Continue
den.get_method()
# Review methodology.md
# Refine if needed, then proceed
```

### 保存您的工作流程

记录完整的管道：

```python
# Save workflow script
with open("research_workflow.py", "w") as f:
    f.write("""
from denario import Denario, Journal

den = Denario(project_dir="./project")
den.set_data_description("...")
den.get_idea()
den.get_method()
den.get_results()
den.get_paper(journal=Journal.APS)
""")
```

### 使用版本控制

追踪研究进展：

```bash
cd project_dir
git init
git add .
git commit -m "Initial data description"

# After each stage
git add .
git commit -m "Generated research idea"
# ... continue committing after each stage
```