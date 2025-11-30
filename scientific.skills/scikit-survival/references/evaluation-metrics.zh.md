<!-- 此文件由机器翻译自 evaluation-metrics.md -->

# 生存模型的评估指标

## 概述

评估生存模型需要考虑审查数据的专门指标。 scikit-survival 提供了三大类指标：
1. 一致性指数（C-index）
2. 时间依赖性 ROC 和 AUC
3. 布赖尔分数

## 一致性指数（C-index）

### 它衡量什么

一致性指数衡量预测风险评分和观察到的事件时间之间的等级相关性。它表示对于一对随机受试者，模型正确排序其生存时间的概率。

**范围**：0 到 1
- 0.5 = 随机预测
- 1.0 = 完美的一致性
- 典型良好性能：0.7-0.8

### 两个实现

#### Harrell 的 C 指数 (concordance_index_censored)

传统的估计器更简单但有局限性。

**何时使用：**
- 低审查率（< 40%）
- 开发过程中快速评估
- 比较同一数据集上的模型

**限制：**
- 由于高审查率而变得越来越有偏见
- 从大约 49% 的审查开始高估了性能

```python
from sksurv.metrics import concordance_index_censored

# Compute Harrell's C-index
result = concordance_index_censored(y_test['event'], y_test['time'], risk_scores)
c_index = result[0]
print(f"Harrell's C-index: {c_index:.3f}")
```

#### Uno 的 C 索引 (concordance_index_ipcw)

纠正审查偏差的审查加权逆概率 (IPCW) 估计器。

**何时使用：**
- 中等到高审查率（> 40%）
- 需要公正的估计
- 比较不同数据集的模型
- 发布结果（更稳健）

**优点：**
- 即使在高审查的情况下也保持稳定
- 更可靠的估计
- 减少偏见

<<<代码块_1>>>

### 在 Harrell's 和 Uno's 之间进行选择

**在以下情况下使用 Uno 的 C 指数：**
- 审查率 > 40%
- 需要最准确的估计
- 比较不同研究的模型
- 出版研究

**在以下情况下使用 Harrell 的 C 指数：**
- 低审查率
- 开发过程中的快速模型比较
- 计算效率至关重要

### 示例比较

<<<代码块_2>>>

## 时间相关的 ROC 和 AUC

### 它衡量什么

时间相关的 AUC 评估特定时间点的模型辨别力。它将在时间 *t* 经历事件的主体与没有经历事件的主体区分开来。

**回答的问题**：“模型预测谁将在时间 t 之前发生事件的效果如何？”

### 何时使用

- 预测特定时间窗口内事件的发生
- 特定时间点的临床决策（例如 5 年生存率）
- 想要评估不同时间范围内的绩效
- 需要区分和时间信息

### 关键函数：cumulative_dynamic_auc

<<<代码块_3>>>

### 解释

- **时间 t** 的 AUC：概率模型正确地将在时间 t 发生过事件的受试者排在没有发生过事件的受试者之上
- **随时间变化的 AUC**：表示模型性能随时间范围的变化
- **平均 AUC**：所有时间点歧视的总体总结

### 示例：比较模型

<<<代码块_4>>>

## 布赖尔分数

### 它衡量什么

Brier 评分将均方误差扩展到带有审查的生存数据。它衡量歧视（排名）和校准（预测概率的准确性）。

**公式**：**(1/n) Σ (S(t|x_i) - I(T_i > t))²**

其中 S(t|x_i) 是受试者 i 在时间 t 的预测生存概率。

**范围**：0 到 1
- 0 = 完美预测
- 越低越好
- 典型良好性能：< 0.2

### 何时使用

- 需要校准评估（不仅仅是排名）
- 想要评估预测概率，而不仅仅是风险评分
- 比较输出生存函数的模型
- 需要概率估计的临床应用

### 关键功能

#### brier_score：单个时间点

<<<代码块_5>>>

#### Integrated_brier_score：跨时间总结

<<<代码块_6>>>

### 解释

- **时间 t 时的 Brier 分数**：时间 t 时预测生存率与实际生存率之间的预期平方差
- **综合 Brier 分数**：不同时间段 Brier 分数的加权平均值
- **较低的值=更好的预测**

### 与空模型的比较

始终与基线进行比较（例如 Kaplan-Meier）：

```python
from sksurv.nonparametric import kaplan_meier_estimator

# Compute Kaplan-Meier baseline
time_km, surv_km = kaplan_meier_estimator(y_train['event'], y_train['time'])

# Predict with KM for each test subject
surv_km_test = [surv_km[time_km <= time_point][-1] if any(time_km <= time_point) else 1.0
                for _ in range(len(X_test))]

bs_km = brier_score(y_train, y_test, surv_km_test, time_point)[1]
bs_model = brier_score(y_train, y_test, surv_at_t, time_point)[1]

print(f"Kaplan-Meier Brier Score: {bs_km:.3f}")
print(f"Model Brier Score: {bs_model:.3f}")
print(f"Improvement: {(bs_km - bs_model) / bs_km * 100:.1f}%")
```

## 使用指标进行交叉验证

### 一致性指数评分器

```python
from sklearn.model_selection import cross_val_score
from sksurv.metrics import as_concordance_index_ipcw_scorer

# Create scorer
scorer = as_concordance_index_ipcw_scorer()

# Perform cross-validation
scores = cross_val_score(model, X, y, cv=5, scoring=scorer)
print(f"Mean C-index: {scores.mean():.3f} (±{scores.std():.3f})")
```

### 综合 Brier Score 评分器

```python
from sksurv.metrics import as_integrated_brier_score_scorer

# Define time points for evaluation
times = np.percentile(y['time'][y['event']], [25, 50, 75])

# Create scorer
scorer = as_integrated_brier_score_scorer(times)

# Perform cross-validation
scores = cross_val_score(model, X, y, cv=5, scoring=scorer)
print(f"Mean IBS: {scores.mean():.3f} (±{scores.std():.3f})")
```

## 使用 GridSearchCV 进行模型选择

```python
from sklearn.model_selection import GridSearchCV
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import as_concordance_index_ipcw_scorer

# Define parameter grid
param_grid = {
    'n_estimators': [100, 200, 300],
    'min_samples_split': [10, 20, 30],
    'max_depth': [None, 10, 20]
}

# Create scorer
scorer = as_concordance_index_ipcw_scorer()

# Perform grid search
cv = GridSearchCV(
    RandomSurvivalForest(random_state=42),
    param_grid,
    scoring=scorer,
    cv=5,
    n_jobs=-1
)
cv.fit(X, y)

print(f"Best parameters: {cv.best_params_}")
print(f"Best C-index: {cv.best_score_:.3f}")
```

## 综合模型评估

### 推荐的评估流程

```python
from sksurv.metrics import (
    concordance_index_censored,
    concordance_index_ipcw,
    cumulative_dynamic_auc,
    integrated_brier_score
)

def evaluate_survival_model(model, X_train, X_test, y_train, y_test):
    """Comprehensive evaluation of survival model"""

    # Get predictions
    risk_scores = model.predict(X_test)
    surv_funcs = model.predict_survival_function(X_test)

    # 1. Concordance Index (both versions)
    c_harrell = concordance_index_censored(y_test['event'], y_test['time'], risk_scores)[0]
    c_uno = concordance_index_ipcw(y_train, y_test, risk_scores)[0]

    # 2. Time-dependent AUC
    times = np.percentile(y_test['time'][y_test['event']], [25, 50, 75])
    auc, mean_auc = cumulative_dynamic_auc(y_train, y_test, risk_scores, times)

    # 3. Integrated Brier Score
    ibs = integrated_brier_score(y_train, y_test, surv_funcs, times)

    # Print results
    print("=" * 50)
    print("Model Evaluation Results")
    print("=" * 50)
    print(f"Harrell's C-index:  {c_harrell:.3f}")
    print(f"Uno's C-index:      {c_uno:.3f}")
    print(f"Mean AUC:           {mean_auc:.3f}")
    print(f"Integrated Brier:   {ibs:.3f}")
    print("=" * 50)

    return {
        'c_harrell': c_harrell,
        'c_uno': c_uno,
        'mean_auc': mean_auc,
        'ibs': ibs,
        'time_auc': dict(zip(times, auc))
    }

# Use the evaluation function
results = evaluate_survival_model(model, X_train, X_test, y_train, y_test)
```

## 选择正确的指标

### 决策指南

**在以下情况下使用 C 指数（Uno's）：**
- 主要目标是排名/歧视
- 不需要校准概率
- 标准生存分析设置
- 最常见的选择

**在以下情况下使用时间相关 AUC：**
- 在特定时间点需要区分
- 特定范围内的临床决策
- 想要了解性能如何随时间变化

**在以下情况下使用 Brier 分数：**
- 需要校准概率估计
- 辨别和校准都很重要
- 需要概率的临床决策
- 想要综合评估

**最佳实践**：报告多个指标以进行综合评估。至少报告：
- Uno 的 C 指数（歧视）
- 综合 Brier 评分（辨别 + 校准）
- 临床相关时间点的时间依赖性 AUC