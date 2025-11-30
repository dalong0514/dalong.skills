<!-- 此文件由机器翻译自 competing-risks.md -->

# 竞争风险分析

## 概述

当受试者经历几种相互排斥的事件（事件类型）之一时，就会出现竞争风险。当一个事件发生时，它会阻止（“竞争”）其他事件的发生。

### 竞争风险的示例

**医学研究：**
- 癌症死亡与心血管疾病死亡与其他原因死亡
- 癌症研究中的复发与不复发的死亡
- 移植患者的不同类型感染

**其他应用：**
- 工作终止：退休、辞职、因故终止
- 设备故障：不同的故障模式
- 客户流失：不同的离开原因

### 关键概念：累积发生函数 (CIF)

**累积发生率函数 (CIF)** 表示在时间 *t* 之前经历特定事件类型的概率，考虑到竞争风险的存在。

**CIF_k(t) = P(T ≤ t，事件类型 = k)**

这与 Kaplan-Meier 估计器不同，后者在存在竞争风险时会高估事件概率。

## 何时使用竞争风险分析

**在以下情况下使用竞争风险：**
- 存在多种互斥的事件类型
- 一个事件的发生会妨碍其他事件的发生
- 需要估计特定事件类型的概率
- 想要了解协变量如何影响不同的事件类型

**请勿在以下情况下使用：**
- 仅感兴趣一种事件类型（标准生存分析）
- 事件不是相互排斥的（使用循环事件方法）
- 竞争事件极其罕见（可以视为审查）

## 具有竞争风险的累积发生率

###cumulative_incidence_competing_risks 函数

估计每种事件类型的累积发生率函数。

```python
from sksurv.nonparametric import cumulative_incidence_competing_risks
from sksurv.datasets import load_leukemia

# Load data with competing risks
X, y = load_leukemia()
# y has event types: 0=censored, 1=relapse, 2=death

# Compute cumulative incidence for each event type
# Returns: time points, CIF for event 1, CIF for event 2, ...
time_points, cif_1, cif_2 = cumulative_incidence_competing_risks(y)

# Plot cumulative incidence functions
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.step(time_points, cif_1, where='post', label='Relapse', linewidth=2)
plt.step(time_points, cif_2, where='post', label='Death in remission', linewidth=2)
plt.xlabel('Time (weeks)')
plt.ylabel('Cumulative Incidence')
plt.title('Competing Risks: Relapse vs Death')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

### 解释

- **时间 t** 的 CIF：在时间 t 之前经历该特定事件的概率
- **所有 CIF 的总和**：经历任何事件（所有原因）的总概率
- **1 - CIF 总和**：无事件且未经审查的概率

## 竞争风险的数据格式

### 使用事件类型创建结构化数组

<<<代码块_1>>>

### 使用事件类型转换数据

<<<代码块_2>>>

## 比较各组之间的累积发生率

### 分层分析

<<<代码块_3>>>

## 具有竞争风险的统计测试

### 格雷测试

使用格雷检验（在生命线等其他软件包中可用）比较组之间的累积发生率函数。

<<<代码块_4>>>

## 竞争风险建模

### 方法 1：特定原因的危害模型

为每种事件类型拟合单独的 Cox 模型，将其他事件类型视为已审查。

<<<代码块_5>>>

**释义：**
- 每个竞赛项目都有单独的模型
- 系数显示对该事件类型的特定原因危害的影响
- 协变量可能会增加一种事件类型的风险，但会降低另一种事件类型的风险

### 方法 2：细灰色次分布风险模型

直接对累积发生率进行建模（不能直接在 scikit-survival 中使用，但可以使用其他包）。

<<<代码块_6>>>

## 实例：完整的竞争风险分析

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sksurv.nonparametric import cumulative_incidence_competing_risks
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.util import Surv

# Simulate competing risks data
np.random.seed(42)
n = 200

# Create features
age = np.random.normal(60, 10, n)
treatment = np.random.choice(['A', 'B'], n)

# Simulate event times and types
# Event types: 0=censored, 1=relapse, 2=death
times = np.random.exponential(100, n)
event_types = np.zeros(n, dtype=int)

# Higher age increases both events, treatment A reduces relapse
for i in range(n):
    if times[i] < 150:  # Event occurred
        # Probability of each event type
        p_relapse = 0.6 if treatment[i] == 'B' else 0.4
        event_types[i] = 1 if np.random.rand() < p_relapse else 2
    else:
        times[i] = 150  # Censored at study end

# Create DataFrame
df = pd.DataFrame({
    'time': times,
    'event_type': event_types,
    'age': age,
    'treatment': treatment
})

# Encode treatment
df['treatment_A'] = (df['treatment'] == 'A').astype(int)

# 1. OVERALL CUMULATIVE INCIDENCE
print("=" * 60)
print("OVERALL CUMULATIVE INCIDENCE")
print("=" * 60)

y_all = Surv.from_arrays(event=(df['event_type'] > 0), time=df['time'])
time_points, cif_relapse, cif_death = cumulative_incidence_competing_risks(y_all)

plt.figure(figsize=(10, 6))
plt.step(time_points, cif_relapse, where='post', label='Relapse', linewidth=2)
plt.step(time_points, cif_death, where='post', label='Death', linewidth=2)
plt.xlabel('Time (days)')
plt.ylabel('Cumulative Incidence')
plt.title('Competing Risks: Relapse vs Death')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()

print(f"5-year relapse incidence: {cif_relapse[-1]:.2%}")
print(f"5-year death incidence: {cif_death[-1]:.2%}")

# 2. STRATIFIED BY TREATMENT
print("\n" + "=" * 60)
print("CUMULATIVE INCIDENCE BY TREATMENT")
print("=" * 60)

for trt in ['A', 'B']:
    mask = df['treatment'] == trt
    y_trt = Surv.from_arrays(
        event=(df.loc[mask, 'event_type'] > 0),
        time=df.loc[mask, 'time']
    )
    time_trt, cif1_trt, cif2_trt = cumulative_incidence_competing_risks(y_trt)
    print(f"\nTreatment {trt}:")
    print(f"  5-year relapse: {cif1_trt[-1]:.2%}")
    print(f"  5-year death: {cif2_trt[-1]:.2%}")

# 3. CAUSE-SPECIFIC MODELS
print("\n" + "=" * 60)
print("CAUSE-SPECIFIC HAZARD MODELS")
print("=" * 60)

X = df[['age', 'treatment_A']]

# Model for relapse (event type 1)
y_relapse = Surv.from_arrays(
    event=(df['event_type'] == 1),
    time=df['time']
)
cox_relapse = CoxPHSurvivalAnalysis()
cox_relapse.fit(X, y_relapse)

print("\nRelapse Model:")
print(f"  Age:        HR = {np.exp(cox_relapse.coef_[0]):.3f}")
print(f"  Treatment A: HR = {np.exp(cox_relapse.coef_[1]):.3f}")

# Model for death (event type 2)
y_death = Surv.from_arrays(
    event=(df['event_type'] == 2),
    time=df['time']
)
cox_death = CoxPHSurvivalAnalysis()
cox_death.fit(X, y_death)

print("\nDeath Model:")
print(f"  Age:        HR = {np.exp(cox_death.coef_[0]):.3f}")
print(f"  Treatment A: HR = {np.exp(cox_death.coef_[1]):.3f}")

print("\n" + "=" * 60)
```

## 重要考虑因素

### 竞争风险审查

- **行政审查**：受试者在研究结束时仍处于危险之中
- **失访**：受试者在事件发生前离开研究
- **竞争事件**：发生了其他事件 - 未针对 CIF 进行审查，但针对特定原因模型进行审查

### 在特定原因模型和次分布模型之间进行选择

**特定原因的危险模型：**
- 更容易解释
- 对危险率的直接影响
- 更好地了解病因
- 可以与 scikit-survival 配合

**细灰色子分布模型：**
- 直接对累积发生率进行建模
- 更好地进行预测和风险分层
- 更适合临床决策
- 需要其他包

### 常见错误

**错误 1**：使用 Kaplan-Meier 来估计特定事件的概率
- **错误**：Kaplan-Meier 对于事件类型 1，将类型 2 视为审查
- **正确**：考虑竞争风险的累积发生率函数

**错误 2**：当竞争风险很大时就忽略它们
- 如果竞争事件率 > 10-20%，应使用竞争风险方法

**错误 3**：混淆特定原因危害和次分布危害
- 他们回答不同的问题
- 针对您的研究问题使用适当的模型

## 总结

**主要功能：**
- `cumulative_incidence_competing_risks`：估计每种事件类型的 CIF
- 针对特定原因的危险安装单独的 Cox 模型
- 使用分层分析来比较组

**最佳实践：**
1. 始终绘制累积发生函数
2. 报告特定事件和总体发生率
3. 在 scikit-survival 中使用特定原因模型
4. 考虑使用 Fine-Gray 模型（其他包）进行预测
5. 明确哪些事件是竞争事件，哪些是审查事件