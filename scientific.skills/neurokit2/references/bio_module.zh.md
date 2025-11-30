<!-- 此文件由机器翻译自 bio_module.md -->

# 多信号集成（生物模块）

## 概述

Bio模块提供统一的功能，可同时处理和分析多个生理信号。它充当协调信号特定处理功能并实现集成多模态分析的包装器。

## 多信号处理

### bio_process()

通过单个函数调用同时处理多个生理信号。

```python
bio_signals, bio_info = nk.bio_process(ecg=None, rsp=None, eda=None, emg=None,
                                       ppg=None, eog=None, sampling_rate=1000)
```

**参数：**
- `ecg`：ECG 信号数组（可选）
- `rsp`：呼吸信号数组（可选）
- `eda`：EDA 信号数组（可选）
- `emg`：EMG 信号阵列（可选）
- `ppg`：PPG 信号数组（可选）
- `eog`：EOG 信号数组（可选）
- `sampling_rate`：采样率（以 Hz 为单位）（信号间必须一致或按信号指定）

**退货：**
- `bio_signals`：包含所有已处理信号的统一数据帧，其列为：
  - 信号特定功能（例如，`ECG_Clean`、`ECG_Rate`、`EDA_Phasic`、`RSP_Rate`）
  - 所有检测到的事件/峰值
  - 衍生措施
- `bio_info`：包含信号特定信息的字典（峰值位置、参数）

**示例：**
<<<代码块_1>>>

**内部工作流程：**
1. 每个信号均由其专用处理函数处理：
   - `ecg_process()` 用于心电图
   - `rsp_process()` 用于呼吸
   - 用于 EDA 的 `eda_process()`
   - `emg_process()` 用于肌电图
   - PPG `ppg_process()`
   - 用于 EOG 的 `eog_process()`
2.结果合并到统一的DataFrame中
3. 计算交叉信号特征（例如，如果 ECG 和 RSP 都存在，则为 RSA）

**优点：**
- 简化的多模式记录 API
- 所有信号的统一时基
- 自动交叉信号特征计算
- 一致的输出格式

## 多信号分析

### bio_analyze()

对处理后的多模态信号进行综合分析。

<<<代码块_2>>>

**参数：**
- `bio_signals`：来自 `bio_process()` 的 DataFrame 或自定义处理的信号
- `sampling_rate`：采样率 (Hz)

**退货：**
- DataFrame 包含所有检测到的信号类型的分析结果：
  - 如果持续时间 ≥ 10 秒，则与间隔相关的指标
  - 事件相关指标（如果持续时间 < 10 秒）
  - 交叉信号指数（例如，如果 ECG + RSP 可用，则为 RSA）

**按信号计算的指标：**
- **ECG**：心率统计、HRV 指数（时间、频率、非线性域）
- **RSP**：呼吸频率统计、RRV、振幅测量
- **EDA**：SCR 计数、振幅、紧张水平、交感指数
- **EMG**：激活计数、幅度统计
- **PPG**：类似于 ECG（心率、HRV）
- **EOG**：眨眼次数、眨眼率

**交叉信号指标：**
- **RSA（呼吸性窦性心律失常）**：如果存在 ECG + RSP
- **心肺耦合**：相位同步指数
- **多模式唤醒**：组合自主神经指数

**示例：**
<<<代码块_3>>>

## 交叉信号特性

当多个信号一起处理时，NeuroKit2 可以计算集成特征：

### 呼吸性窦性心律失常 (RSA)

当心电图和呼吸信号同时存在时自动计算。

<<<代码块_4>>>

**计算：**
- 通过呼吸进行高频 HRV 调制
- 需要同步心电图 R 峰值和呼吸信号
- 方法：Porges-Bohrer 或峰谷法

**释义：**
- RSA 越高：副交感神经（迷走神经）影响越大
- 心肺耦合标记
- 健康指标和情绪调节能力

### 心电图衍生呼吸 (EDR)

如果呼吸信号不可用，NeuroKit2 可以根据心电图进行估计：

<<<代码块_5>>>

**用例：**
- 当无法直接测量时估计呼吸
- 交叉验证呼吸测量结果

### Cardio-EDA 集成

同步心脏和皮肤电活动：

<<<代码块_6>>>

## 实际工作流程

### 完整的多模态录音分析

```python
import neurokit2 as nk
import pandas as pd

# 1. Load multi-modal physiological data
ecg = load_ecg()        # Your data loading function
rsp = load_rsp()
eda = load_eda()
emg = load_emg()

# 2. Process all signals simultaneously
bio_signals, bio_info = nk.bio_process(
    ecg=ecg,
    rsp=rsp,
    eda=eda,
    emg=emg,
    sampling_rate=1000
)

# 3. Visualize processed signals
import matplotlib.pyplot as plt

fig, axes = plt.subplots(4, 1, figsize=(15, 12), sharex=True)

# ECG
axes[0].plot(bio_signals.index / 1000, bio_signals['ECG_Clean'])
axes[0].set_ylabel('ECG')
axes[0].set_title('Multi-Modal Physiological Recording')

# Respiration
axes[1].plot(bio_signals.index / 1000, bio_signals['RSP_Clean'])
axes[1].set_ylabel('Respiration')

# EDA
axes[2].plot(bio_signals.index / 1000, bio_signals['EDA_Phasic'])
axes[2].set_ylabel('EDA (Phasic)')

# EMG
axes[3].plot(bio_signals.index / 1000, bio_signals['EMG_Amplitude'])
axes[3].set_ylabel('EMG Amplitude')
axes[3].set_xlabel('Time (s)')

plt.tight_layout()
plt.show()

# 4. Analyze all signals
results = nk.bio_analyze(bio_signals, sampling_rate=1000)

# 5. Extract key metrics
print("Heart Rate (mean):", results['ECG_Rate_Mean'])
print("HRV (RMSSD):", results['HRV_RMSSD'])
print("Breathing Rate:", results['RSP_Rate_Mean'])
print("SCRs (count):", results['SCR_Peaks_N'])
print("RSA:", results['RSA'])
```

### 事件相关的多模态分析

```python
# 1. Process signals
bio_signals, bio_info = nk.bio_process(ecg=ecg, rsp=rsp, eda=eda, sampling_rate=1000)

# 2. Detect events
events = nk.events_find(trigger_channel, threshold=0.5)

# 3. Create epochs for all signals
epochs = nk.epochs_create(bio_signals, events, sampling_rate=1000,
                          epochs_start=-1.0, epochs_end=10.0,
                          event_labels=event_labels,
                          event_conditions=event_conditions)

# 4. Signal-specific event-related analysis
ecg_eventrelated = nk.ecg_eventrelated(epochs)
rsp_eventrelated = nk.rsp_eventrelated(epochs)
eda_eventrelated = nk.eda_eventrelated(epochs)

# 5. Merge results
all_results = pd.merge(ecg_eventrelated, rsp_eventrelated,
                       left_index=True, right_index=True)
all_results = pd.merge(all_results, eda_eventrelated,
                       left_index=True, right_index=True)

# 6. Statistical comparison by condition
all_results['Condition'] = event_conditions
condition_means = all_results.groupby('Condition').mean()
```

### 不同的采样率

处理具有不同本机采样率的信号：

```python
# ECG at 1000 Hz, EDA at 100 Hz
bio_signals, bio_info = nk.bio_process(
    ecg=ecg_1000hz,
    eda=eda_100hz,
    sampling_rate=1000  # Target sampling rate
)
# EDA will be automatically resampled to 1000 Hz internally
```

或者单独处理并合并：

```python
# Process at native sampling rates
ecg_signals, ecg_info = nk.ecg_process(ecg, sampling_rate=1000)
eda_signals, eda_info = nk.eda_process(eda, sampling_rate=100)

# Resample to common rate
eda_resampled = nk.signal_resample(eda_signals, sampling_rate=100,
                                   desired_sampling_rate=1000)

# Merge manually
bio_signals = pd.concat([ecg_signals, eda_resampled], axis=1)
```

## 用例和应用

### 综合心理生理学研究

捕捉生理唤醒的多个维度：

- **心脏**：定向、注意力、情绪效价
- **呼吸**：唤醒、压力、情绪调节
- **EDA**：交感神经唤醒，情绪强度
- **肌电图**：肌肉紧张、面部表情、惊吓

**示例：情感图片观看**
- 心电图：观看图片时心率减慢（注意）
- EDA：SCR 反映情绪唤醒强度
- RSP：屏住呼吸或变化反映了情感投入
- 面部肌电图：皱眉肌（皱眉）、颧肌（微笑）的效价

### 压力和放松评估

多模态标记提供聚合证据：

- **增加压力**：↑ HR、↓ HRV、↑ EDA、↑ 呼吸频率、↑ 肌肉张力
- **放松**：↓ HR、↑ HRV、↓ EDA、↓ 呼吸频率、缓慢呼吸、↓ 肌肉紧张

**干预效果：**
- 比较干预前与干预后的多模式指数
- 确定哪些模式响应特定技术

### 临床评估

**焦虑症：**
- 提高基线 EDA、HR
- 对压力源的过度反应
- HRV、呼吸变异性降低

**抑郁症：**
- 改变自主平衡（↓ HRV）
- EDA 反应迟钝
- 呼吸模式不规则

**创伤后应激障碍：**
- 过度唤醒（↑ HR，↑ EDA 基线）
- 过度惊吓（肌电图）
- 修改 RSA

### 人机交互

不引人注目的用户状态评估：

- **认知负荷**：↓ HRV，↑ EDA，抑制眨眼
- **沮丧**：↑ HR、↑ EDA、↑ 肌肉紧张
- **参与度**：适度的唤醒，同步的反应
- **无聊**：低唤醒度、不规则模式

### 运动表现和恢复

监控训练负荷和恢复：

- **静息 HRV**：每日监测过度训练
- **EDA**：交感神经激活和压力
- **呼吸**：运动/恢复期间的呼吸模式
- **多模式集成**：全面的恢复评估

## 多模式记录的优点

**收敛效度：**
- 多个指数收敛于同一构造（例如，唤醒）
- 比单一措施更稳健

**区分效度：**
- 不同的信号在特定条件下会分离
- 心电图反映交感神经和副交感神经
- EDA主要反映同情心

**系统集成：**
- 了解全身生理协调
- 交叉信号耦合指标（RSA、一致性）

**冗余和稳健性：**
- 如果其中一个信号质量较差，可以使用其他信号
- 跨模式交叉验证结果

**更丰富的解读：**
- HR 减速 + SCR 增加 = 唤醒定向
- HR 加速 + 无 SCR = 无交感神经兴奋的心脏反应

## 注意事项

### 硬件和同步

- **同一设备**：信号本质上同步
- **不同设备**：需要共同的触发器/时间戳
  - 使用硬件触发器来标记同时发生的事件
  - 基于事件标记的软件对齐
  - 验证同步质量（互相关冗余信号）

### 跨模式的信号质量

- 并非所有信号都具有相同的质量
- 根据研究问题确定优先级
- 记录每个信号的质量问题

### 计算成本

- 处理多个信号会增加计算时间
- 考虑对大型数据集进行批量处理
- 适当降低采样以减少负载

### 分析复杂性

- 更多信号=更多变量=更多统计比较
- 未经纠正，存在第一类错误（误报）的风险
- 使用多变量方法或预先注册的分析

### 解释

- 避免过度解释复杂的多模式模式
- 生理学理论基础
- 在强烈主张之前重复研究结果

## 参考文献

- Berntson, G. G.、Cacioppo, J. T. 和 Quigley, K. S. (1993)。呼吸性窦性心律失常：自主神经起源、生理机制和心理生理学影响。心理生理学，30(2), 183-196。
- Cacioppo, J. T.、Tassinary, L. G. 和 Berntson, G.（编）。 （2017）。心理生理学手册（第四版）。剑桥大学出版社。
- Kreibig, S.D. (2010)。情绪中的自主神经系统活动：综述。生物心理学，84(3), 394-421。
- Laborde, S.、Mosley, E. 和 Thayer, J. F. (2017)。心理生理学研究中的心率变异性和心脏迷走神经张力——实验计划、数据分析和数据报告的建议。心理学前沿，8, 213。