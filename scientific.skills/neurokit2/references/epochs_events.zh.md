<!-- 此文件由机器翻译自 epochs_events.md -->

# 纪元和事件相关分析

## 概述

事件相关分析检查对特定刺激或事件的时间锁定的生理反应。 NeuroKit2 提供了用于所有信号类型的事件检测、纪元创建、平均和事件相关特征提取的工具。

## 事件检测

### events_find()

根据阈值交叉或变化自动检测信号中的事件/触发器。

```python
events = nk.events_find(event_channel, threshold=0.5, threshold_keep='above',
                        duration_min=1, inter_min=0)
```

**参数：**
- `threshold`：检测阈值
- `threshold_keep`：`'above'` 或 `'below'` 阈值
- `duration_min`：要保留的最短事件持续时间（样本）
- `inter_min`：事件之间的最小间隔（样本）

**退货：**
- 字典：
  - `'onset'`：事件开始索引
  - `'offset'`：事件偏移索引（如果适用）
  - `'duration'`：事件持续时间
  - `'label'`：事件标签（如果有多个事件类型）

**常见用例：**

**实验中的 TTL 触发器：**
<<<代码块_1>>>

**按下按钮：**
<<<代码块_2>>>

**状态变化：**
<<<代码块_3>>>

### events_plot()

可视化相对于信号的事件时序。

<<<代码块_4>>>

**显示：**
- 信号追踪
- 事件标记（垂直线或阴影区域）
- 事件标签

**用例：**
- 验证事件检测的准确性
- 检查事件的时间分布
- 划时代前的质量控制

## 纪元创建

### epochs_create()

围绕事件创建数据纪元（段）以进行事件相关分析。

<<<代码块_5>>>

**参数：**
- `data`：带有信号或单个信号的 DataFrame
- `events`：来自 `events_find()` 的事件索引或字典
- `sampling_rate`：信号采样率（Hz）
- `epochs_start`：相对于事件的开始时间（秒，负数 = 之前）
- `epochs_end`：相对于事件的结束时间（秒，正=之后）
- `event_labels`：每个事件的标签列表（可选）
- `event_conditions`：每个事件的条件名称列表（可选）
- `baseline_correction`：如果为真，则从每个时期减去基线平均值

**退货：**
- 数据帧字典，每个时期一个
- 每个数据帧包含与事件相关的时间信号数据（事件开始时索引 = 0）
- 包括 `'Label'` 和 `'Condition'` 列（如果提供）

**典型纪元窗口：**
- **视觉 ERP**：-0.2 至 1.0 秒（200 毫秒基线，刺激后 1 秒）
- **心脏定向**：-1.0 至 10 秒（捕获预期和响应）
- **肌电图惊吓**：-0.1 至 0.5 秒（短暂响应）
- **EDA SCR**：-1.0 至 10 秒（1-3 秒延迟，缓慢恢复）

### 活动标签和条件

按类型和实验条件组织活动：

<<<代码块_6>>>

**访问纪元：**
```python
# Epoch by number
epoch_1 = epochs['1']

# Filter by condition
positive_epochs = {k: v for k, v in epochs.items() if v['Condition'][0] == 'positive'}
```

### 基线校正

从时期中删除刺激前基线以隔离与事件相关的变化：

**自动（纪元创建期间）：**
```python
epochs = nk.epochs_create(data, events, sampling_rate=1000,
                          epochs_start=-0.5, epochs_end=2.0,
                          baseline_correction=True)  # Subtracts mean of entire baseline
```

**手动（纪元创建后）：**
```python
# Subtract baseline period mean
baseline_start = -0.5  # seconds
baseline_end = 0.0     # seconds

for key, epoch in epochs.items():
    baseline_mask = (epoch.index >= baseline_start) & (epoch.index < baseline_end)
    baseline_mean = epoch[baseline_mask].mean()
    epochs[key] = epoch - baseline_mean
```

**何时校正基线：**
- **ERP**：始终（隔离与事件相关的更改）
- **心脏/EDA**：通常（消除个体间基线差异）
- **绝对测量**：有时不需要（例如，分析绝对幅度）

## 时代分析和可视化

### epochs_plot()

可视化单个或平均时期。

```python
nk.epochs_plot(epochs, column='ECG_Rate', condition=None, show=True)
```

**参数：**
- `column`：要绘制哪个信号列
- `condition`：仅绘制特定条件（可选）

**显示：**
- 个人纪元痕迹（半透明）
- 跨时期的平均值（粗线）
- 可选：阴影错误（SEM 或 SD）

**使用案例：**
- 可视化与事件相关的响应
- 比较条件
- 识别异常时期

### epochs_average()

通过统计数据计算跨时期的总平均值。

```python
average_epochs = nk.epochs_average(epochs, output='dict')
```

**参数：**
- `output`：`'dict'`（默认）或`'df'`（DataFrame）

**退货：**
- 字典或数据框：
  - `'Mean'`：每个时间点的历元平均值
  - `'SD'`：标准差
  - `'SE'`：平均值的标准误差
  - `'CI_lower'`、`'CI_upper'`：95% 置信区间

**用例：**
- 计算事件相关电位 (ERP)
- 总平均心脏/EDA/EMG 反应
- 集团层面分析

**特定条件平均：**
```python
# Separate averages by condition
positive_epochs = {k: v for k, v in epochs.items() if v['Condition'][0] == 'positive'}
negative_epochs = {k: v for k, v in epochs.items() if v['Condition'][0] == 'negative'}

avg_positive = nk.epochs_average(positive_epochs)
avg_negative = nk.epochs_average(negative_epochs)
```

### epochs_to_df()

将 epochs 字典转换为统一的 DataFrame。

```python
epochs_df = nk.epochs_to_df(epochs)
```
**退货：**
- 堆叠所有纪元的单个 DataFrame
- 包括 `'Epoch'`、`'Time'`、`'Label'`、`'Condition'` 列
- 促进使用 pandas/seaborn 进行统计分析和绘图

**用例：**
- 为混合效应模型准备数据
- 使用seaborn/plotly绘图
- 导出到R或统计软件

### epochs_to_array()

将历元转换为 3D NumPy 数组。

```python
epochs_array = nk.epochs_to_array(epochs, column='ECG_Rate')
```

**退货：**
- 3D 数组：(n_epochs, n_timepoints, n_columns)

**用例：**
- 机器学习输入（划时代的功能）
- 基于定制阵列的分析
- 阵列数据的统计测试

## 信号特定事件相关分析

NeuroKit2 为每种信号类型提供专门的事件相关分析：

### 心电图事件相关
```python
ecg_epochs = nk.epochs_create(ecg_signals, events, sampling_rate=1000,
                              epochs_start=-1, epochs_end=10)
ecg_results = nk.ecg_eventrelated(ecg_epochs)
```

**计算指标：**
- `ECG_Rate_Baseline`：事件前的心率
- `ECG_Rate_Min/Max`：纪元期间的最小/最大速率
- `ECG_Phase_*`：事件发生时的心脏阶段
- 评价跨时间窗口的动态

### EDA 事件相关
```python
eda_epochs = nk.epochs_create(eda_signals, events, sampling_rate=100,
                              epochs_start=-1, epochs_end=10)
eda_results = nk.eda_eventrelated(eda_epochs)
```

**计算指标：**
- `EDA_SCR`：存在 SCR（二进制）
- `SCR_Amplitude`：最大 SCR 幅度
- `SCR_Latency`：SCR 发作时间
- `SCR_RiseTime`、`SCR_RecoveryTime`
- `EDA_Tonic`：平均滋补水平

### RSP 事件相关
```python
rsp_epochs = nk.epochs_create(rsp_signals, events, sampling_rate=100,
                              epochs_start=-0.5, epochs_end=5)
rsp_results = nk.rsp_eventrelated(rsp_epochs)
```

**计算指标：**
- `RSP_Rate_Mean`：平均呼吸频率
- `RSP_Amplitude_Mean`：平均呼吸深度
- `RSP_Phase`：事件时的呼吸阶段
- 速率/幅度动态

### 肌电图事件相关
```python
emg_epochs = nk.epochs_create(emg_signals, events, sampling_rate=1000,
                              epochs_start=-0.1, epochs_end=1.0)
emg_results = nk.emg_eventrelated(emg_epochs)
```

**计算指标：**
- `EMG_Activation`：存在激活
- `EMG_Amplitude_Mean/Max`：幅度统计
- `EMG_Onset_Latency`：激活开始的时间
- `EMG_Bursts`：激活突发数

### EOG 事件相关
```python
eog_epochs = nk.epochs_create(eog_signals, events, sampling_rate=500,
                              epochs_start=-0.5, epochs_end=2.0)
eog_results = nk.eog_eventrelated(eog_epochs)
```

**计算指标：**
- `EOG_Blinks_N`：纪元期间的闪烁次数
- `EOG_Rate_Mean`：闪烁率
- 时间眨眼分布

### PPG 活动相关
```python
ppg_epochs = nk.epochs_create(ppg_signals, events, sampling_rate=100,
                              epochs_start=-1, epochs_end=10)
ppg_results = nk.ppg_eventrelated(ppg_epochs)
```

**计算指标：**
- 类似于心电图：速率动态、相位信息

## 实际工作流程

### 完整的事件相关分析管道

```python
import neurokit2 as nk

# 1. Process physiological signals
ecg_signals, ecg_info = nk.ecg_process(ecg, sampling_rate=1000)
eda_signals, eda_info = nk.eda_process(eda, sampling_rate=100)

# 2. Align sampling rates if needed
eda_signals_resampled = nk.signal_resample(eda_signals, sampling_rate=100,
                                           desired_sampling_rate=1000)

# 3. Merge signals into single DataFrame
signals = pd.concat([ecg_signals, eda_signals_resampled], axis=1)

# 4. Detect events
events = nk.events_find(trigger_channel, threshold=0.5)

# 5. Add event labels and conditions
event_labels = ['trial1', 'trial2', 'trial3', ...]
event_conditions = ['condition_A', 'condition_B', 'condition_A', ...]

# 6. Create epochs
epochs = nk.epochs_create(signals, events, sampling_rate=1000,
                          epochs_start=-1.0, epochs_end=5.0,
                          event_labels=event_labels,
                          event_conditions=event_conditions,
                          baseline_correction=True)

# 7. Signal-specific event-related analysis
ecg_results = nk.ecg_eventrelated(epochs)
eda_results = nk.eda_eventrelated(epochs)

# 8. Merge results
results = pd.merge(ecg_results, eda_results, left_index=True, right_index=True)

# 9. Statistical analysis by condition
results['Condition'] = event_conditions
condition_comparison = results.groupby('Condition').mean()
```

### 处理多种事件类型

```python
# Different event types with different markers
event_type1 = nk.events_find(trigger_ch1, threshold=0.5)
event_type2 = nk.events_find(trigger_ch2, threshold=0.5)

# Combine events with labels
all_events = np.concatenate([event_type1['onset'], event_type2['onset']])
event_labels = ['type1'] * len(event_type1['onset']) + ['type2'] * len(event_type2['onset'])

# Sort by time
sort_idx = np.argsort(all_events)
all_events = all_events[sort_idx]
event_labels = [event_labels[i] for i in sort_idx]

# Create epochs
epochs = nk.epochs_create(signals, all_events, sampling_rate=1000,
                          epochs_start=-0.5, epochs_end=3.0,
                          event_labels=event_labels)

# Separate by type
type1_epochs = {k: v for k, v in epochs.items() if v['Label'][0] == 'type1'}
type2_epochs = {k: v for k, v in epochs.items() if v['Label'][0] == 'type2'}
```

### 质量控制和伪品剔除

```python
# Remove epochs with excessive noise or artifacts
clean_epochs = {}
for key, epoch in epochs.items():
    # Example: reject if EDA amplitude too high (movement artifact)
    if epoch['EDA_Phasic'].abs().max() < 5.0:  # Threshold
        # Example: reject if heart rate change too large (invalid)
        if epoch['ECG_Rate'].max() - epoch['ECG_Rate'].min() < 50:
            clean_epochs[key] = epoch

print(f"Kept {len(clean_epochs)}/{len(epochs)} epochs")

# Analyze clean epochs
results = nk.ecg_eventrelated(clean_epochs)
```

## 统计考虑因素

### 样本量
- **ERP/平均**：每个条件最少 20-30 次以上试验
- **单独试验分析**：混合效应模型处理可变试验计数
- **组比较**：功率分析的试点数据

### 时间窗口选择
- **先验假设**：根据文献预先注册时间窗口
- **探索性**：使用完整纪元，正确进行多重比较
- **避免**：根据观察到的数据选择窗口（圆形）

### 基线期
- 应该没有预期效果
- 足够的持续时间进行稳定估计（典型值 500-1000 毫秒）
- 快速动态的时间较短（例如，惊吓：100 毫秒就足够了）

### 条件比较
- 受试者内设计的重复测量方差分析
- 不平衡数据的混合效应模型
- 非参数比较的排列检验
- 纠正多重比较（时间点/信号）

## 常见应用

**认知心理学：**
- P300 ERP分析
- 错误相关消极性（ERN）
- 注意力眨眼
- 工作记忆负荷效应

**情感神经科学：**
- 情绪图片查看（EDA、HR、面部肌电图）
- 恐惧调节（心率减速、SCR）
- 效价/唤醒维度

**临床研究：**
- 惊吓反应（眼轮匝肌肌电图）
- 定向响应（心率减速）
- 预期和预测错误

**心理生理学：**
- 心脏防御反应
- 定向反应与防御反应
- 情绪时呼吸的变化

**人机交互：**
- 活动期间的用户参与度
- 意外/违反预期
- 任务事件期间的认知负荷

## 参考文献

- 拉克，S.J. (2014)。事件相关电位技术简介（第二版）。麻省理工学院出版社。
- Bradley, M. M. 和 Lang, P. J. (2000)。测量情绪：行为、感觉和生理。载于 R. D. Lane 和 L. Nadel（主编），《情感认知神经科学》（第 242-276 页）。牛津大学出版社。
- Boucsein, W. (2012)。皮肤电活动（第二版）。施普林格。
- Gratton, G.、Coles, M. G. 和 Donchin, E. (1983)。一种离线去除眼部伪影的新方法脑电图和临床神经生理学，55(4), 468-484。