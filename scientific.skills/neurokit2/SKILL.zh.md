<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：神经套件2
描述：用于分析生理数据的综合生物信号处理工具包，包括 ECG、EEG、EDA、RSP、PPG、EMG 和 EOG 信号。在处理心血管信号、大脑活动、皮肤电反应、呼吸模式、肌肉活动或眼球运动时使用此技能。适用于心率变异性分析、事件相关电位、复杂性测量、自主神经系统评估、心理生理学研究和多模态生理信号集成。
---

# 神经套件2

## 概述

NeuroKit2 是一个用于处理和分析生理信号（生物信号）的综合 Python 工具包。使用此技能处理心血管、神经、自主、呼吸和肌肉信号，以进行心理生理学研究、临床应用和人机交互研究。

## 何时使用此技能

在处理以下问题时应用此技能：
- **心脏信号**：ECG、PPG、心率变异性 (HRV)、脉搏分析
- **脑信号**：EEG 频段、微观状态、复杂性、源定位
- **自主信号**：皮肤电活动 (EDA/GSR)、皮肤电导反应 (SCR)
- **呼吸信号**：呼吸频率、呼吸变异性 (RRV)、每次呼吸量
- **肌肉信号**：肌电波幅、肌肉激活检测
- **眼动追踪**：EOG、眨眼检测和分析
- **多模态集成**：同时处理多个生理信号
- **复杂性分析**：熵度量、分形维数、非线性动力学

## 核心能力

### 1. 心脏信号处理 (ECG/PPG)

处理心电图和光电体积描记信号以进行心血管分析。有关详细工作流程，请参阅 `references/ecg_cardiac.md`。

**主要工作流程：**
- 心电图处理流程：清洗→R峰检测→勾画→质量评估
- 跨时间、频率和非线性域的 HRV 分析
- PPG脉冲分析和质量评估
- 心电图衍生的呼吸提取

**主要功能：**
```python
import neurokit2 as nk

# Complete ECG processing pipeline
signals, info = nk.ecg_process(ecg_signal, sampling_rate=1000)

# Analyze ECG data (event-related or interval-related)
analysis = nk.ecg_analyze(signals, sampling_rate=1000)

# Comprehensive HRV analysis
hrv = nk.hrv(peaks, sampling_rate=1000)  # Time, frequency, nonlinear domains
```

### 2.心率变异性分析

根据心脏信号计算综合 HRV 指标。有关所有索引和特定于域的分析，请参阅 `references/hrv.md`。

**支持的域：**
- **时域**：SDNN、RMSSD、pNN50、SDSD 和派生指标
- **频域**：ULF、VLF、LF、HF、VHF 功率和比率
- **非线性域**：庞加莱图 (SD1/SD2)、熵度量、分形维数
- **专业**：呼吸性窦性心律失常（RSA）、复发定量分析（RQA）

**主要功能：**
<<<代码块_1>>>

### 3. 脑信号分析 (EEG)

分析脑电图信号的频率功率、复杂性和微观状态模式。请参阅 `references/eeg.md` 了解详细的工作流程和跨国公司集成。

**主要能力：**
- 频段功率分析（Delta、Theta、Alpha、Beta、Gamma）
- 渠道质量评估和重新参考
- 来源本地化（sLORETA、MNE）
- 微观状态分割和过渡动力学
- 全球场功率和相异性测量

**主要功能：**
<<<代码块_2>>>

### 4. 皮肤电活动 (EDA)

处理皮肤电导信号以进行自主神经系统评估。有关详细工作流程，请参阅 `references/eda.md`。

**主要工作流程：**
- 信号分解为主音和相位分量
- 皮肤电导反应（SCR）检测和分析
- 交感神经系统指数计算
- 自相关和变化点检测

**主要功能：**
<<<代码块_3>>>

### 5.呼吸信号处理（RSP）

分析呼吸模式和呼吸变异性。有关详细工作流程，请参阅`references/rsp.md`。

**主要能力：**
- 呼吸频率计算和变异性分析
- 呼吸幅度和对称性评估
- 每次呼吸量（fMRI 应用）
- 呼吸幅度变异性（RAV）

**主要功能：**
<<<代码块_4>>>

### 6.肌电图（EMG）

处理肌肉活动信号以进行激活检测和幅度分析。请参阅 `references/emg.md` 了解工作流程。

**主要功能：**
<<<代码块_5>>>

### 7. 眼电图 (EOG)

分析眼球运动和眨眼模式。请参阅 `references/eog.md` 了解工作流程。

**主要功能：**
<<<代码块_6>>>

### 8. 一般信号处理

对任何信号应用滤波、分解和变换操作。有关综合实用程序，请参阅 `references/signal_processing.md`。

**关键操作：**
- 滤波（低通、高通、带通、带阻）
- 分解（EMD、SSA、小波）
- 峰值检测和校正
- 功率谱密度估计
- 信号插值和重采样
- 自相关和同步分析

**主要功能：**
```python
# Filtering
filtered = nk.signal_filter(signal, sampling_rate=1000, lowcut=0.5, highcut=40)

# Peak detection
peaks = nk.signal_findpeaks(signal)

# Power spectral density
psd = nk.signal_psd(signal, sampling_rate=1000)
```

### 9. 复杂性和熵分析

计算非线性动力学、分形维数和信息论测量。请参阅 `references/complexity.md` 了解所有可用指标。

**可用措施：**
- **熵**：香农、近似、样本、排列、谱、模糊、多尺度
- **分形维数**：Katz、Higuchi、Petrosian、Sevcik、相关维数
- **非线性动力学**：Lyapunov 指数、Lempel-Ziv 复杂度、递归量化
- **DFA**：去趋势波动分析、多重分形 DFA
- **信息论**：费舍尔信息、互信息

**主要功能：**
```python
# Multiple complexity metrics at once
complexity_indices = nk.complexity(signal, sampling_rate=1000)

# Specific measures
apen = nk.entropy_approximate(signal)
dfa = nk.fractal_dfa(signal)
lyap = nk.complexity_lyapunov(signal, sampling_rate=1000)
```

### 10. 事件相关分析

围绕刺激事件创建纪元并分析生理反应。请参阅 `references/epochs_events.md` 了解工作流程。

**主要能力：**
- 从事件标记创建纪元
- 与事件相关的平均和可视化
- 基线校正选项
- 带置信区间的总平均计算

**主要功能：**
```python
# Find events in signal
events = nk.events_find(trigger_signal, threshold=0.5)

# Create epochs around events
epochs = nk.epochs_create(signals, events, sampling_rate=1000,
                          epochs_start=-0.5, epochs_end=2.0)

# Average across epochs
grand_average = nk.epochs_average(epochs)
```

### 11. 多信号集成

同时处理多个生理信号并统一输出。请参阅 `references/bio_module.md` 了解集成工作流程。

**主要功能：**
```python
# Process multiple signals at once
bio_signals, bio_info = nk.bio_process(
    ecg=ecg_signal,
    rsp=rsp_signal,
    eda=eda_signal,
    emg=emg_signal,
    sampling_rate=1000
)

# Analyze all processed signals
bio_analysis = nk.bio_analyze(bio_signals, sampling_rate=1000)
```

## 分析模式

NeuroKit2 根据数据持续时间自动在两种分析模式之间进行选择：

**事件相关分析**（< 10 秒）：
- 分析刺激锁定反应
- 基于纪元的分割
- 适用于离散试验的实验范式

**间隔相关分析**（≥ 10 秒）：
- 描述长时间内的生理模式特征
- 静息状态或持续活动
- 适用于基线测量和长期监测

大多数 `*_analyze()` 函数会自动选择适当的模式。

## 安装

```bash
uv pip install neurokit2
```

对于开发版本：
```bash
uv pip install https://github.com/neuropsychology/NeuroKit/zipball/dev
```

## 常见工作流程

### 快速入门：心电图分析
```python
import neurokit2 as nk

# Load example data
ecg = nk.ecg_simulate(duration=60, sampling_rate=1000)

# Process ECG
signals, info = nk.ecg_process(ecg, sampling_rate=1000)

# Analyze HRV
hrv = nk.hrv(info['ECG_R_Peaks'], sampling_rate=1000)

# Visualize
nk.ecg_plot(signals, info)
```

### 多模态分析
```python
# Process multiple signals
bio_signals, bio_info = nk.bio_process(
    ecg=ecg_signal,
    rsp=rsp_signal,
    eda=eda_signal,
    sampling_rate=1000
)

# Analyze all signals
results = nk.bio_analyze(bio_signals, sampling_rate=1000)
```

### 事件相关潜力
```python
# Find events
events = nk.events_find(trigger_channel, threshold=0.5)

# Create epochs
epochs = nk.epochs_create(processed_signals, events,
                          sampling_rate=1000,
                          epochs_start=-0.5, epochs_end=2.0)

# Event-related analysis for each signal type
ecg_epochs = nk.ecg_eventrelated(epochs)
eda_epochs = nk.eda_eventrelated(epochs)
```

## 参考文献

该技能包括按信号类型和分析方法组织的综合参考文档：

- **ecg_cardiac.md**：ECG/PPG 处理、R 峰值检测、描绘、质量评估
- **hrv.md**：所有领域的心率变异性指数
- **eeg.md**：脑电图分析、频段、微观状态、源定位
- **eda.md**：皮肤电活动处理和 SCR 分析
- **rsp.md**：呼吸信号处理和变异性
- **ppg.md**：光电体积描记法信号分析
- **emg.md**：肌电图处理和激活检测
- **eog.md**：眼电图和眨眼分析
- **signal_processing.md**：通用信号实用程序和转换
- **complexity.md**：熵、分形和非线性度量
- **epochs_events.md**：事件相关分析和纪元创建
- **bio_module.md**：多信号集成工作流程

使用读取工具根据需要加载特定参考文件以访问详细的功能文档和参数。

## 其他资源

- 官方文档：https://neuropsychology.github.io/NeuroKit/
- GitHub 存储库：https://github.com/neuropsychology/NeuroKit
- 出版物：Makowski 等人。 （2021）。 NeuroKit2：用于神经生理信号处理的 Python 工具箱。行为研究方法。 https://doi.org/10.3758/s13428-020-01516-y