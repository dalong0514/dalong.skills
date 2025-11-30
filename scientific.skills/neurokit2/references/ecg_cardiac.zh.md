<!-- 此文件由机器翻译自 ecg_cardiac.md -->

# 心电图和心脏信号处理

## 概述

处理心电图 (ECG) 和光电体积描记法 (PPG) 信号以进行心血管分析。该模块提供用于 R 峰值检测、波形描绘、质量评估和心率分析的综合工具。

## 主要加工流程

### ecg_process()

完整的自动化心电图处理管道，可协调多个步骤。

```python
signals, info = nk.ecg_process(ecg_signal, sampling_rate=1000, method='neurokit')
```

**管道步骤：**
1.信号清洗（去噪）
2. R峰值检测
3.心率计算
4. 质量评估
5. QRS 轮廓（P、Q、S、T 波）
6. 心动相位测定

**退货：**
- `signals`：带有清理后的心电图、峰值、速率、质量、心动相位的 DataFrame
- `info`：具有 R 峰值位置和处理参数的字典

**常用方法：**
- `'neurokit'`（默认）：综合 NeuroKit2 管道
- `'biosppy'`：基于 BioSPPy 的处理
- `'pantompkins1985'`：泛汤普金斯算法
- `'hamilton2002'`、`'elgendi2010'`、`'engzeemod2012'`：替代方法

## 预处理函数

### ecg_clean()

使用特定于方法的滤波去除原始心电图信号中的噪声。

<<<代码块_1>>>

**方法：**
- `'neurokit'`：高通巴特沃斯滤波器 (0.5 Hz) + 电力线滤波
- `'biosppy'`：0.67-45 Hz 之间的 FIR 滤波
- `'pantompkins1985'`：带通 5-15 Hz + 基于导数
- `'hamilton2002'`：带通 8-16 Hz
- `'elgendi2010'`：带通 8-20 Hz
- `'engzeemod2012'`：FIR 带通 0.5-40 Hz

**关键参数：**
- `powerline`：消除 50 或 60 Hz 电力线噪声（默认值：50）

### ecg_peaks()

通过可选的伪影校正检测 ECG 信号中的 R 峰值。

<<<代码块_2>>>

**可用方法（13+ 算法）：**
- `'neurokit'`：针对可靠性进行优化的混合方法
- `'pantompkins1985'`：经典的 Pan-Tompkins 算法
- `'hamilton2002'`：汉密尔顿的自适应阈值
- `'christov2004'`：Christov 的自适应方法
- `'gamboa2008'`：甘博亚的方法
- `'elgendi2010'`：Elgendi 的两条移动平均线
- `'engzeemod2012'`：修改的 Engelse-Zeelenberg
- `'kalidas2017'`：基于 XQRS
- `'martinez2004'`、`'rodrigues2021'`、`'koka2022'`、`'promac'`：高级方法

**伪影校正：**
设置 `correct_artifacts=True` 以应用 Lipponen 和 Tarvainen (2019) 校正：
- 检测异位节拍、长/短间隔、错过的节拍
- 使用基于阈值的检测和可配置的参数

**退货：**
- 带有包含峰值索引的 `'ECG_R_Peaks'` 键的字典

### ecg_delineate()

识别 P、Q、S、T 波及其起始/偏移。

<<<代码块_3>>>

**方法：**
- `'dwt'`（默认）：基于离散小波变换的检测
- `'peak'`：R 峰值周围的简单峰值检测
- `'cwt'`：连续小波变换（Martinez 等人，2004）

**检测到的成分：**
- P 波：`ECG_P_Peaks`、`ECG_P_Onsets`、`ECG_P_Offsets`
- Q 波：`ECG_Q_Peaks`
- S 波：`ECG_S_Peaks`
- T 波：`ECG_T_Peaks`、`ECG_T_Onsets`、`ECG_T_Offsets`
- QRS 复合波：起始点和偏移点

**退货：**
- `waves`：包含所有波浪索引的字典
- `waves_peak`：具有峰值幅度的字典

### 心电图质量()

评估心电图信号完整性和质量。

<<<代码块_4>>>

**方法：**
- `'averageQRS'`（默认）：模板匹配相关性（Zhao & Zhu，2018）
  - 每次心跳返回质量分数 0-1
  - 阈值：>0.6 = 良好质量
- `'zhao2018'`：使用峰度、功率谱分布的多索引方法

**使用案例：**
- 识别低质量信号段
- 在分析之前过滤掉嘈杂的心跳
- 验证 R 峰值检测精度

## 分析函数

### ecg_analyze()

自动选择事件相关或间隔相关模式的高级分析。

<<<代码块_5>>>

**模式选择：**
- 持续时间 < 10 秒 → 事件相关分析
- 持续时间≥10秒→间隔相关分析

**退货：**
具有适合分析模式的心脏指标的数据帧。

### ecg_eventlated()

分析刺激锁定心电图时期的事件相关反应。

<<<代码块_6>>>

**计算指标：**
- `ECG_Rate_Baseline`：刺激前的平均心率
- `ECG_Rate_Min/Max`：纪元期间的最小/最大心率
- `ECG_Phase_Atrial/Ventricular`：刺激开始时的心脏阶段
- 跨纪元时间窗口的评级动态

**用例：**
具有离散试验的实验范式（例如刺激呈现、任务事件）。

### ecg_intervallated()

分析静息状态或长时间的连续心电图记录。

```python
results = nk.ecg_intervalrelated(signals, sampling_rate=1000)
```

**计算指标：**
- `ECG_Rate_Mean`：间隔内的平均心率
- 全面的 HRV 指标（委托给 `hrv()` 函数）
  - 时域：SDNN、RMSSD、pNN50 等
  - 频域：LF、HF、LF/HF 比
  - 非线性：庞加莱、熵、分形测度

**最短持续时间：**
- 基本费率：任何期限
- HRV 频率指标：建议 ≥60 秒，最佳 1-5 分钟

## 实用函数

### 心电图率()

根据 R 峰值间隔计算瞬时心率。

```python
heart_rate = nk.ecg_rate(peaks, sampling_rate=1000, desired_length=None)
```

**方法：**
- 计算连续 R 峰值之间的节拍间隔 (IBI)
- 转换为每分钟节拍 (BPM)：60 / IBI
- 如果指定了`desired_length`，则进行插值以匹配信号长度

**退货：**
- 瞬时心率值数组

### ecg_phase()

确定心房和心室收缩/舒张期。

```python
cardiac_phase = nk.ecg_phase(ecg_cleaned, rpeaks, delineate_info)
```

**计算的阶段：**
- `ECG_Phase_Atrial`：心房收缩 (1) 与舒张 (0)
- `ECG_Phase_Ventricular`：心室收缩 (1) 与舒张 (0)
- `ECG_Phase_Completion_Atrial/Ventricular`：阶段完成百分比 (0-1)

**用例：**
- 心脏锁定刺激呈现
- 心理生理学实验对心动周期事件进行计时

### ecg_segment()

提取个体心跳进行形态分析。

```python
heartbeats = nk.ecg_segment(ecg_cleaned, rpeaks, sampling_rate=1000)
```

**退货：**
- 纪元字典，每个纪元包含一次心跳
- 以 R 峰值为中心，具有可配置的前/后窗口
- 对于逐次心跳形态比较很有用

### ecg_invert()

自动检测并纠正反转的心电图信号。

```python
corrected_ecg, is_inverted = nk.ecg_invert(ecg_signal, sampling_rate=1000)
```

**方法：**
- 分析 QRS 复合极性
- 如果主要为负，则翻转信号
- 返回校正后的信号和反转状态

### ecg_rsp()

提取心电图衍生呼吸 (EDR) 作为呼吸代理信号。

```python
edr_signal = nk.ecg_rsp(ecg_cleaned, sampling_rate=1000, method='vangent2019')
```

**方法：**
- `'vangent2019'`：带通滤波 0.1-0.4 Hz
- `'charlton2016'`：带通 0.15-0.4 Hz
- `'soni2019'`：带通 0.08-0.5 Hz

**用例：**
- 当直接呼吸信号不可用时估计呼吸
- 多模态生理分析

## 模拟和可视化

### ecg_simulate()

生成合成心电图信号以进行测试和验证。

```python
synthetic_ecg = nk.ecg_simulate(duration=10, sampling_rate=1000, heart_rate=70, method='ecgsyn', noise=0.01)
```

**方法：**
- `'ecgsyn'`：现实动力学模型（McSharry 等人，2003）
  - 模拟 P-QRS-T 复合波形态
  - 生理上合理的波形
- `'simple'`：更快的基于小波的近似
  - 类高斯 QRS 复合波
  - 不太现实但计算效率高

**关键参数：**
- `heart_rate`：平均 BPM（默认值：70）
- `heart_rate_std`：心率变异幅度（默认值：1）
- `noise`：高斯噪声级别（默认值：0.01）
- `random_state`：可重复性的种子

### ecg_plot()

可视化处理后的心电图以及检测到的 R 峰值和信号质量。

```python
nk.ecg_plot(signals, info)
```

**显示：**
- 原始和清理后的心电图信号
- 检测到的 R 峰重叠
- 心率追踪
- 信号质量指标

## ECG 特定注意事项

### 采样率建议
- **最小值**：基本 R 峰值检测为 250 Hz
- **推荐**：500-1000 Hz 用于波形描绘
- **高分辨率**：2000+ Hz 用于详细形态分析

### 信号持续时间要求
- **R 峰值检测**：任何持续时间（至少 2 次心跳）
- **基础心率**：≥10秒
- **HRV时域**：≥60秒
- **HRV 频域**：1-5 分钟（最佳）
- **超低频HRV**：≥24小时

### 常见问题和解决方案

**R 峰检测不良：**
- 尝试不同的方法：`method='pantompkins1985'` 通常很稳健
- 确保足够的采样率（≥250 Hz）
- 检查心电图倒置：使用 `ecg_invert()`
- 应用伪影校正：`correct_artifacts=True`

**噪音信号：**
- 根据噪音类型使用适当的清洁方法
- 如果在美国/欧洲之外，请调整电源线频率
- 在分析之前考虑信号质量评估

**缺少波形分量：**
- 提高采样率（≥500 Hz 用于描绘）
- 尝试不同的描绘方法（`'dwt'`、`'peak'`、`'cwt'`）
- 使用`ecg_quality()`验证信号质量

## 与其他信号集成
### ECG + RSP：呼吸性窦性心律失常
```python
# Process both signals
ecg_signals, ecg_info = nk.ecg_process(ecg, sampling_rate=1000)
rsp_signals, rsp_info = nk.rsp_process(rsp, sampling_rate=1000)

# Compute RSA
rsa = nk.hrv_rsa(ecg_info['ECG_R_Peaks'], rsp_signals['RSP_Clean'], sampling_rate=1000)
```

### 多模式集成
```python
# Process multiple signals at once
bio_signals, bio_info = nk.bio_process(
    ecg=ecg_signal,
    rsp=rsp_signal,
    eda=eda_signal,
    sampling_rate=1000
)
```

## 参考文献

- Pan, J. 和 Tompkins, W. J. (1985)。实时 QRS 检测算法。 IEEE 生物医学工程汇刊，32(3), 230-236。
- 汉密尔顿，P. (2002)。开源心电图分析。心脏病学计算机，101-104。
- Martinez, J. P.、Almeida, R.、Olmos, S.、Rocha, A. P. 和 Laguna, P. (2004)。基于小波的心电图轮廓仪：标准数据库的评估。 IEEE 生物医学工程汇刊，51(4), 570-581。
- Lipponen, J. A. 和 Tarvainen, M. P. (2019)。使用新颖的心跳分类进行心率变异性时间序列伪影校正的稳健算法。医学工程与技术杂志，43(3), 173-181。