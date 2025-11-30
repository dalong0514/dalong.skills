<!-- 此文件由机器翻译自 eda.md -->

# 皮肤电活动 (EDA) 分析

## 概述

皮肤电活动 (EDA)，也称为皮肤电反应 (GSR) 或皮肤电导 (SC)，测量皮肤的电导，反映交感神经系统的唤醒和汗腺活动。 EDA 广泛应用于心理生理学、情感计算和测谎。

## 主要加工流程

### eda_process()

自动处理原始 EDA 信号，返回主音/相位分解和 SCR 特征。

```python
signals, info = nk.eda_process(eda_signal, sampling_rate=100, method='neurokit')
```

**管道步骤：**
1.信号净化（低通滤波）
2. 强直相分解
3.皮肤电导反应（SCR）检测
4. SCR特征提取（起始、峰值、幅度、上升/恢复时间）

**退货：**
- `signals`：数据帧具有：
  - `EDA_Clean`：过滤后的信号
  - `EDA_Tonic`：缓慢变化的基线
  - `EDA_Phasic`：快速变化的响应
  - `SCR_Onsets`、`SCR_Peaks`、`SCR_Height`：响应标记
  - `SCR_Amplitude`、`SCR_RiseTime`、`SCR_RecoveryTime`：响应特征
- `info`：带有处理参数的字典

**方法：**
- `'neurokit'`：cvxEDA 分解 + Neurokit 峰值检测
- `'biosppy'`：中值平滑 + biosppy 方法

## 预处理函数

### eda_clean()

通过低通滤波去除噪声。

<<<代码块_1>>>

**方法：**
- `'neurokit'`：低通巴特沃斯滤波器（3 Hz 截止）
- `'biosppy'`：低通巴特沃斯滤波器（5 Hz 截止）

**自动跳过：**
- 如果采样率 < 7 Hz，则跳过清洁（已经是低通）

**理由：**
- EDA 频率内容通常为 0-3 Hz
- 消除高频噪声和运动伪影
- 保留慢速 SCR（典型上升时间 1-3 秒）

### eda_phasic()

将 EDA 分解为强直（慢基线）和阶段（快速反应）组件。

<<<代码块_2>>>

**方法：**

**1. cvxEDA（默认，推荐）：**
<<<代码块_3>>>
- 凸优化方法（Greco 等人，2016）
- 稀疏相位驱动模型
- 最准确的生理学
- 计算密集但卓越的分解

**2.中值平滑：**
<<<代码块_4>>>
- 具有可配置窗口的中值滤波器
- 快速、简单
- 不如 cvxEDA 准确

**3.高通滤波（Biopac 的 Acqknowledge）：**
<<<代码块_5>>>
- 高通滤波器 (0.05 Hz) 提取相位
- 快速计算
- 由减法得出的补药

**4. SparsEDA：**
<<<代码块_6>>>
- 稀疏反卷积方法
- 替代优化方法

**退货：**
- `tonic`：缓慢变化的皮肤电导水平 (SCL)
- `phasic`：快速皮肤电导响应 (SCR)

**生理解释：**
- **补品 (SCL)**：基线唤醒、一般激活、水合作用
- **阶段性 (SCR)**：与事件相关的反应、定向、情绪反应

### eda_peaks()

检测相位分量中的皮肤电导反应 (SCR)。

```python
peaks, info = nk.eda_peaks(eda_phasic, sampling_rate=100, method='neurokit',
                           amplitude_min=0.1)
```

**方法：**
- `'neurokit'`：针对可靠性进行了优化，可配置阈值
- `'gamboa2008'`：Gamboa 算法
- `'kim2004'`：Kim 的方法
- `'vanhalem2020'`：Van Halem 方法
- `'nabian2018'`：Nabian 算法

**关键参数：**
- `amplitude_min`：最小 SCR 幅度（默认值：0.1 µS）
  - 太低：噪声导致误报
  - 太高：错过小但有效的回应
- `rise_time_max`：最大上升时间（默认值：2 秒）
- `rise_time_min`：最小上升时间（默认值：0.01 秒）

**退货：**
- 字典：
  - `SCR_Onsets`：SCR 开始的索引
  - `SCR_Peaks`：峰值幅度索引
  - `SCR_Height`：基线上方的峰值高度
  - `SCR_Amplitude`：起始峰值幅度
  - `SCR_RiseTime`：开始到峰值的持续时间
  - `SCR_RecoveryTime`：峰值到恢复持续时间（50% 衰减）

**SCR 时序约定：**
- **延迟**：刺激后 1-3 秒（典型）
- **上升时间**：0.5-3秒
- **恢复时间**：2-10 秒（恢复 50%）
- **最小振幅**：0.01-0.05 µS（检测阈值）

### eda_fixpeaks()

正确检测到的 SCR 峰值（当前为 EDA 的占位符）。

```python
corrected_peaks = nk.eda_fixpeaks(peaks)
```

**注：** 由于动态较慢，对 EDA 的重要性不如心脏信号。

## 分析函数

### eda_analyze()

根据数据持续时间自动选择适当的分析类型。

```python
analysis = nk.eda_analyze(signals, sampling_rate=100)
```

**模式选择：**
- 持续时间 < 10 秒 → `eda_eventrelated()`
- 持续时间 ≥ 10 秒 → `eda_intervalrelated()`

**退货：**
- 具有适合分析模式的 EDA 指标的 DataFrame

### eda_eventlated()

分析刺激锁定 EDA 时期的事件相关响应。

```python
results = nk.eda_eventrelated(epochs)
```

**计算指标（每个时期）：**
- `EDA_SCR`：存在 SCR（二进制：0 或 1）
- `SCR_Amplitude`：历元期间的最大 SCR 幅度
- `SCR_Magnitude`：平均阶段性活动
- `SCR_Peak_Amplitude`：起始到峰值幅度
- `SCR_RiseTime`：从发病到达到峰值的时间
- `SCR_RecoveryTime`：恢复 50% 的时间
- `SCR_Latency`：从刺激到 SCR 发作的延迟
- `EDA_Tonic`：时期内的平均补品水平

**典型参数：**
- 纪元持续时间：刺激后 0-10 秒
- 基线：刺激前 -1 至 0 秒
- 预期 SCR 延迟：1-3 秒

**使用案例：**
- 情绪刺激处理（图像、声音）
- 认知负荷评估（心算）
- 预期和预测错误
- 定向反应

### eda_intervallated()

分析扩展 EDA 记录的整体唤醒和激活模式。

```python
results = nk.eda_intervalrelated(signals, sampling_rate=100)
```

**计算指标：**
- `SCR_Peaks_N`：检测到的 SCR 数量
- `SCR_Peaks_Amplitude_Mean`：平均 SCR 幅度
- `EDA_Tonic_Mean`、`EDA_Tonic_SD`：补品级别统计
- `EDA_Sympathetic`：交感神经系统指数
- `EDA_SympatheticN`：标准化交感神经指数
- `EDA_Autocorrelation`：时间结构（滞后 4 秒）
- `EDA_Phasic_*`：相位分量的平均值、SD、最小值、最大值

**录音时长：**
- **最短**：10 秒
- **推荐**：稳定 SCR 速率需要 60 秒以上
- **交感指数**：需要≥64秒

**使用案例：**
- 静息状态唤醒评估
- 压力水平监测
- 基线交感神经活动
- 长期情感状态

## 专业分析功能

### eda_sympathetic()

从频段 (0.045-0.25 Hz) 得出交感神经系统活动。

```python
sympathetic = nk.eda_sympathetic(signals, sampling_rate=100, method='posada',
                                  show=False)
```

**方法：**
- `'posada'`：Posada-Quintero 方法 (2016)
  - 0.045-0.25 Hz 频段的光谱功率
  - 针对其他自主措施进行验证
- `'ghiasi'`：Ghiasi 方法 (2018)
  - 基于频率的替代方法

**要求：**
- **最短持续时间**：64 秒
- 足够目标频段的频率分辨率

**退货：**
- `EDA_Sympathetic`：交感指数（绝对）
- `EDA_SympatheticN`：标准化交感神经指数 (0-1)

**释义：**
- 较高的值：增加交感神经兴奋
- 反映强直交感神经活动，而不是阶段性反应
- 补充 SCR 分析

**使用案例：**
- 压力评估
- 随时间的唤醒监测
- 认知负荷测量
- 补充 HRV 以实现自主平衡

### eda_autocor()

计算自相关以评估 EDA 信号的时间结构。

```python
autocorr = nk.eda_autocor(eda_phasic, sampling_rate=100, lag=4)
```

**参数：**
- `lag`：时间延迟（以秒为单位）（默认值：4 秒）

**释义：**
- 高自相关性：持续、缓慢变化的信号
- 低自相关：快速、不相关的波动
- 反映 SCR 的时间规律性

**用例：**
- 评估信号质量
- 表征反应模式
- 区分持续性唤醒与短暂性唤醒

### eda_changepoints()

检测 EDA 信号均值和方差的突变。

```python
changepoints = nk.eda_changepoints(eda_phasic, penalty=10000, show=False)
```

**方法：**
- 基于惩罚的分割
- 识别状态之间的转换

**参数：**
- `penalty`：控制灵敏度（默认值：10,000）
  - 更高的惩罚：更少、更稳健的变化点
  - 较低的惩罚：对微小变化更敏感

**退货：**
- 检测到的变化点的索引
- 可选的分段可视化

**使用案例：**
- 识别连续监控中的状态转换
- 按唤醒水平细分数据
- 检测实验中的相变
- 自动纪元定义

## 可视化

### eda_plot()

创建已处理 EDA 的静态或交互式可视化。

```python
nk.eda_plot(signals, info, static=True)
```

**显示：**
- 原始和清理过的 EDA 信号
- 补品和阶段成分
- 检测到 SCR 起始、峰值和恢复
- 交感神经指数时程（如果计算的话）

**交互模式（`static=False`）：**
- 基于情节的交互式探索
- 缩放、平移、悬停以查看详细信息
- 导出为图像格式

## 模拟与测试

### eda_simulate()

生成具有可配置参数的合成 EDA 信号。

```python
synthetic_eda = nk.eda_simulate(duration=10, sampling_rate=100, scr_number=3,
                                noise=0.01, drift=0.01)
```

**参数：**
- `duration`：信号长度（以秒为单位）
- `sampling_rate`：采样频率 (Hz)
- `scr_number`：要包含的 SCR 数量
- `noise`：高斯噪声级别
- `drift`：慢基线漂移幅度
- `random_state`：可重复性的种子

**退货：**
- 具有真实 SCR 形态的合成 EDA 信号

**使用案例：**
- 算法测试和验证
- 教育示范
- 方法比较

## 实际考虑

### 采样率建议
- **最小值**：10 Hz（足以满足慢速 SCR）
- **标准**：20-100 Hz（大多数商业系统）
- **高分辨率**：1000 Hz（研究级，过采样）

### 录音时长
- **SCR 检测**：≥10 秒（取决于刺激）
- **与事件相关**：每次试验通常需要 10-20 秒
- **与间隔相关**：稳定估计≥60 秒
- **交感指数**：≥64秒（频率分辨率）

### 电极放置
- **标准站点**：
  - 手掌：远端/中指骨（手指）
  足底：脚底
- **高密度**：鱼际/小鱼际隆起
- **避免**：多毛皮肤、汗腺密度低的区域
- **双边**：左手与右手（通常相似）

### 信号质量问题

**平坦信号（无变化）：**
- 检查电极接触和凝胶
- 验证在汗腺丰富的区域是否正确放置
- 允许 5-10 分钟的适应期

**噪音过大：**
- 运动伪影：最大限度地减少参与者的运动
- 电气干扰：检查接地、屏蔽
- 热效应：控制室温

**基线漂移：**
- 正常：几分钟内缓慢变化
- 过多：电极极化、接触不良
- 解决方案：使用`eda_phasic()`来分离主音漂移

**无回应者：**
- ~5-10% 的人口具有最低限度的 EDA
- 遗传/生理变异
- 不表示设备故障

### 最佳实践

**预处理工作流程：**
```python
# 1. Clean signal
cleaned = nk.eda_clean(eda_raw, sampling_rate=100, method='neurokit')

# 2. Decompose tonic/phasic
tonic, phasic = nk.eda_phasic(cleaned, sampling_rate=100, method='cvxeda')

# 3. Detect SCRs
signals, info = nk.eda_peaks(phasic, sampling_rate=100, amplitude_min=0.05)

# 4. Analyze
analysis = nk.eda_analyze(signals, sampling_rate=100)
```

**事件相关工作流程：**
```python
# 1. Process signal
signals, info = nk.eda_process(eda_raw, sampling_rate=100)

# 2. Find events
events = nk.events_find(trigger_channel, threshold=0.5)

# 3. Create epochs (-1 to 10 seconds around stimulus)
epochs = nk.epochs_create(signals, events, sampling_rate=100,
                          epochs_start=-1, epochs_end=10)

# 4. Event-related analysis
results = nk.eda_eventrelated(epochs)

# 5. Statistical analysis
# Compare SCR amplitude across conditions
```

## 临床和研究应用

**情感和情感科学：**
- 情绪的唤醒维度（不是效价）
- 情感图片查看
- 音乐引发的情绪
- 恐惧调节

**认知过程：**
- 脑力负荷和努力
- 注意和警惕
- 决策和不确定性
- 错误处理

**临床人群：**
- 焦虑症：基线升高、反应过度
- PTSD：恐惧调节、消退缺陷
- 自闭症：非典型的唤醒模式
- 精神病态：恐惧反应减少

**应用的设置：**
- 测谎（测谎仪）
- 用户体验研究
- 驾驶员监控
- 现实环境中的压力评估

**神经影像整合：**
- fMRI：EDA 与杏仁核、岛叶活动相关
- 脑成像期间的并发记录
- 自主大脑耦合

## 解释指南

**可控硅振幅：**
- **0.01-0.05 µS**：小但可检测
- **0.05-0.2 µS**：中等响应
- **>0.2 µS**：大响应
- **上下文相关**：主体内标准化

**可控硅频率：**
- **休息**：每分钟 1-3 次 SCR（典型）
- **有压力**：每分钟 >5 次 SCR
- **非特异性 SCR**：自发的（无可识别的刺激）

**补品SCL：**
- **范围**：2-20 µS（个体间差异很大）
- **受试者内变化**比绝对水平更容易解释
- **增加**：觉醒、压力、认知负荷
- **减少**：放松、习惯

## 参考文献

- Boucsein, W. (2012)。皮肤电活动（第二版）。施普林格科学与商业媒体。
- Greco, A.、Valenza, G. 和 Scilingo, E. P. (2016)。 cvxEDA：皮肤电活动处理的凸优化方法。 IEEE 生物医学工程汇刊，63(4), 797-804。
- Posada-Quintero, H. F.、Florian, J. P.、Orjuela-Cañón, A. D.、Aljama-Corrales, T.、Charleston-Villalobos, S. 和 Chon, K. H. (2016)。用于评估交感神经功能的皮肤电活动的功率谱密度分析。生物医学工程年鉴，44(10), 3124-3135。
- Dawson, M. E.、Schell, A. M. 和 Filion, D. L. (2017)。皮肤电系统。见《心理生理学手册》（第 217-243 页）。剑桥大学出版社。