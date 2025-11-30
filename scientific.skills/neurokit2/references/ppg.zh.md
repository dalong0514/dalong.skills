<!-- 此文件由机器翻译自 ppg.md -->

# 光电体积描记法 (PPG) 分析

## 概述

光电体积描记法 (PPG) 使用光学传感器测量微血管组织中的血容量变化。 PPG 广泛应用于可穿戴设备、脉搏血氧计以及用于心率、脉搏特征和心血管评估的临床监测仪。

## 主要加工流程

### ppg_process()

自动化 PPG 信号处理管道。

```python
signals, info = nk.ppg_process(ppg_signal, sampling_rate=100, method='elgendi')
```

**管道步骤：**
1.信号清洗（过滤）
2. 收缩期峰值检测
3.心率计算
4. 信号质量评估

**退货：**
- `signals`：数据帧具有：
  - `PPG_Clean`：过滤后的 PPG 信号
  - `PPG_Peaks`：收缩期峰值标记
  - `PPG_Rate`：瞬时心率 (BPM)
  - `PPG_Quality`：信号质量指示器
- `info`：具有峰值索引和参数的字典

**方法：**
- `'elgendi'`：Elgendi 等人。 (2013) 算法（默认，稳健）
- `'nabian2018'`：Nabian 等人。 （2018）方法

## 预处理函数

### ppg_clean()

准备原始 PPG 信号以进行峰值检测。

<<<代码块_1>>>

**方法：**

**1. Elgendi（默认）：**
- 巴特沃斯带通滤波器（0.5-8 Hz）
- 消除基线漂移和高频噪声
- 针对峰值检测可靠性进行了优化

**2.纳比安2018：**
- 替代过滤方法
- 不同的频率特性

**PPG信号特征：**
- **收缩峰值**：快速上冲，尖锐峰值（心脏射血）
- **重搏切迹**：次要峰（主动脉瓣关闭）
- **基线**：由于呼吸、运动、灌注而缓慢漂移

### ppg_peaks()

检测 PPG 信号的收缩峰值。

<<<代码块_2>>>

**方法：**
- `'elgendi'`：具有动态阈值的两条移动平均线
- `'bishop'`：Bishop 算法
- `'nabian2018'`：Nabian 的方法
- `'scipy'`：简单的 scipy 峰值检测

**伪影校正：**
- 设置`correct_artifacts=True`进行生理合理性检查
- 根据节拍间间隔异常值去除虚假峰值

**退货：**
- 带有包含峰值索引的 `'PPG_Peaks'` 键的字典

**典型的节拍间隔：**
- 成人静息：600-1200 毫秒（50-100 BPM）
- 运动员：可以更长（心动过缓）
- 压力/锻炼：较短（<600 ms，>100 BPM）

### ppg_findpeaks()

通过算法比较进行低电平峰值检测。

<<<代码块_3>>>

**用例：**
- 自定义参数调整
- 算法测试
- 研究方法开发

## 分析函数

### ppg_analyze()

自动选择事件相关或间隔相关分析。

<<<代码块_4>>>

**模式选择：**
- 持续时间 < 10 秒 → 与事件相关
- 持续时间 ≥ 10 秒 → 与间隔相关

### ppg_eventlated()

分析 PPG 对离散事件/刺激的反应。

<<<代码块_5>>>

**计算指标（每个时期）：**
- `PPG_Rate_Baseline`：事件前的心率
- `PPG_Rate_Min/Max`：纪元期间的最小/最大心率
- 跨纪元时间窗口的评级动态

**使用案例：**
- 心血管对情绪刺激的反应
- 认知负荷评估
- 应激反应范例

### ppg_intervallated()

分析扩展 PPG 录音。

<<<代码块_6>>>

**计算指标：**
- `PPG_Rate_Mean`：平均心率
- 心率变异性 (HRV) 指标
  - 委托给 `hrv()` 函数
  - 时域、频域和非线性域

**录音时长：**
- 最短：基本费率 60 秒
- HRV 分析：建议 2-5 分钟

**使用案例：**
- 静息态心血管评估
- 可穿戴设备数据分析
- 长期心率监测

## 质量评估

### ppg_quality()

评估信号质量和可靠性。

```python
quality = nk.ppg_quality(ppg_signal, sampling_rate=100, method='averageQRS')
```

**方法：**

**1.平均QRS（默认）：**
- 模板匹配方法
- 将每个脉冲与平均模板相关联
- 返回每个节拍的质量分数 0-1
- 阈值：>0.6 = 可接受的质量

**2.差异：**
- 地形差异测量
- 检测形态变化

**使用案例：**
- 识别损坏的段
- 在分析前过滤低质量数据
- 验证峰值检测精度

**常见质量问题：**
- 运动伪影：信号突然变化
- 传感器接触不良：振幅低、噪音大
- 血管收缩：信号幅度降低（寒冷、压力）

## 实用函数

### ppg_segment()

提取单个脉冲进行形态分析。

```python
pulses = nk.ppg_segment(cleaned_ppg, peaks, sampling_rate=100)
```

**退货：**
- 脉搏时期词典，每个时期都以收缩期峰值为中心
- 启用脉冲间比较
- 跨条件的形态分析

**使用案例：**
- 脉搏波分析
- 动脉僵硬度代理
- 血管老化评估

### ppg_methods()

记录分析中使用的预处理方法。

```python
methods_info = nk.ppg_methods(method='elgendi')
```

**退货：**
- 记录处理管道的字符串
- 适用于出版物中的方法部分

## 模拟和可视化

### ppg_simulate()

生成用于测试的合成 PPG 信号。

```python
synthetic_ppg = nk.ppg_simulate(duration=60, sampling_rate=100, heart_rate=70,
                                noise=0.1, random_state=42)
```

**参数：**
- `heart_rate`：平均 BPM（默认值：70）
- `heart_rate_std`：HRV 幅度
- `noise`：高斯噪声级别
- `random_state`：再现性种子

**使用案例：**
- 算法验证
- 参数优化
- 教育示范

### ppg_plot()

可视化处理后的 PPG 信号。

```python
nk.ppg_plot(signals, info, static=True)
```

**显示：**
- 原始和清洁的 PPG 信号
- 检测到收缩期峰值
- 瞬时心率轨迹
- 信号质量指标

## 实际考虑

### 采样率建议
- **最小值**：20 Hz（基本心率）
- **标准**：50-100 Hz（大多数可穿戴设备）
- **高分辨率**：200-500 Hz（研究、脉搏波分析）
- **过多**：>1000 Hz（PPG 不需要）

### 录音时长
- **心率**：≥10 秒（几次心跳）
- **HRV 分析**：最少 2-5 分钟
- **长期监控**：几小时到几天（可穿戴设备）

### 传感器放置

**常用站点：**
- **指尖**：最高信号质量，最常见
- **耳垂**：运动伪影较少，临床使用
- **手腕**：可穿戴设备（智能手表）
- **额头**：反射模式，医疗监测

**透射率与反射率：**
- **透射率**：光穿过组织（指尖、耳垂）
  - 更高的信号质量
  - 减少运动伪影
- **反射率**：从组织（手腕、前额）反射的光
  - 更容易受到噪音影响
  - 方便可穿戴设备

### 常见问题和解决方案

**低信号幅度：**
- 灌注不良：温暖双手，增加血液流动
- 传感器接触：调整位置，清洁皮肤
- 血管收缩：环境温度、压力

**运动伪影：**
- 可穿戴设备的主要问题
- 自适应滤波、基于加速度计的校正
- 模板匹配、异常值拒绝

**基线漂移：**
- 呼吸调节（正常）
- 运动或压力变化
- 高通滤波或去趋势

**缺失的峰：**
- 信号质量低：检查传感器接触情况
- 算法参数：调整阈值
- 尝试替代检测方法

### 最佳实践

**标准工作流程：**
```python
# 1. Clean signal
cleaned = nk.ppg_clean(ppg_raw, sampling_rate=100, method='elgendi')

# 2. Detect peaks with artifact correction
peaks, info = nk.ppg_peaks(cleaned, sampling_rate=100, correct_artifacts=True)

# 3. Assess quality
quality = nk.ppg_quality(cleaned, sampling_rate=100)

# 4. Comprehensive processing (alternative)
signals, info = nk.ppg_process(ppg_raw, sampling_rate=100)

# 5. Analyze
analysis = nk.ppg_analyze(signals, sampling_rate=100)
```

** PPG 的 HRV：**
```python
# Process PPG signal
signals, info = nk.ppg_process(ppg_raw, sampling_rate=100)

# Extract peaks and compute HRV
hrv_indices = nk.hrv(info['PPG_Peaks'], sampling_rate=100)

# PPG-derived HRV is valid but may differ slightly from ECG-derived HRV
# Differences due to pulse arrival time, vascular properties
```

## 临床和研究应用

**可穿戴健康监测：**
- 消费类智能手表和健身追踪器
- 连续心率监测
- 睡眠追踪和活动评估

**临床监测：**
- 脉搏血氧饱和度（SpO2 + 心率）
- 围手术期监测
- 重症监护心率评估

**心血管评估：**
- 脉搏波分析
- 动脉僵硬度代理（脉搏到达时间）
- 血管老化指数

**自主功能：**
- PPG 的 HRV (PPG-HRV)
- 压力和恢复监测
- 脑力负荷评估

**远程患者监护：**
- 远程医疗应用
- 家庭健康追踪
- 慢性病管理

**情感计算：**
- 从生理信号中识别情绪
- 用户体验研究
- 人机交互

## PPG 与 ECG

**PPG的优势：**
- 非侵入式，无电极
- 方便长期监测
- 低成本、小型化
- 适用于可穿戴设备

**PPG 的缺点：**
- 更容易受到运动伪影的影响
- 灌注不良时信号质量较低
- 心脏脉冲到达时间延迟
- 无法评估心电活动

**HRV 比较：**
- PPG-HRV 通常对时/频域有效
- 由于脉冲传输时间的变化，可能略有不同
- 临床 HRV（如有）首选心电图
- PPG 可接受研究和消费应用

## 解释指南

**来自 PPG 的心率：**
- 与心电图得出的心率相同的解释
- 对于速率计算来说，轻微的延迟（脉冲到达时间）可以忽略不计
- 运动伪影更常见：通过信号质量进行验证

**脉冲幅度：**
- 反映外周灌注
- 增加：血管舒张、温暖
- 减少：血管收缩、寒冷、压力、接触不良

**脉冲形态：**
- 收缩期峰值：心脏射血
- 重搏切迹：主动脉瓣关闭，动脉顺应性
- 老化/僵硬：重搏切迹更早、更明显

## 参考文献

- Elgendi, M. (2012)。指尖光电体积描记图信号的分析。当前心脏病学评论，8(1), 14-25。
- Elgendi, M.、Norton, I.、Brearley, M.、Abbott, D. 和 Schuurmans, D. (2013)。热带条件下应急响应人员测量的加速光电体积描记图中的收缩峰值检测。 PloS 一，8(10)，e76585。
- 艾伦，J. (2007)。光电体积描记法及其在临床生理测量中的应用。生理测量，28(3)，R1。
- Tamura, T.、Maeda, Y.、Sekine, M. 和 Yoshida, M. (2014)。可穿戴光电体积描记传感器——过去和现在。电子学，3(2), 282-302。