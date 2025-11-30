<!-- 此文件由机器翻译自 eeg.md -->

# 脑电图分析和微观状态

## 概述

分析脑电图 (EEG) 信号的频带功率、通道质量评估、源定位和微观状态识别。 NeuroKit2 与 MNE-Python 集成，以实现全面的脑电图处理工作流程。

## 核心脑电图功能

### eeg_power()

计算指定通道的标准频段的功率。

```python
power = nk.eeg_power(eeg_data, sampling_rate=250, channels=['Fz', 'Cz', 'Pz'],
                     frequency_bands={'Delta': (0.5, 4),
                                     'Theta': (4, 8),
                                     'Alpha': (8, 13),
                                     'Beta': (13, 30),
                                     'Gamma': (30, 45)})
```

**标准频段：**
- **Delta (0.5-4 Hz)**：深度睡眠，无意识过程
- **Theta (4-8 Hz)**：困倦、冥想、记忆编码
- **Alpha (8-13 Hz)**：放松的清醒状态，闭上眼睛
- **Beta (13-30 Hz)**：积极思考、专注、焦虑
- **伽玛 (30-45 Hz)**：认知处理、结合

**退货：**
- 每个通道的功率值×频带组合的数据帧
- 列：`Channel_Band`（例如“Fz_Alpha”、“Cz_Beta”）

**使用案例：**
- 静息状态分析
- 认知状态分类
- 睡眠分期
- 冥想或神经反馈监测

### eeg_badchannels()

使用统计异常值检测来识别有问题的通道。

<<<代码块_1>>>

**检测方法：**
- 跨渠道的标准差异常值
- 与其他渠道的相关性
- 平坦或死通道
- 噪音过大的通道

**参数：**
- `bad_threshold`：异常值检测的 Z 分数阈值（默认值：2）

**退货：**
- 被识别为有问题的频道名称列表

**用例：**
- 分析前的质量控制
- 自动不良通道拒绝
- 插值或排除决策

### eeg_rereference()

重新表达相对于不同参考点的电压测量值。

<<<代码块_2>>>

**参考类型：**
- `'average'`：平均参考（所有电极的平均值）
- `'REST'`：参比电极标准化技术
- `'bipolar'`：电极对之间的差分记录
- 特定通道名称：使用单电极作为参考

**常用参考：**
- **平均参考**：高密度脑电图最常见
- **相连乳突**：传统临床脑电图
- **Vertex (Cz)**：有时用于 ERP 研究
- **REST**：近似无穷大参考

**退货：**
- 重新参考脑电图数据

### eeg_gfp()

计算全局场功率 - 每个时间点所有电极的标准偏差。

<<<代码块_3>>>

**释义：**
- 高 GFP：跨区域的强而同步的大脑活动
- 低 GFP：活动弱或不同步
- GFP 峰：稳定形貌的点，用于微状态检测

**使用案例：**
- 识别稳定地形模式的时期
- 选择微观状态分析的时间点
- 事件相关电位（ERP）可视化

### eeg_diss()

测量电场配置之间的地形差异。

<<<代码块_4>>>

**方法：**
- 基于 GFP：标准化差异
- 空间相关性
- 余弦距离

**用例：**
- 比较条件之间的地形
- 微观状态转变分析
- 模板匹配

## 源本地化

### eeg_source()

执行源重建以根据头皮记录估计大脑水平的活动。

<<<代码块_5>>>

**方法：**
- `'sLORETA'`：标准化低分辨率电磁断层扫描
  - 点源零定位误差
  - 良好的空间分辨率
- `'MNE'`：最小范数估计
  - 快速、完善
  - 偏向于肤浅的来源
- `'dSPM'`：动态统计参数映射
  - 标准化跨国公司
- `'eLORETA'`：精确的 LORETA
  - 提高定位精度

**要求：**
- 正演模型（前导场矩阵）
- 共同注册的电极位置
- 头部模型（边界元或球形）

**退货：**
- 源空间活动估计

### eeg_source_extract()

从特定的大脑解剖区域提取活动。

<<<代码块_6>>>

**区域选项：**
- 标准地图集：Desikan-Killiany、Destrieux、AAL
- 自定义投资回报率
- 布罗德曼地区

**退货：**
- 每个地区的时间序列
- 体素的平均或主成分

**使用案例：**
- 感兴趣区域分析
- 功能连接
- 源级统计

## 微观状态分析

微观状态是稳定大脑结构的短暂（80-120 ms）周期，代表协调的神经网络。通常有 4-7 个具有不同功能的微状态类别（通常标记为 A、B、C、D）。

### microstates_segment()

使用聚类算法识别和提取微观状态。
```python
microstates = nk.microstates_segment(eeg_data, n_microstates=4, sampling_rate=250,
                                      method='kmod', normalize=True)
```

**方法：**
- `'kmod'`（默认）：针对 EEG 地形优化的修改 k 均值
  - 极性不变聚类
  - 在微观状态文献中最常见
- `'kmeans'`：标准 k 均值聚类
- `'kmedoids'`：K-medoids（对异常值更稳健）
- `'pca'`：主成分分析
- `'ica'`：独立成分分析
- `'aahc'`：原子化和凝聚层次聚类

**参数：**
- `n_microstates`：微状态类别的数量（通常为 4-7）
- `normalize`：标准化拓扑（推荐：True）
- `n_inits`：随机初始化的数量（增加稳定性）

**退货：**
- 字典：
  - `'maps'`：微观状态模板拓扑
  - `'labels'`：每个时间点的微观状态标签
  - `'gfp'`：全局场功率
  - `'gev'`：全局解释方差

### microstates_findnumber()

估计最佳微观状态数。

```python
optimal_k = nk.microstates_findnumber(eeg_data, show=True)
```

**标准：**
- **全局解释方差 (GEV)**：解释的方差百分比
  - 肘部法：在GEV曲线中找到“膝盖”
  - 通常实现 70-80% GEV
- **Krzanowski-Lai (KL) 标准**：平衡拟合和简约性的统计测量
  - 最大KL表示最佳k

**典型范围：** 4-7 个微观状态
- 4 种微观状态：经典 A、B、C、D 状态
- 5-7 个微观状态：更细粒度的分解

### microstates_classify()

根据前后和左右通道值重新排序微观状态。

```python
classified = nk.microstates_classify(microstates)
```

**目的：**
- 标准化跨学科的微观状态标签
- 匹配常规 A、B、C、D 地形：
  - **A**：左右方向，顶枕骨
  - **B**：左右方向，额颞叶
  - **C**：前后方向，额中央
  - **D**：前-中、前-后（C 的逆）

**退货：**
- 重新排序微观状态图和标签

### microstates_clean()

预处理脑电图数据以进行微观状态提取。

```python
cleaned_eeg = nk.microstates_clean(eeg_data, sampling_rate=250)
```

**预处理步骤：**
- 带通滤波（通常为 2-20 Hz）
- 文物拒绝
- 不良通道插值
- 重新参考平均值

**理由：**
- 微观状态反映大规模网络活动
- 高频和低频伪影可能会扭曲地形

### microstates_peaks()

识别 GFP 峰以进行微状态分析。

```python
peak_indices = nk.microstates_peaks(eeg_data, sampling_rate=250)
```

**目的：**
- 通常在 GFP 峰值处分析微观状态
- 峰值代表最大、稳定地形活动的时刻
- 减少计算负载和噪声敏感性

**退货：**
- GFP局部最大值的指数

### microstates_static()

计算各个微观状态的时间特性。

```python
static_metrics = nk.microstates_static(microstates)
```

**指标：**
- **持续时间（毫秒）**：每个微状态花费的平均时间
  - 典型：80-120 毫秒
  - 体现稳定性和持久性
- **发生次数（每秒）**：微观状态出现的频率
  - 进入每个状态的频率
- **覆盖率 (%)**：每个微观状态总时间的百分比
  - 相对优势
- **全局解释方差 (GEV)**：每个类别解释的方差
  - 模板拟合质量

**退货：**
- 包含每个微状态类指标的 DataFrame

**释义：**
- 持续时间的变化：改变网络稳定性
- 发生的变化：改变状态动态
- 覆盖范围的变化：特定网络的主导地位

### microstates_dynamic()

分析微观状态之间的转变模式。

```python
dynamic_metrics = nk.microstates_dynamic(microstates)
```

**指标：**
- **转移矩阵**：从状态 i 转移到状态 j 的概率
  - 揭示优先顺序
- **转换率**：总体转换频率
  - 速率更高：切换更快速
- **熵**：转变的随机性
  - 高熵：不可预测的切换
  - 低熵：定型序列
- **马尔可夫测试**：转换是否依赖于历史？

**退货：**
- 带有转换统计的字典

**使用案例：**
- 识别临床人群中的异常微状态序列
- 网络动态性和灵活性
- 依赖于状态的信息处理

### microstates_plot()

可视化微观状态拓扑和时间进程。

```python
nk.microstates_plot(microstates, eeg_data)
```

**显示：**
- 每个微观状态类别的地形图
- 带有微状态标签的 GFP 迹线
- 显示状态序列的转换图
- 统计汇总

## 跨国公司整合实用程序

### mne_data()

从 MNE-Python 访问示例数据集。
```python
raw = nk.mne_data(dataset='sample', directory=None)
```

**可用数据集：**
- `'sample'`：多模态 (MEG/EEG) 示例
- `'ssvep'`：稳态视觉诱发电位
- `'eegbci'`：运动想象 BCI 数据集

### mne_to_df() / mne_to_dict()

将 MNE 对象转换为 NeuroKit 兼容格式。

```python
df = nk.mne_to_df(raw)
data_dict = nk.mne_to_dict(epochs)
```

**用例：**
- 在 NeuroKit2 中处理 MNE 处理的数据
- 在格式之间转换以进行分析

### mne_channel_add() / mne_channel_extract()

管理跨国公司对象中的各个渠道。

```python
# Extract specific channels
subset = nk.mne_channel_extract(raw, ['Fz', 'Cz', 'Pz'])

# Add derived channels
raw_with_eog = nk.mne_channel_add(raw, new_channel_data, ch_name='EOG')
```

### mne_crop()

按时间或样本修剪录音。

```python
cropped = nk.mne_crop(raw, tmin=10, tmax=100)
```

### mne_templateMRI()

为源本地化提供模板剖析。

```python
subjects_dir = nk.mne_templateMRI()
```

**用例：**
- 无需单独 MRI 进行源分析
- 集团级源码本地化
- fsaverage 模板大脑

### eeg_simulate()

生成用于测试的合成脑电图信号。

```python
synthetic_eeg = nk.eeg_simulate(duration=60, sampling_rate=250, n_channels=32)
```

## 实际考虑

### 采样率建议
- **最小值**：100 Hz 用于基本功率分析
- **标准**：大多数应用为 250-500 Hz
- **高分辨率**：1000+ Hz 可提供详细的时间动态

### 录音时长
- **功率分析**：稳定估计≥2 分钟
- **微观状态**：≥2-5 分钟，优选更长的时间
- **静息状态**：典型 3-10 分钟
- **与事件相关**：取决于试验次数（每个条件≥30 次试验）

### 工件管理
- **眨眼**：通过 ICA 或回归消除
- **肌肉伪影**：高通滤波器（≥1 Hz）或手动抑制
- **不良通道**：在分析之前检测并插值
- **线路噪声**：陷波滤波器，50/60 Hz

### 最佳实践

**功率分析：**
```python
# 1. Clean data
cleaned = nk.signal_filter(eeg_data, sampling_rate=250, lowcut=0.5, highcut=45)

# 2. Identify and interpolate bad channels
bad = nk.eeg_badchannels(cleaned, sampling_rate=250)
# Interpolate bad channels using MNE

# 3. Re-reference
rereferenced = nk.eeg_rereference(cleaned, reference='average')

# 4. Compute power
power = nk.eeg_power(rereferenced, sampling_rate=250, channels=channel_list)
```

**微观状态工作流程：**
```python
# 1. Preprocess
cleaned = nk.microstates_clean(eeg_data, sampling_rate=250)

# 2. Determine optimal number of states
optimal_k = nk.microstates_findnumber(cleaned, show=True)

# 3. Segment microstates
microstates = nk.microstates_segment(cleaned, n_microstates=optimal_k,
                                     sampling_rate=250, method='kmod')

# 4. Classify to standard labels
microstates = nk.microstates_classify(microstates)

# 5. Compute temporal metrics
static = nk.microstates_static(microstates)
dynamic = nk.microstates_dynamic(microstates)

# 6. Visualize
nk.microstates_plot(microstates, cleaned)
```

## 临床和研究应用

**认知神经科学：**
- 注意力、工作记忆、执行功能
- 语言处理
- 感官知觉

**临床人群：**
- 癫痫：癫痫发作检测、定位
- 阿尔茨海默病：脑电图减慢、微观状态改变
- 精神分裂症：微观状态改变，尤其是 C 状态
- ADHD：θ/β 比率增加
- 抑郁症：额叶α不对称

**意识研究：**
- 麻醉监测
- 意识障碍
- 睡眠分期

**神经反馈：**
- 实时频段训练
- 阿尔法增强放松
- 专注 Beta 增强

## 参考文献

- 米歇尔，C.M.，和科尼格，T.（2018）。脑电图微状态作为研究全脑神经元网络时间动态的工具：综述。神经影像，180, 577-593。
- Pascual-Marqui, R. D.、Michel, C. M. 和 Lehmann, D. (1995)。将脑电活动分割为微观状态：模型估计和验证。 IEEE 生物医学工程汇刊，42(7), 658-665。
- Gramfort, A.、Luessi, M.、Larson, E.、Engemann, D. A.、Strohmeier, D.、Brodbeck, C., ... & Hämäläinen, M. (2013)。使用 MNE-Python 进行脑磁图和脑电图数据分析。神经科学前沿，7, 267。