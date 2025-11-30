<!-- 此文件由机器翻译自 signal_processing.md -->

# 通用信号处理

## 概述

NeuroKit2 提供适用于任何时间序列数据的全面信号处理实用程序。这些函数支持适用于所有信号类型的滤波、转换、峰值检测、分解和分析操作。

## 预处理函数

### signal_filter()

应用频域滤波来消除噪声或隔离频段。

```python
filtered = nk.signal_filter(signal, sampling_rate=1000, lowcut=None, highcut=None,
                            method='butterworth', order=5)
```

**滤波器类型（通过低切/高切组合）：**

**低通**（仅高切）：
<<<代码块_1>>>
- 删除高切以上的频率
- 平滑信号，消除高频噪声

**高通**（仅低切）：
<<<代码块_2>>>
- 删除低切以下的频率
- 消除基线漂移、直流偏移

**带通**（低切和高切）：
<<<代码块_3>>>
- 保留低切和高切之间的频率
- 隔离特定频段

**带阻/缺口**（电源线移除）：
<<<代码块_4>>>
- 消除 50 或 60 Hz 电力线噪声
- 窄陷波滤波器

**方法：**
- `'butterworth'`（默认）：平滑频率响应，平坦通带
- `'bessel'`：线性相位，最小振铃
- `'chebyshev1'`：更陡峭的滚降，通带中的纹波
- `'chebyshev2'`：陡峭的滚降，阻带中的纹波
- `'elliptic'`：两个频段中最陡的滚降和波纹
- `'powerline'`：50/60 Hz 陷波滤波器

**订购参数：**
- 更高阶：过渡更陡峭，振铃更大
- 低阶：过渡更平缓，振铃更少
- 典型：生理信号 2-5

### signal_sanitize()

删除无效值（NaN、inf）并可选择进行插值。

<<<代码块_5>>>

**使用案例：**
- 处理缺失的数据点
- 删除标记为 NaN 的工件
- 为需要连续数据的算法准备信号

### signal_resample()

更改信号的采样率（上采样或下采样）。

<<<代码块_6>>>

**方法：**
- `'interpolation'`：三次样条插值
- `'FFT'`：频域重采样
- `'poly'`：多相滤波（最适合下采样）

**使用案例：**
- 匹配多模式录音的采样率
- 减少数据大小（下采样）
- 提高时间分辨率（上采样）

### signal_fillmissing()

插入缺失或无效的数据点。

```python
filled = nk.signal_fillmissing(signal, method='linear')
```

**方法：**
- `'linear'`：线性插值
- `'nearest'`：最近邻居
- `'pad'`：向前/向后填充
- `'cubic'`：三次样条
- `'polynomial'`：多项式拟合

## 转换函数

### signal_detrend()

从信号中消除缓慢趋势。

```python
detrended = nk.signal_detrend(signal, method='polynomial', order=1)
```

**方法：**
- `'polynomial'`：拟合并减去多项式（阶数 1 = 线性）
- `'loess'`：局部加权回归
- `'tarvainen2002'`：平滑先于去趋势

**使用案例：**
- 消除基线漂移
- 分析前稳定均值
- 为平稳性假设算法做好准备

### signal_decompose()

将信号分解为组成部分。

```python
components = nk.signal_decompose(signal, sampling_rate=1000, method='emd')
```

**方法：**

**经验模态分解 (EMD)：**
```python
components = nk.signal_decompose(signal, sampling_rate=1000, method='emd')
```
- 数据自适应分解为本征模态函数 (IMF)
- 每个 IMF 代表不同的频率内容（从高到低）
- 没有预定义的基函数

**奇异谱分析 (SSA)：**
```python
components = nk.signal_decompose(signal, method='ssa')
```
- 分解为趋势、振荡和噪音
- 基于轨迹矩阵特征值分解

**小波分解：**
- 时频表示
- 时间和频率均本地化

**退货：**
- 带有分量信号的字典
- 趋势、振荡成分、残差

**使用案例：**
- 隔离生理节律
- 将信号与噪声分开
- 多尺度分析

### signal_recompose()

从分解的成分中重建信号。

```python
reconstructed = nk.signal_recompose(components, indices=[1, 2, 3])
```

**用例：**
- 分解后选择性重建
- 删除特定的 IMF 或组件
- 自适应过滤

### signal_binarize()

根据阈值将连续信号转换为二进制 (0/1)。

```python
binary = nk.signal_binarize(signal, method='threshold', threshold=0.5)
```

**方法：**
- `'threshold'`：简单阈值
- `'median'`：基于中位数
- `'mean'`：基于均值
- `'quantile'`：基于百分位

**用例：**
- 从连续信号中检测事件
- 触发提取
- 国家分类

### 信号失真()

添加受控噪声或伪影以进行测试。

```python
distorted = nk.signal_distort(signal, sampling_rate=1000, noise_amplitude=0.1,
                              noise_frequency=50, artifacts_amplitude=0.5)
```

**参数：**
- `noise_amplitude`：高斯噪声级别
- `noise_frequency`：正弦干扰（例如电力线）
- `artifacts_amplitude`：随机尖峰伪影
- `artifacts_number`：要添加的工件数量

**使用案例：**
- 算法稳健性测试
- 预处理方法评估
- 真实的数据模拟

### signal_interpolate()

在新时间点插入信号或填补空白。

```python
interpolated = nk.signal_interpolate(x_values, y_values, x_new=None, method='quadratic')
```

**方法：**
- `'linear'`、`'quadratic'`、`'cubic'`：多项式插值
- `'nearest'`：最近邻居
- `'monotone_cubic'`：保持单调性

**用例：**
- 将不规则样本转换为规则网格
- 上采样以实现可视化
- 将信号与不同的时基对齐

### signal_merge()

组合具有不同采样率的多个信号。

```python
merged = nk.signal_merge(signal1, signal2, time1=None, time2=None, sampling_rate=None)
```

**用例：**
- 多模态信号集成
- 合并来自不同设备的数据
- 基于时间戳同步

### signal_flatline()

识别持续信号的周期（伪影或传感器故障）。

```python
flatline_mask = nk.signal_flatline(signal, duration=5.0, sampling_rate=1000)
```

**退货：**
- 二进制掩码，其中 True 表示平坦周期
- 持续时间阈值可防止正常稳定性出现误报

### 信号噪声()

向信号中添加各种类型的噪声。

```python
noisy = nk.signal_noise(signal, sampling_rate=1000, noise_type='gaussian',
                        amplitude=0.1)
```

**噪音类型：**
- `'gaussian'`：白噪声
- `'pink'`：1/f 噪声（常见于生理信号）
- `'brown'`：布朗（随机游走）
- `'powerline'`：正弦干扰 (50/60 Hz)

### signal_surrogate()

生成保留某些属性的代理信号。

```python
surrogate = nk.signal_surrogate(signal, method='IAAFT')
```

**方法：**
- `'IAAFT'`：迭代幅度调整傅里叶变换
  - 保留幅度分布和功率谱
- `'random_shuffle'`：随机排列（零假设检验）

**用例：**
- 非线性测试
- 统计检验的零假设生成

## 峰值检测和校正

### signal_findpeaks()

检测信号中的局部最大值（峰值）。

```python
peaks_dict = nk.signal_findpeaks(signal, height_min=None, height_max=None,
                                 relative_height_min=None, relative_height_max=None)
```

**关键参数：**
- `height_min/max`：绝对幅度阈值
- `relative_height_min/max`：相对于信号范围 (0-1)
- `threshold`：最低显着性
- `distance`：峰值之间的最小样本

**退货：**
- 字典：
  - `'Peaks'`：峰值指数
  - `'Height'`：峰值幅度
  - `'Distance'`：峰间间隔

**使用案例：**
- 任何信号的通用峰值检测
- R 峰值、呼吸峰值、脉搏峰值
- 事件检测

### signal_fixpeaks()

纠正检测到的伪影和异常峰。

```python
corrected = nk.signal_fixpeaks(peaks, sampling_rate=1000, iterative=True,
                               method='Kubios', interval_min=None, interval_max=None)
```

**方法：**
- `'Kubios'`：Kubios HRV 软件方法（默认）
- `'Malik1996'`：工作组标准 (1996)
- `'Kamath1993'`：卡马斯的方法

**更正：**
- 删除生理上不合理的间隔
- 插入缺失的峰值
- 删除额外检测到的峰值（重复）

**用例：**
- R-R间隔中的伪影校正
- 提高 HRV 分析质量
- 呼吸或脉搏峰值校正

## 分析函数

### signal_rate()

根据事件发生（峰值）计算瞬时速率。

```python
rate = nk.signal_rate(peaks, sampling_rate=1000, desired_length=None)
```

**方法：**
- 计算事件间隔
- 转换为每分钟事件数
- 插值以匹配所需的长度

**用例：**
- R 峰值心率
- 呼吸峰值的呼吸频率
- 任何周期性事件发生率

### signal_period()

查找信号中的主导周期/频率。

```python
period = nk.signal_period(signal, sampling_rate=1000, method='autocorrelation',
                          show=False)
```

**方法：**
- `'autocorrelation'`：自相关函数中的峰值
- `'powerspectraldensity'`：频谱峰值

**退货：**
- 以样本或秒为单位的周期
- 频率（1/周期），单位 Hz

**用例：**
- 检测主导节奏
- 估计基频
- 呼吸频率、心率估计

### 信号相位()

计算信号的瞬时相位。

```python
phase = nk.signal_phase(signal, method='hilbert')
```

**方法：**
- `'hilbert'`：希尔伯特变换（解析信号）
- `'wavelet'`：基于小波的相位

**退货：**
- 以弧度表示的相位（-π 到 π）或 0 到 1（标准化）

**使用案例：**
- 锁相分析
- 同步措施
- 相位幅度耦合

### signal_psd()

计算功率谱密度。

```python
psd, freqs = nk.signal_psd(signal, sampling_rate=1000, method='welch',
                           max_frequency=None, show=False)
```

**方法：**
- `'welch'`：韦尔奇周期图（加窗 FFT，默认）
- `'multitapers'`：多锥度方法（卓越的谱估计）
- `'lomb'`：Lomb-Scargle（不均匀采样数据）
- `'burg'`：自回归（参数）

**退货：**
- `psd`：每个频率的功率（单位²/Hz）
- `freqs`：频率箱 (Hz)

**用例：**
- 频率内容分析
- HRV频域
- 光谱特征

### 信号功率()

计算特定频段的功率。

```python
power_dict = nk.signal_power(signal, sampling_rate=1000, frequency_bands={
    'VLF': (0.003, 0.04),
    'LF': (0.04, 0.15),
    'HF': (0.15, 0.4)
}, method='welch')
```

**退货：**
- 包含每频段绝对和相对功率的字典
- 峰值频率

**用例：**
- HRV频率分析
- 脑电图带功率
- 节奏量化

### signal_autocor()

计算自相关函数。

```python
autocorr = nk.signal_autocor(signal, lag=1000, show=False)
```

**释义：**
- 滞后时的高自相关性：信号重复每个滞后样本
- 周期性信号：在周期的倍数处出现峰值
- 随机信号：快速衰减到零

**使用案例：**
- 检测周期性
- 评估时间结构
- 信号记忆

### signal_zerocrossings()

计算过零（符号变化）。

```python
n_crossings = nk.signal_zerocrossings(signal)
```

**释义：**
- 更多交叉：更高频率的内容
- 与主频率相关（粗略估计）

**用例：**
- 简单的频率估计
- 信号规律性评估

### signal_changepoints()

检测信号属性的突然变化（均值、方差）。

```python
changepoints = nk.signal_changepoints(signal, penalty=10, method='pelt', show=False)
```

**方法：**
- `'pelt'`：修剪精确线性时间（快速、精确）
- `'binseg'`：二进制分割（更快，近似）

**参数：**
- `penalty`：控制灵敏度（越高=变化点越少）

**退货：**
- 检测到的变化点的索引
- 变化点之间的段

**使用案例：**
- 将信号分段为状态
- 检测转变（例如睡眠阶段、觉醒状态）
- 自动纪元定义

### signal_synchrony()

评估两个信号之间的同步。

```python
sync = nk.signal_synchrony(signal1, signal2, method='correlation')
```

**方法：**
- `'correlation'`：皮尔逊相关性
- `'coherence'`：频域相干性
- `'mutual_information'`：信息理论测量
- `'phase'`：锁相值

**使用案例：**
- 心脑耦合
- 脑间同步
- 多渠道协调

### signal_smooth()

应用平滑来减少噪音。

```python
smoothed = nk.signal_smooth(signal, method='convolution', kernel='boxzen', size=10)
```

**方法：**
- `'convolution'`：应用内核（boxcar、Gaussian 等）
- `'median'`：中值滤波器（对异常值具有鲁棒性）
- `'savgol'`：Savitzky-Golay 滤波器（保留峰值）
- `'loess'`：局部加权回归

**内核类型（用于卷积）：**
- `'boxcar'`：简单移动平均线
- `'gaussian'`：高斯加权平均值
- `'hann'`、`'hamming'`、`'blackman'`：窗口函数

**使用案例：**
- 降噪
- 趋势提取
- 可视化增强

### signal_timeFrequency()

时频表示（频谱图）。

```python
tf, time, freq = nk.signal_timefrequency(signal, sampling_rate=1000, method='stft',
                                        max_frequency=50, show=False)
```

**方法：**
- `'stft'`：短时傅立叶变换
- `'cwt'`：连续小波变换

**退货：**
- `tf`：时频矩阵（各时频点功率）
- `time`：时间段
- `freq`：频率箱

**使用案例：**
- 非平稳信号分析
- 随时间变化的频率内容
- EEG/MEG时频分析

## 模拟

### signal_simulate()

生成各种合成信号以进行测试。

```python
signal = nk.signal_simulate(duration=10, sampling_rate=1000, frequency=[5, 10],
                            amplitude=[1.0, 0.5], noise=0.1)
```

**信号类型：**
- 正弦振荡（指定频率）
- 多频率成分
- 高斯噪声
- 组合

**使用案例：**
- 算法测试
- 方法验证
- 教育示范

## 可视化

### signal_plot()

可视化信号和可选标记。

```python
nk.signal_plot(signal, sampling_rate=1000, peaks=None, show=True)
```

**特点：**
- 以秒为单位的时间轴
- 峰值标记
- 信号阵列的多个子图

## 实用技巧

**选择过滤器参数：**
- **低切**：设置低于感兴趣的最低频率
- **Highcut**：设置高于感兴趣的最高频率
- **顺序**：从 2-5 开始，如果过渡太慢则增加
- **方法**：巴特沃斯是安全默认值

**处理边缘效应：**
- 滤波在信号边缘引入伪影
- 过滤前填充信号，然后修剪
- 或丢弃最初/最后几秒

**处理差距：**
- 小间隙：`signal_fillmissing()` 带插值
- 大间隙：分段信号，单独分析
- 将间隙标记为 NaN，谨慎使用插值

**组合操作：**
```python
# Typical preprocessing pipeline
signal = nk.signal_sanitize(raw_signal)  # Remove invalid values
signal = nk.signal_filter(signal, sampling_rate=1000, lowcut=0.5, highcut=40)  # Bandpass
signal = nk.signal_detrend(signal, method='polynomial', order=1)  # Remove linear trend
```

**性能考虑：**
- 滤波：基于 FFT 的方法对于长信号更快
- 重采样：在管道早期进行下采样以加快速度
- 大型数据集：如果内存有限，则分块处理

## 参考文献
- Virtanen, P. 等人。 （2020）。 SciPy 1.0：Python 中科学计算的基本算法。自然方法，17(3), 261-272。
- Tarvainen, M. P.、Ranta-aho, P. O. 和 Karjalainen, P. A. (2002)。一种应用于 HRV 分析的先进去趋势方法。 IEEE 生物医学工程汇刊，49(2), 172-175。
- Huang, N.E. 等人（1998）。用于非线性和非平稳时间序列分析的经验模态分解和希尔伯特谱。伦敦皇家学会会议记录 A，454（1971），903-995。