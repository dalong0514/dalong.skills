<!-- 此文件由机器翻译自 complexity.md -->

# 复杂度和熵分析

## 概述

复杂性度量量化时间序列信号的不规则性、不可预测性和多尺度结构。 NeuroKit2 提供全面的熵、分形维数和非线性动力学测量，用于评估生理信号的复杂性。

## 主要功能

### 复杂度()

同时计算多个复杂性指标以进行探索性分析。

```python
complexity_indices = nk.complexity(signal, sampling_rate=1000, show=False)
```

**退货：**
- DataFrame 具有跨类别的众多复杂性度量：
  - 熵指数
  - 分形维数
  - 非线性动力学测量
  - 信息论指标

**用例：**
- 探索性分析以确定相关措施
- 全面的信号表征
- 信号之间的比较研究

## 参数优化

在计算复杂性度量之前，应确定最佳嵌入参数：

### 复杂度_延迟()

确定相空间重建的最佳时间延迟 (τ)。

<<<代码块_1>>>

**方法：**
- `'fraser1986'`：互信息第一个最小值
- `'theiler1990'`：自相关第一个过零
- `'casdagli1991'`：曹的方法

**用途：** 在熵、吸引子重建中嵌入延迟

### 复杂度_维度()

确定最佳嵌入尺寸 (m)。

<<<代码块_2>>>

**方法：**
- `'afn'`：平均错误最近邻
- `'fnn'`：错误的最近邻居
- `'correlation'`：相关维度饱和度

**用途：** 熵计算、相空间重建

### 复杂性_容差()

确定熵测量的最佳容差 (r)。

<<<代码块_3>>>

**方法：**
- `'sd'`：基于标准差（0.1-0.25 × SD 典型值）
- `'maxApEn'`：最大化 ApEn
- `'recurrence'`：基于复发率

**用途：** 近似熵、样本熵

###复杂度_k()

确定 Higuchi 分形维数的最佳 k 参数。

<<<代码块_4>>>

**用于：** Higuchi 分形维数计算

## 熵测量

熵量化了随机性、不可预测性和信息内容。

### entropy_shannon()

香农熵 - 经典信息论测度。

<<<代码块_5>>>

**释义：**
- 更高：更随机，更不可预测
- 较低：更规律、可预测
- 单位：位（信息）

**使用案例：**
- 一般随机性评估
- 信息内容
- 信号不规则

### 熵_近似()

近似熵 (ApEn) - 模式的规律性。

<<<代码块_6>>>

**参数：**
- `delay`：时间延迟 (τ)
- `dimension`：嵌入尺寸（米）
- `tolerance`：相似度阈值 (r)

**释义：**
- 较低的 ApEn：更规则、自相似的模式
- 更高的 ApEn：更复杂、不规则
- 对信号长度敏感（建议≥100-300点）

**生理应用：**
- HRV：心脏病时 ApEn 降低
- 脑电图：神经系统疾病中的 ApEn 改变

### 熵_样本()

样本熵 (SampEn) - 改进的 ApEn。

```python
sampen = nk.entropy_sample(signal, delay=1, dimension=2, tolerance='sd')
```

**相对于 ApEn 的优势：**
- 较少依赖信号长度
- 录音之间更加一致
- 无自我匹配偏差

**释义：**
- 与 ApEn 相同但更可靠
- 大多数应用中的首选

**典型值：**
- HRV：0.5-2.5（取决于具体情况）
- 脑电图：0.3-1.5

### entropy_multiscale()

多尺度熵 (MSE) - 跨时间尺度的复杂性。

```python
mse = nk.entropy_multiscale(signal, scale=20, dimension=2, tolerance='sd',
                            method='MSEn', show=False)
```

**方法：**
- `'MSEn'`：多尺度样本熵
- `'MSApEn'`：多尺度近似熵
- `'CMSE'`：复合多尺度熵
- `'RCMSE'`：精炼复合多尺度熵

**释义：**
- 不同粗粒度尺度的熵
- 健康/复杂系统：跨多个尺度的高熵
- 有病/更简单的系统：熵减少，尤其是在更大的尺度上

**使用案例：**
- 区分真正的复杂性和随机性
- 白噪声：跨尺度恒定
- 粉红噪声/复杂性：跨尺度的结构化变化

### 熵_模糊()

模糊熵 - 使用模糊隶属函数。

```python
fuzzen = nk.entropy_fuzzy(signal, delay=1, dimension=2, tolerance='sd', r=0.2)
```

**优点：**
- 对于噪声信号更稳定
- 模式匹配的模糊边界
- 短信号下性能更佳

### 熵_排列()

排列熵 - 基于序数模式。

```python
perment = nk.entropy_permutation(signal, delay=1, dimension=3)
```

**方法：**
- 将信号编码为序数模式（排列）
- 计算模式频率
- 对噪声和非平稳性具有鲁棒性

**释义：**
- 较低：更规则的顺序结构
- 更高：更随机的排序

**使用案例：**
- 脑电图分析
- 麻醉深度监测
- 快速计算

### 熵谱()

谱熵 - 基于功率谱。

```python
spec_ent = nk.entropy_spectral(signal, sampling_rate=1000, bands=None)
```

**方法：**
- 功率谱的归一化香农熵
- 量化频率分布规律

**释义：**
- 0：单频（纯音）
- 1：白噪声（平坦频谱）

**使用案例：**
- EEG：频谱分布随状态变化
- 麻醉监测

### 熵_svd()

奇异值分解熵。

```python
svd_ent = nk.entropy_svd(signal, delay=1, dimension=2)
```

**方法：**
- 轨迹矩阵上的 SVD
- 奇异值分布的熵

**使用案例：**
- 吸引子复杂性
- 确定性与随机动力学

### 熵_微分()

微分熵 - 香农熵的连续模拟。

```python
diff_ent = nk.entropy_differential(signal)
```

**用于：** 连续概率分布

### 其他熵测量

**Tsallis 熵：**
```python
tsallis = nk.entropy_tsallis(signal, q=2)
```
- 参数为 q 的广义熵
- q=1 减少为香农熵

**Rényi 熵：**
```python
renyi = nk.entropy_renyi(signal, alpha=2)
```
- 参数为 α 的广义熵

**额外的专门熵：**
- `entropy_attention()`：注意力熵
- `entropy_grid()`：基于网格的熵
- `entropy_increment()`：增加熵
- `entropy_slope()`：斜率熵
- `entropy_dispersion()`：色散熵
- `entropy_symbolicdynamic()`：符号动力学熵
- `entropy_range()`：范围熵
- `entropy_phase()`：相熵
- `entropy_quadratic()`、`entropy_cumulative_residual()`、`entropy_rate()`：特殊变体

## 分形维数测量

分形维数表征了自相似性和粗糙度。

### fractal_katz()

Katz 分形维数 - 波形复杂性。

```python
kfd = nk.fractal_katz(signal)
```

**释义：**
- 1：直线
- >1：增加粗糙度和复杂性
- 典型范围：1.0-2.0

**优点：**
- 简单、快速的计算
- 无需参数调整

### fractal_higuchi()

Higuchi 分形维数 - 自相似性。

```python
hfd = nk.fractal_higuchi(signal, k_max=10)
```

**方法：**
- 从原始时间序列构建 k 个新时间序列
- 根据长度比例关系估计尺寸

**释义：**
- 更高的 HFD：更复杂、不规则
- 较低的 HFD：更平滑、更规则

**使用案例：**
- 脑电图复杂性
- 心率变异性分析
- 癫痫检测

### fractal_petrosian()

Petrosian 分形维数 - 快速估计。

```python
pfd = nk.fractal_petrosian(signal)
```

**优点：**
- 快速计算
- 直接计算（无需曲线拟合）

### fractal_sevcik()

Sevcik 分形维数 - 归一化波形复杂性。

```python
sfd = nk.fractal_sevcik(signal)
```

### fractal_nld()

归一化长度密度 - 基于曲线长度的测量。

```python
nld = nk.fractal_nld(signal)
```

### fractal_psdslope()

功率谱密度斜率 - 频域分形测量。

```python
slope = nk.fractal_psdslope(signal, sampling_rate=1000)
```

**方法：**
- 对数-对数功率谱的线性拟合
- 斜率 β 与分形维数有关

**释义：**
- β ≈ 0：白噪声（随机）
- β ≈ -1：粉红噪声（1/f，复数）
- β ≈ -2：布朗噪声（布朗运动）

### fractal_hurst()

赫斯特指数 - 长程依赖性。

```python
hurst = nk.fractal_hurst(signal, show=False)
```

**释义：**
- H < 0.5：反持久（均值恢复）
- H = 0.5：随机游走（白噪声）
- H > 0.5：持久（趋势、长记忆）

**使用案例：**
- 评估长期相关性
- 金融时间序列
- 心率变异性分析

### fractal_correlation()

相关维度 - 吸引子维度。

```python
corr_dim = nk.fractal_correlation(signal, delay=1, dimension=10, radius=64)
```

**方法：**
- Grassberger-Procaccia 算法
- 估计相空间中吸引子的维数

**释义：**
- 低维：确定性、低维混沌
- 高维：高维混沌或噪声

### fractal_dfa()

去趋势波动分析 - 缩放指数。

```python
dfa_alpha = nk.fractal_dfa(signal, multifractal=False, q=2, show=False)
```

**释义：**
- α < 0.5：反相关
- α = 0.5：不相关（白噪声）
- α = 1.0：1/f 噪声（粉红噪声，健康的复杂性）
- α = 1.5：布朗噪声
- α > 1.0：持续的长程相关性

**HRV 应用：**
- α1（短期，4-11次）：反映自主调节
- α2（长期，>11 次心跳）：长期相关性
- α1 降低：心脏病理学

### fractal_mfdfa()

多重分形 DFA - 多尺度分形特性。

```python
mfdfa_results = nk.fractal_mfdfa(signal, q=None, show=False)
```

**方法：**
- 将 DFA 扩展到多个 q 阶
- 表征多重分形谱

**退货：**
- 广义赫斯特指数 h(q)
- 多重分形谱f(α)
- 宽度表示多重分形强度

**使用案例：**
- 检测多重分形结构
- 健康与疾病中的 HRV 多重分形
- 脑电图多尺度动力学

### fractal_tmf()

多重分形非线性 - 与单分形的偏差。

```python
tmf = nk.fractal_tmf(signal)
```

**释义：**
- 量化与简单缩放的偏差
- 更高：更多的多重分形结构

### 分形密度()

密度分形维数。

```python
density_fd = nk.fractal_density(signal)
```

### fractal_linelength()

线长度 - 总变化量度。

```python
linelength = nk.fractal_linelength(signal)
```

**用例：**
- 简单复杂度代理
- 脑电图癫痫发作检测

## 非线性动力学

###复杂性_lyapunov()

最大李雅普诺夫指数 - 混沌和发散。

```python
lyap = nk.complexity_lyapunov(signal, delay=None, dimension=None,
                              sampling_rate=1000, show=False)
```

**释义：**
- λ < 0：稳定的定点
- λ = 0：周期轨道
- λ > 0：混乱（附近的轨迹呈指数发散）

**使用案例：**
- 检测生理信号中的混乱
- HRV：正李亚普诺夫表明非线性动力学
- EEG：癫痫检测（癫痫发作前 λ 减小）

###复杂度_lempelziv()

Lempel-Ziv Complexity - 算法复杂性。

```python
lz = nk.complexity_lempelziv(signal, symbolize='median')
```

**方法：**
- 计算不同模式的数量
- 粗粒度的随机性测量

**释义：**
- 较低：重复的、可预测的模式
- 更高：多样化、不可预测的模式

**使用案例：**
- 脑电图：意识水平、麻醉
- HRV：自主复杂性

### 复杂度_rqa()

递归量化分析 - 相空间递归。

```python
rqa_indices = nk.complexity_rqa(signal, delay=1, dimension=3, tolerance='sd')
```

**指标：**
- **复发率 (RR)**：复发状态的百分比
- **决定论 (DET)**：行中重复点的百分比
- **层流 (LAM)**：垂直结构的百分比（层流状态）
- **陷印时间 (TT)**：平均垂直线长度
- **最长对角线/垂直**：系统可预测性
- **熵 (ENTR)**：线长度分布的香农熵

**释义：**
- 高DET：确定性动力学
- 高 LAM：系统陷入特定状态
- 低 RR：随机、非循环动态

**使用案例：**
- 检测系统动力学的转变
- 生理状态变化
- 非线性时间序列分析

### Complexity_hjorth()

Hjorth 参数 - 时域复杂度。

```python
hjorth = nk.complexity_hjorth(signal)
```

**指标：**
- **活动**：信号方差
- **迁移率**：导数与信号的标准差的比例
- **复杂性**：随导数改变流动性

**使用案例：**
- 脑电图特征提取
- 癫痫检测
- 信号表征

### 复杂度去相关()

去相关时间 - 记忆持续时间。

```python
decorr_time = nk.complexity_decorrelation(signal, show=False)
```

**释义：**
- 自相关降至阈值以下的时间滞后
- 较短：快速波动、记忆力短
- 更长：波动缓慢，记忆力长

### 复杂度_相对粗糙度()

相对粗糙度 - 平滑度测量。

```python
roughness = nk.complexity_relativeroughness(signal)
```

## 信息论

### 渔民信息()

费希尔信息 - 秩序的衡量。

```python
fisher = nk.fisher_information(signal, delay=1, dimension=2)
```

**释义：**
- 高：有序、结构化
- 低：无序、随机

**使用案例：**
- 与香农熵结合（Fisher-Shannon 平面）
- 表征系统复杂性

### Fishershannon_information()

费希尔香农信息产品。

```python
fs = nk.fishershannon_information(signal)
```

**方法：**
- Fisher信息和Shannon熵的乘积
- 表征有序-无序平衡

### 相互信息()

互信息 - 变量之间共享信息。

```python
mi = nk.mutual_information(signal1, signal2, method='knn')
```

**方法：**
- `'knn'`：k-最近邻（非参数）
- `'kernel'`：核密度估计
- `'binning'`：基于直方图

**使用案例：**
- 信号之间的耦合
- 特征选择
- 非线性相关性

## 实际考虑

### 信号长度要求

|测量 |最小长度|最佳长度 |
|--------|-------------|----------------|
|香农熵 | 50 | 50 200+ |
| ApEn、SampEn | 100-300 | 500-1000 |
|多尺度熵 | 500 | 500每个规模 1000+
| DFA | 500 | 500 1000+ |
|李亚普诺夫| 1000 | 1000 5000+ |
|相关维度| 1000 | 1000 5000+ |

### 参数选择

**一般准则：**
- 首先使用参数优化功能
- 或使用常规默认值：
  - 延迟 (τ)：HRV 为 1，EEG 为自相关第一最小值
  - 尺寸（米）：2-3（典型值）
  - 公差 (r)：0.2 × SD 通用

**灵敏度：**
- 结果可以是参数敏感的
- 报告使用的参数
- 考虑敏感性分析

### 标准化和预处理

**标准化：**
- 许多措施对信号幅度敏感
- 经常推荐 Z 分数标准化
- 消除趋势可能是必要的

**平稳性：**
- 一些措施假设平稳
- 检查统计测试（例如 ADF 测试）
- 分段非平稳信号

### 解释

**取决于上下文：**
- 没有普遍的“好”或“坏”复杂性
- 比较受试者内或组间
- 考虑生理背景

**复杂性与随机性：**
- 最大熵≠最大复杂度
- 真正的复杂性：结构化的可变性
- 白噪声：高熵但低复杂度（MSE 区分）

## 应用

**心血管：**
- HRV 复杂性：减少心脏病、衰老
- DFA α1：心肌梗死后的预后标志物

**神经科学：**
- 脑电图复杂性：意识、麻醉深度
- 熵：阿尔茨海默病、癫痫、睡眠阶段
- 排列熵：麻醉监测

**心理学：**
- 抑郁、焦虑的复杂性丧失
- 在压力下增加规律性

**老化：**
- 跨系统老化带来的“复杂性损失”
- 降低多尺度复杂性

**关键转变：**
- 状态转换之前复杂性发生变化
- 预警信号（严重减速）

## 参考文献

-平卡斯，S.M.（1991）。近似熵作为系统复杂性的度量。美国国家科学院院刊，88(6), 2297-2301。
- Richman, J. S. 和 Moorman, J. R. (2000)。使用近似熵和样本熵进行生理时间序列分析。美国生理学杂志 - 心脏和循环生理学，278(6)，H2039-H2049。
- 彭 C. K. 等人（1995）。非平稳心跳时间序列中缩放指数和交叉现象的量化。混沌，5(1)，82-87。
- Costa, M.、Goldberger, A. L. 和 Peng, C. K. (2005)。生物信号的多尺度熵分析。物理评论 E，71(2)，021906。
- Grassberger, P. 和 Procaccia, I. (1983)。测量奇异吸引子的奇异程度。物理学 D：非线性现象，9(1-2), 189-208。