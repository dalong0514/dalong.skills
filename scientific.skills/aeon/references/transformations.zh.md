<!-- 此文件由机器翻译自 transformations.md -->

# 转换

Aeon 提供了广泛的转换功能，用于从时间序列数据进行预处理、特征提取和表示学习。

## 转换类型

Aeon 的区别在于：
- **CollectionTransformers**：转换多个时间序列（集合）
- **SeriesTransformers**：转换单个时间序列

## 收集变形金刚

### 基于卷积的特征提取

使用随机内核快速、可扩展的特征生成：

- `RocketTransformer` - 随机卷积核
- `MiniRocketTransformer` - 简化的 ROCKET 以提高速度
- `MultiRocketTransformer` - 增强型 ROCKET 变体
- `HydraTransformer` - 多分辨率扩张卷积
- `MultiRocketHydraTransformer` - 结合 ROCKET 和 Hydra
- `ROCKETGPU` - GPU 加速变体

**使用时机**：任何 ML 算法都需要快速、可扩展的功能、强大的基线性能。

### 统计特征提取

基于时间序列特征的领域无关特征：

- `Catch22` - 22 个规范时间序列特征
- `TSFresh` - 全面的自动特征提取（100 多个特征）
- `TSFreshRelevant` - 通过相关性过滤进行特征提取
- `SevenNumberSummary` - 描述性统计（平均值、标准差、分位数）

**何时使用**：需要可解释的功能、与领域无关的方法或提供传统的机器学习。

### 基于字典的表示

离散表示的符号近似：

- `SAX` - 符号聚合近似
- `PAA` - 分段聚合近似
- `SFA` - 符号傅里叶近似
- `SFAFast` - 优化的 SFA
- `SFAWhole` - 整个系列的 SFA（无窗口）
- `BORF` - 感受野袋

**使用时**：需要离散/符号表示、降维、可解释性。

### 基于 Shapelet 的特征

判别子序列提取：

- `RandomShapeletTransform` - 随机判别 shapelet
- `RandomDilatedShapeletTransform` - 用于多尺度的扩张 shapelet
- `SAST` - 可扩展且准确的子序列变换
- `RSAST` - 随机 SAST

**使用时**：需要可解释的判别模式、相位不变特征。

### 基于间隔的特征

时间间隔的统计摘要：

- `RandomIntervals` - 随机间隔的特征
- `SupervisedIntervals` - 监督间隔选择
- `QUANTTransformer` - 基于分位数的区间特征

**使用时间**：本地化到特定窗口的预测模式。

### 预处理转换

数据准备和标准化：

- `MinMaxScaler` - 缩放至 [0, 1] 范围
- `Normalizer` - Z 归一化（零均值，单位方差）
- `Centerer` - 中心到零均值
- `SimpleImputer` - 填充缺失值
- `DownsampleTransformer` - 降低时间分辨率
- `Tabularizer` - 将时间序列转换为表格格式

**使用时**：需要标准化、缺失值处理、格式转换。

### 专业化转型

先进的分析方法：

- `MatrixProfile` - 计算距离剖面以发现模式
- `DWTTransformer` - 离散小波变换
- `AutocorrelationFunctionTransformer` - ACF 计算
- `Dobin` - 使用邻居的基于距离的异常值 BasIs
- `SignatureTransformer` - 路径签名方法
- `PLATransformer` - 分段线性逼近

### 类不平衡处理

- `ADASYN` - 自适应合成采样
- `SMOTE` - 合成少数过采样
- `OHIT` - 时间序列高度不平衡的过采样

**使用时间**：类别不平衡的分类。

### 管道组成

- `CollectionTransformerPipeline` - 链接多个变压器

## 系列变压器

转换单个时间序列（例如，用于预测中的预处理）。

### 统计分析

- `AutoCorrelationSeriesTransformer` - 自相关
- `StatsModelsACF` - 使用 statsmodels 的 ACF
- `StatsModelsPACF` - 部分自相关

### 平滑和过滤

- `ExponentialSmoothing` - 指数加权移动平均线
- `MovingAverage` - 简单或加权移动平均线
- `SavitzkyGolayFilter` - 多项式平滑
- `GaussianFilter` - 高斯核平滑
- `BKFilter` - Baxter-King 带通滤波器
- `DiscreteFourierApproximation` - 基于傅里叶的过滤

**何时使用**：需要降噪、趋势提取或频率过滤。
### 降维

- `PCASeriesTransformer` - 主成分分析
- `PlASeriesTransformer` - 分段线性逼近

### 转换

- `BoxCoxTransformer` - 方差稳定
- `LogTransformer` - 对数缩放
- `ClaSPTransformer` - 分类分数配置文件

### 管道组成

- `SeriesTransformerPipeline` - 链式系列变压器

## 快速入门：特征提取

```python
from aeon.transformations.collection.convolution_based import RocketTransformer
from aeon.classification.sklearn import RotationForest
from aeon.datasets import load_classification

# Load data
X_train, y_train = load_classification("GunPoint", split="train")
X_test, y_test = load_classification("GunPoint", split="test")

# Extract ROCKET features
rocket = RocketTransformer()
X_train_features = rocket.fit_transform(X_train)
X_test_features = rocket.transform(X_test)

# Use with any sklearn classifier
clf = RotationForest()
clf.fit(X_train_features, y_train)
accuracy = clf.score(X_test_features, y_test)
```

## 快速入门：预处理管道

<<<代码块_1>>>

## 快速入门：系列平滑

<<<代码块_2>>>

## 算法选择

### 对于特征提取：
- **速度 + 性能**：MiniRocketTransformer
- **可解释性**：Catch22、TSFresh
- **降维**：PAA、SAX、PCA
- **判别模式**：Shapelet 变换
- **全面的功能**：TSFresh（具有更长的运行时间）

### 对于预处理：
- **标准化**：标准化器、MinMaxScaler
- **平滑**：移动平均、SavitzkyGolayFilter
- **缺失值**：SimpleImputer
- **频率分析**：DWTTransformer、傅里叶方法

### 对于符号表示：
- **快速逼近**：PAA
- **基于字母表**：SAX
- **基于频率**：SFA、SFAFast

## 最佳实践

1. **仅适合训练数据**：避免数据泄漏
   <<<代码块_3>>>

2. **管道组成**：复杂工作流程的链式变压器
   <<<代码块_4>>>

3. **特征选择**：TSFresh可以生成很多特征；考虑选择
   <<<代码块_5>>>

4. **内存考虑因素**：一些 Transformer 在大型数据集上占用大量内存
   - 使用 MiniRocket 代替 ROCKET 来提高速度
   - 考虑对很长的系列进行下采样
   - 使用ROCKETGPU进行GPU加速

5. **领域知识**：选择与领域匹配的转换：
   - 周期性数据：基于傅里叶的方法
   - 噪声数据：平滑滤波器
   - 尖峰检测：小波变换