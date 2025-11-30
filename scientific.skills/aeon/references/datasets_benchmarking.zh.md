<!-- 此文件由机器翻译自 datasets_benchmarking.md -->

# 数据集和基准测试

Aeon 提供了用于加载数据集和基准时间序列算法的综合工具。

## 数据集加载

### 特定任务加载器

**分类数据集**：
```python
from aeon.datasets import load_classification

# Load train/test split
X_train, y_train = load_classification("GunPoint", split="train")
X_test, y_test = load_classification("GunPoint", split="test")

# Load entire dataset
X, y = load_classification("GunPoint")
```

**回归数据集**：
<<<代码块_1>>>

**预测数据集**：
<<<代码块_2>>>

**异常检测数据集**：
<<<代码块_3>>>

### 文件格式加载器

**从 .ts 文件加载**：
<<<代码块_4>>>

**从 .tsf 文件加载**：
<<<代码块_5>>>

**从 ARFF 文件加载**：
<<<代码块_6>>>

**从 TSV 文件加载**：
```python
from aeon.datasets import load_from_tsv_file

data = load_from_tsv_file("path/to/data.tsv")
```

**加载时间评估 CSV**：
```python
from aeon.datasets import load_from_timeeval_csv_file

X, y = load_from_timeeval_csv_file("path/to/timeeval.csv")
```

### 写入数据集

**写入.ts格式**：
```python
from aeon.datasets import write_to_ts_file

write_to_ts_file(X, "output.ts", y=y, problem_name="MyDataset")
```

**写入 ARFF 格式**：
```python
from aeon.datasets import write_to_arff_file

write_to_arff_file(X, "output.arff", y=y)
```

## 内置数据集

Aeon 包含多个用于快速测试的基准数据集：

### 分类
- `ArrowHead` - 形状分类
- `GunPoint` - 手势识别
- `ItalyPowerDemand` - 能源需求
- `BasicMotions` - 运动分类
- 以及 UCR/UEA 档案中的 100 多个内容

### 回归
- `Covid3Month` - 新冠肺炎预测
- 来自 Monash TSER 档案的各种数据集

### 细分
- 时间序列分割数据集
- 人类活动数据
- 传感器数据收集

### 特别收藏
- `RehabPile` - 康复数据（分类和回归）

## 数据集元数据

获取有关数据集的信息：

```python
from aeon.datasets import get_dataset_meta_data

metadata = get_dataset_meta_data("GunPoint")
print(metadata)
# {'n_train': 50, 'n_test': 150, 'length': 150, 'n_classes': 2, ...}
```

## 基准测试工具

### 加载已发布的结果

访问预先计算的基准测试结果：

```python
from aeon.benchmarking import get_estimator_results

# Get results for specific algorithm on dataset
results = get_estimator_results(
    estimator_name="ROCKET",
    dataset_name="GunPoint"
)

# Get all available estimators for a dataset
estimators = get_available_estimators("GunPoint")
```

### 重采样策略

创建可重复的训练/测试分割：

```python
from aeon.benchmarking import stratified_resample

# Stratified resampling maintaining class distribution
X_train, X_test, y_train, y_test = stratified_resample(
    X, y,
    random_state=42,
    test_size=0.3
)
```

### 性能指标

时间序列任务的专门指标：

**异常检测指标**：
```python
from aeon.benchmarking.metrics.anomaly_detection import (
    range_precision,
    range_recall,
    range_f_score,
    range_roc_auc_score
)

# Range-based metrics for window detection
precision = range_precision(y_true, y_pred, alpha=0.5)
recall = range_recall(y_true, y_pred, alpha=0.5)
f1 = range_f_score(y_true, y_pred, alpha=0.5)
auc = range_roc_auc_score(y_true, y_scores)
```

**聚类指标**：
```python
from aeon.benchmarking.metrics.clustering import clustering_accuracy

# Clustering accuracy with label matching
accuracy = clustering_accuracy(y_true, y_pred)
```

**细分指标**：
```python
from aeon.benchmarking.metrics.segmentation import (
    count_error,
    hausdorff_error
)

# Number of change points difference
count_err = count_error(y_true, y_pred)

# Maximum distance between predicted and true change points
hausdorff_err = hausdorff_error(y_true, y_pred)
```

### 统计测试

算法比较的事后分析：

```python
from aeon.benchmarking import (
    nemenyi_test,
    wilcoxon_test
)

# Nemenyi test for multiple algorithms
results = nemenyi_test(scores_matrix, alpha=0.05)

# Pairwise Wilcoxon signed-rank test
stat, p_value = wilcoxon_test(scores_alg1, scores_alg2)
```

## 基准集合

### UCR/UEA 时间序列档案

访问综合基准存储库：

```python
# Classification: 112 univariate + 30 multivariate datasets
X_train, y_train = load_classification("Chinatown", split="train")

# Automatically downloads from timeseriesclassification.com
```

### 莫纳什预测档案

```python
# Load forecasting datasets
y = load_forecasting("nn5_daily", return_X_y=False)
```

### 已发布的基准测试结果

主要比赛的预先计算结果：

- 2017 年单变量烘焙赛
- 2021 多元分类
- 2023 年单变量烘焙赛

## 工作流程示例

完整的基准测试工作流程：

```python
from aeon.datasets import load_classification
from aeon.classification.convolution_based import RocketClassifier
from aeon.benchmarking import get_estimator_results
from sklearn.metrics import accuracy_score
import numpy as np

# Load dataset
dataset_name = "GunPoint"
X_train, y_train = load_classification(dataset_name, split="train")
X_test, y_test = load_classification(dataset_name, split="test")

# Train model
clf = RocketClassifier(n_kernels=10000, random_state=42)
clf.fit(X_train, y_train)
y_pred = clf.predict(X_test)

# Evaluate
accuracy = accuracy_score(y_test, y_pred)
print(f"Accuracy: {accuracy:.4f}")

# Compare with published results
published = get_estimator_results("ROCKET", dataset_name)
print(f"Published ROCKET accuracy: {published['accuracy']:.4f}")
```

## 最佳实践

### 1. 使用标准分割

为了重现性，请使用提供的训练/测试分割：

```python
# Good: Use standard splits
X_train, y_train = load_classification("GunPoint", split="train")
X_test, y_test = load_classification("GunPoint", split="test")

# Avoid: Creating custom splits
X, y = load_classification("GunPoint")
X_train, X_test, y_train, y_test = train_test_split(X, y)
```

### 2.设置随机种子

确保再现性：

```python
clf = RocketClassifier(random_state=42)
results = stratified_resample(X, y, random_state=42)
```

### 3. 报告多个指标

不要依赖单一指标：

```python
from sklearn.metrics import accuracy_score, f1_score, precision_score

accuracy = accuracy_score(y_test, y_pred)
f1 = f1_score(y_test, y_pred, average='weighted')
precision = precision_score(y_test, y_pred, average='weighted')
```

### 4.交叉验证

对于小数据集的稳健评估：

```python
from sklearn.model_selection import cross_val_score

scores = cross_val_score(
    clf, X_train, y_train,
    cv=5,
    scoring='accuracy'
)
print(f"CV Accuracy: {scores.mean():.4f} (+/- {scores.std():.4f})")
```

### 5. 与基线比较

始终与简单基线进行比较：

```python
from aeon.classification.distance_based import KNeighborsTimeSeriesClassifier

# Simple baseline: 1-NN with Euclidean distance
baseline = KNeighborsTimeSeriesClassifier(n_neighbors=1, distance="euclidean")
baseline.fit(X_train, y_train)
baseline_acc = baseline.score(X_test, y_test)

print(f"Baseline: {baseline_acc:.4f}")
print(f"Your model: {accuracy:.4f}")
```

### 6. 统计意义

测试改进是否具有统计显着性：

```python
from aeon.benchmarking import wilcoxon_test

# Run on multiple datasets
accuracies_alg1 = [0.85, 0.92, 0.78, 0.88]
accuracies_alg2 = [0.83, 0.90, 0.76, 0.86]

stat, p_value = wilcoxon_test(accuracies_alg1, accuracies_alg2)
if p_value < 0.05:
    print("Difference is statistically significant")
```

## 数据集发现

查找符合条件的数据集：

```python
# List all available classification datasets
from aeon.datasets import get_available_datasets

datasets = get_available_datasets("classification")
print(f"Found {len(datasets)} classification datasets")

# Filter by properties
univariate_datasets = [
    d for d in datasets
    if get_dataset_meta_data(d)['n_channels'] == 1
]
```