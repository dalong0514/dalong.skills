<!-- 此文件由机器翻译自 signal_processing.md -->

# 信号处理

## 概述

PyOpenMS 提供用于处理原始质谱数据的算法，包括平滑、过滤、取峰、质心、归一化和反卷积。

## 算法模式

大多数信号处理算法都遵循标准模式：

```python
import pyopenms as ms

# 1. Create algorithm instance
algo = ms.AlgorithmName()

# 2. Get and modify parameters
params = algo.getParameters()
params.setValue("parameter_name", value)
algo.setParameters(params)

# 3. Apply to data
algo.filterExperiment(exp)  # or filterSpectrum(spec)
```

## 平滑

### 高斯滤波器

应用高斯平滑来减少噪声：

<<<代码块_1>>>

### Savitzky-Golay 滤波器

保留峰形状的多项式平滑：

<<<代码块_2>>>

## 峰值选取和质心

### 高分辨率峰值选取器

检测高分辨率数据中的峰值：

<<<代码块_3>>>

### CWT 峰值选取器

基于连续小波变换的取峰：

<<<代码块_4>>>

## 标准化

### 标准化器

标准化光谱内的峰强度：

<<<代码块_5>>>

## 峰值过滤

### 门槛割草机

删除低于强度阈值的峰值：

<<<代码块_6>>>

### 窗户割草机

仅在滑动窗口中保留最高峰：

```python
# Create window mower
window_mower = ms.WindowMower()

# Configure parameters
params = window_mower.getParameters()
params.setValue("windowsize", 50.0)  # Window size in m/z
params.setValue("peakcount", 2)  # Keep top N peaks per window
window_mower.setParameters(params)

# Apply filter
window_mower.filterExperiment(exp)
```

### N 个最大的山峰

仅保留 N 个最强烈的峰：

```python
# Create N largest filter
n_largest = ms.NLargest()

# Configure parameters
params = n_largest.getParameters()
params.setValue("n", 200)  # Keep 200 most intense peaks
n_largest.setParameters(params)

# Apply filter
n_largest.filterExperiment(exp)
```

## 基线减少

### 形态过滤器

使用形态学操作删除基线：

```python
# Create morphological filter
morph_filter = ms.MorphologicalFilter()

# Configure parameters
params = morph_filter.getParameters()
params.setValue("struc_elem_length", 3.0)  # Structuring element size
params.setValue("method", "tophat")  # Method: "tophat", "bothat", "erosion", "dilation"
morph_filter.setParameters(params)

# Apply filter
morph_filter.filterExperiment(exp)
```

## 频谱合并

### 光谱合并

将多个光谱合并为一个：

```python
# Create merger
merger = ms.SpectraMerger()

# Configure parameters
params = merger.getParameters()
params.setValue("average_gaussian:spectrum_type", "profile")
params.setValue("average_gaussian:rt_FWHM", 5.0)  # RT window
merger.setParameters(params)

# Merge spectra
merger.mergeSpectraBlockWise(exp)
```

## 反卷积

### 电荷反卷积

确定电荷状态并转换为中性质量：

```python
# Create feature deconvoluter
deconvoluter = ms.FeatureDeconvolution()

# Configure parameters
params = deconvoluter.getParameters()
params.setValue("charge_min", 1)
params.setValue("charge_max", 4)
params.setValue("potential_charge_states", "1,2,3,4")
deconvoluter.setParameters(params)

# Apply deconvolution
feature_map_out = ms.FeatureMap()
deconvoluter.compute(exp, feature_map, feature_map_out, ms.ConsensusMap())
```

### 同位素解卷积

删除同位素模式：

```python
# Create isotope wavelet transform
isotope_wavelet = ms.IsotopeWaveletTransform()

# Configure parameters
params = isotope_wavelet.getParameters()
params.setValue("max_charge", 3)
params.setValue("intensity_threshold", 10.0)
isotope_wavelet.setParameters(params)

# Apply transformation
isotope_wavelet.transform(exp)
```

## 保留时间对齐

### 地图对齐

调整多次运行的保留时间：

```python
# Create map aligner
aligner = ms.MapAlignmentAlgorithmPoseClustering()

# Load multiple experiments
exp1 = ms.MSExperiment()
exp2 = ms.MSExperiment()
ms.MzMLFile().load("run1.mzML", exp1)
ms.MzMLFile().load("run2.mzML", exp2)

# Create reference
reference = ms.MSExperiment()

# Align experiments
transformations = []
aligner.align(exp1, exp2, transformations)

# Apply transformation
transformer = ms.MapAlignmentTransformer()
transformer.transformRetentionTimes(exp2, transformations[0])
```

## 质量校准

### 内部校准

使用已知参考质量校准质量轴：

```python
# Create internal calibration
calibration = ms.InternalCalibration()

# Set reference masses
reference_masses = [500.0, 1000.0, 1500.0]  # Known m/z values

# Calibrate
calibration.calibrate(exp, reference_masses)
```

## 质量控制

### 频谱统计

计算质量指标：

```python
# Get spectrum
spec = exp.getSpectrum(0)

# Calculate statistics
mz, intensity = spec.get_peaks()

# Total ion current
tic = sum(intensity)

# Base peak
base_peak_intensity = max(intensity)
base_peak_mz = mz[intensity.argmax()]

print(f"TIC: {tic}")
print(f"Base peak: {base_peak_mz} m/z at {base_peak_intensity}")
```

## 频谱预处理管道

### 完整的预处理示例

```python
import pyopenms as ms

def preprocess_experiment(input_file, output_file):
    """Complete preprocessing pipeline."""

    # Load data
    exp = ms.MSExperiment()
    ms.MzMLFile().load(input_file, exp)

    # 1. Smooth with Gaussian filter
    gaussian = ms.GaussFilter()
    gaussian.filterExperiment(exp)

    # 2. Pick peaks
    picker = ms.PeakPickerHiRes()
    exp_picked = ms.MSExperiment()
    picker.pickExperiment(exp, exp_picked)

    # 3. Normalize intensities
    normalizer = ms.Normalizer()
    params = normalizer.getParameters()
    params.setValue("method", "to_TIC")
    normalizer.setParameters(params)
    normalizer.filterExperiment(exp_picked)

    # 4. Filter low-intensity peaks
    mower = ms.ThresholdMower()
    params = mower.getParameters()
    params.setValue("threshold", 10.0)
    mower.setParameters(params)
    mower.filterExperiment(exp_picked)

    # Save processed data
    ms.MzMLFile().store(output_file, exp_picked)

    return exp_picked

# Run pipeline
exp_processed = preprocess_experiment("raw_data.mzML", "processed_data.mzML")
```

## 最佳实践

### 参数优化

代表性数据的测试参数：

```python
# Try different Gaussian widths
widths = [0.1, 0.2, 0.5]

for width in widths:
    exp_test = ms.MSExperiment()
    ms.MzMLFile().load("test_data.mzML", exp_test)

    gaussian = ms.GaussFilter()
    params = gaussian.getParameters()
    params.setValue("gaussian_width", width)
    gaussian.setParameters(params)
    gaussian.filterExperiment(exp_test)

    # Evaluate quality
    # ... add evaluation code ...
```

### 保留原始数据

保留原始数据进行比较：

```python
# Load original
exp_original = ms.MSExperiment()
ms.MzMLFile().load("data.mzML", exp_original)

# Create copy for processing
exp_processed = ms.MSExperiment(exp_original)

# Process copy
gaussian = ms.GaussFilter()
gaussian.filterExperiment(exp_processed)

# Original remains unchanged
```

### 配置文件与质心数据

处理前检查数据类型：

```python
# Check if spectrum is centroided
spec = exp.getSpectrum(0)

if spec.isSorted():
    # Likely centroided
    print("Centroid data")
else:
    # Likely profile
    print("Profile data - apply peak picking")
```