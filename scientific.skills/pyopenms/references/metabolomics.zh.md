<!-- 此文件由机器翻译自 metabolomics.md -->

# 代谢组学工作流程

## 概述

PyOpenMS 提供用于非靶向代谢组学分析的专用工具，包括针对小分子优化的特征检测、加合物分组、化合物识别以及与代谢组学数据库的集成。

## 非靶向代谢组学管道

### 完整的工作流程

```python
import pyopenms as ms

def metabolomics_pipeline(input_files, output_dir):
    """
    Complete untargeted metabolomics workflow.

    Args:
        input_files: List of mzML file paths (one per sample)
        output_dir: Directory for output files
    """

    # Step 1: Peak picking and feature detection
    feature_maps = []

    for mzml_file in input_files:
        print(f"Processing {mzml_file}...")

        # Load data
        exp = ms.MSExperiment()
        ms.MzMLFile().load(mzml_file, exp)

        # Peak picking if needed
        if not exp.getSpectrum(0).isSorted():
            picker = ms.PeakPickerHiRes()
            exp_picked = ms.MSExperiment()
            picker.pickExperiment(exp, exp_picked)
            exp = exp_picked

        # Feature detection
        ff = ms.FeatureFinder()
        params = ff.getParameters("centroided")

        # Metabolomics-specific parameters
        params.setValue("mass_trace:mz_tolerance", 5.0)  # ppm, tighter for metabolites
        params.setValue("mass_trace:min_spectra", 5)
        params.setValue("isotopic_pattern:charge_low", 1)
        params.setValue("isotopic_pattern:charge_high", 2)  # Mostly singly charged

        features = ms.FeatureMap()
        ff.run("centroided", exp, features, params, ms.FeatureMap())

        features.setPrimaryMSRunPath([mzml_file.encode()])
        feature_maps.append(features)

        print(f"  Detected {features.size()} features")

    # Step 2: Adduct detection and grouping
    print("Detecting adducts...")
    adduct_grouped_maps = []

    adduct_detector = ms.MetaboliteAdductDecharger()
    params = adduct_detector.getParameters()
    params.setValue("potential_adducts", "[M+H]+,[M+Na]+,[M+K]+,[M+NH4]+,[M-H]-,[M+Cl]-")
    params.setValue("charge_min", 1)
    params.setValue("charge_max", 1)
    adduct_detector.setParameters(params)

    for fm in feature_maps:
        fm_out = ms.FeatureMap()
        adduct_detector.compute(fm, fm_out, ms.ConsensusMap())
        adduct_grouped_maps.append(fm_out)

    # Step 3: RT alignment
    print("Aligning retention times...")
    aligner = ms.MapAlignmentAlgorithmPoseClustering()

    params = aligner.getParameters()
    params.setValue("max_num_peaks_considered", 1000)
    params.setValue("pairfinder:distance_MZ:max_difference", 10.0)
    params.setValue("pairfinder:distance_MZ:unit", "ppm")
    aligner.setParameters(params)

    aligned_maps = []
    transformations = []
    aligner.align(adduct_grouped_maps, aligned_maps, transformations)

    # Step 4: Feature linking
    print("Linking features...")
    grouper = ms.FeatureGroupingAlgorithmQT()

    params = grouper.getParameters()
    params.setValue("distance_RT:max_difference", 60.0)  # seconds
    params.setValue("distance_MZ:max_difference", 5.0)  # ppm
    params.setValue("distance_MZ:unit", "ppm")
    grouper.setParameters(params)

    consensus_map = ms.ConsensusMap()
    grouper.group(aligned_maps, consensus_map)

    print(f"Created {consensus_map.size()} consensus features")

    # Step 5: Gap filling (fill missing values)
    print("Filling gaps...")
    # Gap filling not directly available in Python API
    # Would use TOPP tool FeatureFinderMetaboIdent

    # Step 6: Export results
    consensus_file = f"{output_dir}/consensus.consensusXML"
    ms.ConsensusXMLFile().store(consensus_file, consensus_map)

    # Export to CSV for downstream analysis
    df = consensus_map.get_df()
    csv_file = f"{output_dir}/metabolite_table.csv"
    df.to_csv(csv_file, index=False)

    print(f"Results saved to {output_dir}")

    return consensus_map

# Run pipeline
input_files = ["sample1.mzML", "sample2.mzML", "sample3.mzML"]
consensus = metabolomics_pipeline(input_files, "output")
```

## 加合物检测

### 配置加合物类型

<<<代码块_1>>>

### 访问加合物信息

<<<代码块_2>>>

## 化合物鉴定

### 基于质量的注释

<<<代码块_3>>>

### 基于 MS/MS 的识别

<<<代码块_4>>>

## 数据标准化

### 总离子流 (TIC) 归一化

<<<代码块_5>>>

## 质量控制

### 变异系数 (CV) 过滤

<<<代码块_6>>>

### 空白过滤

```python
# Remove features present in blank samples
blank_cols = [col for col in df.columns if 'Blank' in col]
sample_cols = [col for col in df.columns if 'Sample' in col]

if blank_cols and sample_cols:
    # Calculate mean intensity in blanks and samples
    blank_mean = df[blank_cols].mean(axis=1)
    sample_mean = df[sample_cols].mean(axis=1)

    # Keep features with 3x higher intensity in samples than blanks
    ratio = sample_mean / (blank_mean + 1)  # Add 1 to avoid division by zero
    filtered_df = df[ratio > 3]

    print(f"Features before blank filtering: {len(df)}")
    print(f"Features after blank filtering: {len(filtered_df)}")
```

## 缺失值插补

```python
import pandas as pd
import numpy as np

# Load data
df = consensus_map.get_df()

# Replace zeros with NaN
df = df.replace(0, np.nan)

# Count missing values
missing_per_feature = df.isnull().sum(axis=1)
print(f"Features with >50% missing: {sum(missing_per_feature > len(df.columns)/2)}")

# Simple imputation: replace with minimum value
for col in df.columns:
    if df[col].dtype in [np.float64, np.int64]:
        min_val = df[col].min() / 2  # Half minimum
        df[col].fillna(min_val, inplace=True)
```

## 代谢物表导出

### 创建分析就绪表

```python
import pandas as pd

def create_metabolite_table(consensus_map, output_file):
    """
    Create metabolite quantification table for statistical analysis.
    """

    # Get column headers (file descriptions)
    headers = consensus_map.getColumnHeaders()

    # Initialize data structure
    data = {
        'mz': [],
        'rt': [],
        'feature_id': []
    }

    # Add sample columns
    for map_idx, header in headers.items():
        sample_name = header.label or f"Sample_{map_idx}"
        data[sample_name] = []

    # Extract feature data
    for idx, cons_feature in enumerate(consensus_map):
        data['mz'].append(cons_feature.getMZ())
        data['rt'].append(cons_feature.getRT())
        data['feature_id'].append(f"F{idx:06d}")

        # Initialize intensities
        intensities = {map_idx: 0.0 for map_idx in headers.keys()}

        # Fill in measured intensities
        for handle in cons_feature.getFeatureList():
            map_idx = handle.getMapIndex()
            intensities[map_idx] = handle.getIntensity()

        # Add to data structure
        for map_idx, header in headers.items():
            sample_name = header.label or f"Sample_{map_idx}"
            data[sample_name].append(intensities[map_idx])

    # Create DataFrame
    df = pd.DataFrame(data)

    # Sort by RT
    df = df.sort_values('rt')

    # Save to CSV
    df.to_csv(output_file, index=False)

    print(f"Metabolite table with {len(df)} features saved to {output_file}")

    return df

# Create table
df = create_metabolite_table(consensus_map, "metabolite_table.csv")
```

## 与外部工具集成

### 导出 MetaboAnalyst

```python
def export_for_metaboanalyst(df, output_file):
    """
    Format data for MetaboAnalyst input.

    Requires sample names as columns, features as rows.
    """

    # Transpose DataFrame
    # Remove metadata columns
    sample_cols = [col for col in df.columns if col not in ['mz', 'rt', 'feature_id']]

    # Extract sample data
    sample_data = df[sample_cols]

    # Transpose (samples as rows, features as columns)
    df_transposed = sample_data.T

    # Add feature identifiers as column names
    df_transposed.columns = df['feature_id']

    # Save
    df_transposed.to_csv(output_file)

    print(f"MetaboAnalyst format saved to {output_file}")

# Export
export_for_metaboanalyst(df, "for_metaboanalyst.csv")
```

## 最佳实践

### 样本大小和重复

- 每 5-10 次进样包括 QC 样品（混合样品）
- 运行空白样品以识别污染
- 每组至少使用 3 个生物重复
- 随机化样品注射顺序

### 参数优化

混合 QC 样品的测试参数：

```python
# Test different mass trace parameters
mz_tolerances = [3.0, 5.0, 10.0]
min_spectra_values = [3, 5, 7]

for tol in mz_tolerances:
    for min_spec in min_spectra_values:
        ff = ms.FeatureFinder()
        params = ff.getParameters("centroided")
        params.setValue("mass_trace:mz_tolerance", tol)
        params.setValue("mass_trace:min_spectra", min_spec)

        features = ms.FeatureMap()
        ff.run("centroided", exp, features, params, ms.FeatureMap())

        print(f"tol={tol}, min_spec={min_spec}: {features.size()} features")
```

### 保留时间窗口

根据色谱法调整：

```python
# For 10-minute LC gradient
params.setValue("distance_RT:max_difference", 30.0)  # 30 seconds

# For 60-minute LC gradient
params.setValue("distance_RT:max_difference", 90.0)  # 90 seconds
```