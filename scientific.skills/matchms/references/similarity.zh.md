<!-- 此文件由机器翻译自 similarity.md -->

# Matchms 相似度函数参考

本文档提供了有关 matchms 中可用的所有相似性评分方法的详细信息。

## 概述

Matchms 提供了多种相似性函数来比较质谱。使用 `calculate_scores()` 计算参考光谱集合和查询光谱集合之间的成对相似性。

```python
from matchms import calculate_scores
from matchms.similarity import CosineGreedy

scores = calculate_scores(references=library_spectra,
                         queries=query_spectra,
                         similarity_function=CosineGreedy())
```

## 基于峰值的相似度函数

这些函数根据峰模式（m/z 和强度值）比较质谱。

### 余弦贪心

**描述**：使用快速贪婪匹配算法计算两个光谱之间的余弦相似度。峰在指定的容差内进行匹配，并根据匹配的峰强度计算相似性。

**何时使用**：
- 大型数据集的快速相似度计算
- 通用光谱匹配
- 当速度优先于数学上的最佳匹配时

**参数**：
- `tolerance`（浮点数，默认=0.1）：峰匹配的最大 m/z 差异（道尔顿）
- `mz_power`（浮点数，默认=0.0）：m/z 加权的指数（0 = 无加权）
- `intensity_power`（浮点数，默认=1.0）：强度加权的指数

**示例**：
<<<代码块_1>>>

**输出**：0.0 和 1.0 之间的相似性得分，加上匹配峰值的数量。

---

### 余弦匈牙利语

**描述**：使用匈牙利算法计算余弦相似度以实现最佳峰值匹配。提供数学上最佳的峰值分配，但比 CosineGreedy 慢。

**何时使用**：
- 当需要最佳峰值匹配时
- 高质量的参考库比较
- 需要可重复的、数学上严格的结果的研究

**参数**：
- `tolerance`（浮点型，默认=0.1）：峰匹配的最大 m/z 差异
- `mz_power`（浮点型，默认=0.0）：m/z 加权的指数
- `intensity_power`（浮点型，默认=1.0）：强度加权指数

**示例**：
<<<代码块_2>>>

**输出**：0.0 和 1.0 之间的最佳相似性得分，加上匹配的峰值。

**注意**：比 CosineGreedy 慢；用于较小的数据集或准确性至关重要的情况。

---

### 修正余弦

**描述**：通过考虑母体 m/z 差异来扩展余弦相似性。根据前体质量之间的差异应用质量转移后允许峰匹配。可用于比较相关化合物（同位素、加合物、类似物）的光谱。

**何时使用**：
- 比较不同前体质量的光谱
- 识别结构类似物或衍生物
- 交叉电离模式比较
- 当母体质量差异有意义时

**参数**：
- `tolerance`（浮点数，默认=0.1）：移位后峰匹配的最大 m/z 差异
- `mz_power`（浮点型，默认=0.0）：m/z 加权的指数
- `intensity_power`（浮点数，默认=1.0）：强度加权的指数

**示例**：
<<<代码块_3>>>

**要求**：两个光谱必须具有有效的前体_mz 元数据。

---

### 中性损失余弦

**描述**：根据中性损失模式而不是片段 m/z 值计算相似性。中性损失是通过从前体 m/z 中减去片段 m/z 得出的。对于识别具有相似碎片模式的化合物特别有用。

**何时使用**：
- 比较不同前体质量的碎片模式
- 识别具有相似中性损失特征的化合物
- 常规余弦评分的补充
- 代谢物鉴定和分类

**参数**：
- `tolerance` (float, default=0.1): 匹配的最大中性损失差
- `mz_power`（浮点型，默认=0.0）：损失值权重的指数
- `intensity_power`（浮点型，默认=1.0）：强度加权指数

**示例**：
<<<代码块_4>>>

**要求**：
- 两个光谱必须具有有效的前体_mz 元数据
- 在评分之前使用 `add_losses()` 过滤器计算中性损失

---

## 结构相似函数

这些函数比较分子结构而不是光谱峰。

### 指纹相似度

**描述**：计算源自化学结构（SMILES 或 InChI）的分子指纹之间的相似性。支持多种指纹类型和相似度指标。

**何时使用**：
- 没有光谱数据的结构相似性
- 结合结构和光谱相似性
- 在光谱匹配之前预过滤候选者
- 结构-活性关系研究
**参数**：
- `fingerprint_type`（str，默认=“daylight”）：指纹类型
  - `"daylight"`：日光指纹
  - `"morgan1"`、`"morgan2"`、`"morgan3"`：半径为 1、2 或 3 的摩根指纹
- `similarity_measure`（str，默认=“jaccard”）：相似度度量
  - `"jaccard"`：杰卡德索引（交集/并集）
  - `"dice"`: 骰子系数 (2 * 交集 / (size1 + size2))
  - `"cosine"`：余弦相似度

**示例**：
<<<代码块_5>>>

**要求**：
- Spectra 必须具有有效的 SMILES 或 InChI 元数据
- 使用`add_fingerprint()`过滤器来计算指纹
- 需要rdkit库

---

## 基于元数据的相似度函数

这些函数比较元数据字段而不是光谱或结构数据。

### 元数据匹配

**描述**：比较光谱之间用户定义的元数据字段。支持分类数据的精确匹配和数值数据的基于容差的匹配。

**何时使用**：
- 按实验条件过滤（碰撞能量、保留时间）
- 仪器特定匹配
- 将元数据约束与光谱相似性相结合
- 基于元数据的自定义过滤

**参数**：
- `field` (str)：要比较的元数据字段名称
- `matching_type` (str, default="exact"): 匹配方法
  - `"exact"`：精确的字符串/值匹配
  - `"difference"`：数值的绝对差
  - `"relative_difference"`：数值的相对差异
- `tolerance`（浮点型，可选）：数字匹配的最大差异

**示例（精确匹配）**：
<<<代码块_6>>>

**示例（数字匹配）**：
```python
# Match retention time within 0.5 minutes
similarity_func = MetadataMatch(field="retention_time",
                                matching_type="difference",
                                tolerance=0.5)
scores = calculate_scores(references, queries, similarity_func)
```

**输出**：精确匹配时返回 1.0（匹配）或 0.0（不匹配）。对于数值匹配，返回基于差异的相似度分数。

---

### PrecursorMzMatch

**描述**：基于前体 m/z 值的二元匹配。根据母体质量是否在指定容差范围内匹配，返回 True/False。

**何时使用**：
- 按前体质量预过滤光谱库
- 基于大众的快速候选人选择
- 与其他相似性度量相结合
- 同量异位化合物鉴定

**参数**：
- `tolerance`（浮点数，默认=0.1）：匹配的最大 m/z 差异
- `tolerance_type`（str，默认=“Dalton”）：公差单位
  - `"Dalton"`：绝对质量差
  - `"ppm"`：百万分之一（相对）

**示例**：
```python
from matchms.similarity import PrecursorMzMatch

# Match precursor within 0.1 Da
similarity_func = PrecursorMzMatch(tolerance=0.1, tolerance_type="Dalton")
scores = calculate_scores(references, queries, similarity_func)

# Match precursor within 10 ppm
similarity_func = PrecursorMzMatch(tolerance=10, tolerance_type="ppm")
scores = calculate_scores(references, queries, similarity_func)
```

**输出**：1.0（匹配）或0.0（不匹配）

**要求**：两个光谱必须具有有效的前体_mz 元数据。

---

### 家长质量匹配

**描述**：基于父质量（中性质量）值的二进制匹配。与 PrecursorMzMatch 类似，但使用计算的母体质量而不是母体 m/z。

**何时使用**：
- 比较不同电离模式的光谱
- 与加合物无关的匹配
- 基于中性质量的图书馆搜索

**参数**：
- `tolerance`（浮点数，默认=0.1）：匹配的最大质量差
- `tolerance_type`（str，默认=“道尔顿”）：公差单位（“道尔顿”或“ppm”）

**示例**：
```python
from matchms.similarity import ParentMassMatch

similarity_func = ParentMassMatch(tolerance=0.1, tolerance_type="Dalton")
scores = calculate_scores(references, queries, similarity_func)
```

**输出**：1.0（匹配）或0.0（不匹配）

**要求**：两个光谱必须具有有效的parent_mass元数据。

---

## 组合多个相似度函数

结合多个相似性指标以进行稳健的化合物识别：

```python
from matchms import calculate_scores
from matchms.similarity import CosineGreedy, ModifiedCosine, FingerprintSimilarity

# Calculate multiple similarity scores
cosine_scores = calculate_scores(refs, queries, CosineGreedy())
modified_cosine_scores = calculate_scores(refs, queries, ModifiedCosine())
fingerprint_scores = calculate_scores(refs, queries, FingerprintSimilarity())

# Combine scores with weights
for i, query in enumerate(queries):
    for j, ref in enumerate(refs):
        combined_score = (0.5 * cosine_scores.scores[j, i] +
                         0.3 * modified_cosine_scores.scores[j, i] +
                         0.2 * fingerprint_scores.scores[j, i])
```

## 访问分数结果

`Scores` 对象提供了多种方法来访问结果：

```python
# Get best matches for a query
best_matches = scores.scores_by_query(query_spectrum, sort=True)[:10]

# Get scores as numpy array
score_array = scores.scores

# Get scores as pandas DataFrame
import pandas as pd
df = scores.to_dataframe()

# Filter by threshold
high_scores = [(i, j, score) for i, j, score in scores.to_list() if score > 0.7]

# Save scores
scores.to_json("scores.json")
scores.to_pickle("scores.pkl")
```

## 性能考虑因素

**快速方法**（大型数据集）：
- 余弦贪心
- 前体MzMatch
- 家长质量匹配

**慢方法**（较小的数据集或较高的准确度）：
- 余弦匈牙利语
- ModifiedCosine（比 CosineGreedy 慢）
- 中性损失余弦
- FingerprintSimilarity（需要指纹计算）

**建议**：对于大规模库检索，请使用 PrecursorMzMatch 来预筛选候选者，然后应用 CosineGreedy 或 ModifiedCosine 来筛选结果。

## 常见相似性工作流程

### 标准库匹配
```python
from matchms.similarity import CosineGreedy

scores = calculate_scores(library_spectra, query_spectra,
                         CosineGreedy(tolerance=0.1))
```

### 多指标匹配
```python
from matchms.similarity import CosineGreedy, ModifiedCosine, FingerprintSimilarity

# Spectral similarity
cosine = calculate_scores(refs, queries, CosineGreedy())
modified = calculate_scores(refs, queries, ModifiedCosine())

# Structural similarity
fingerprint = calculate_scores(refs, queries, FingerprintSimilarity())
```

### 前体过滤匹配
```python
from matchms.similarity import PrecursorMzMatch, CosineGreedy

# First filter by precursor mass
mass_filter = calculate_scores(refs, queries, PrecursorMzMatch(tolerance=0.1))

# Then calculate cosine only for matching precursors
cosine_scores = calculate_scores(refs, queries, CosineGreedy())
```

## 进一步阅读

有关详细的 API 文档、参数说明和数学公式，请参阅：
https://matchms.readthedocs.io/en/latest/api/matchms.similarity.html