<!-- 此文件由机器翻译自 data_structures.md -->

# 核心数据结构

## 概述

PyOpenMS 使用带有 Python 绑定的 C++ 对象。了解这些核心数据结构对于有效的数据操作至关重要。

## 光谱和实验对象

### MS实验

完整 LC-MS 实验数据（光谱和色谱图）的容器。

```python
import pyopenms as ms

# Create experiment
exp = ms.MSExperiment()

# Load from file
ms.MzMLFile().load("data.mzML", exp)

# Access properties
print(f"Number of spectra: {exp.getNrSpectra()}")
print(f"Number of chromatograms: {exp.getNrChromatograms()}")

# Get RT range
rts = [spec.getRT() for spec in exp]
print(f"RT range: {min(rts):.1f} - {max(rts):.1f} seconds")

# Access individual spectrum
spec = exp.getSpectrum(0)

# Iterate through spectra
for spec in exp:
    if spec.getMSLevel() == 2:
        print(f"MS2 spectrum at RT {spec.getRT():.2f}")

# Get metadata
exp_settings = exp.getExperimentalSettings()
instrument = exp_settings.getInstrument()
print(f"Instrument: {instrument.getName()}")
```

### MSSpectrum

具有 m/z 和强度阵列的单独质谱。

<<<代码块_1>>>

### MS色谱图

色谱图（TIC、XIC 或 SRM 离子对）。

<<<代码块_2>>>

## 特征对象

### 功能

检测到的色谱峰具有 2D 空间范围 (RT-m/z)。

<<<代码块_3>>>

### 特征图

单次 LC-MS 运行的特征集合。

<<<代码块_4>>>

### 共识功能

跨多个样本链接的功能。

<<<代码块_5>>>

### 共识图

跨样本的共识特征的收集。

<<<代码块_6>>>

## 识别对象

### 肽鉴定

单一光谱的识别结果。

```python
# Load identifications
protein_ids = []
peptide_ids = []
ms.IdXMLFile().load("identifications.idXML", protein_ids, peptide_ids)

# Access peptide identification
peptide_id = peptide_ids[0]

# Spectrum metadata
print(f"RT: {peptide_id.getRT():.2f}")
print(f"m/z: {peptide_id.getMZ():.4f}")

# Identification metadata
print(f"Identifier: {peptide_id.getIdentifier()}")
print(f"Score type: {peptide_id.getScoreType()}")
print(f"Higher score better: {peptide_id.isHigherScoreBetter()}")

# Get peptide hits
hits = peptide_id.getHits()
print(f"Number of hits: {len(hits)}")

for hit in hits:
    print(f"  Sequence: {hit.getSequence().toString()}")
    print(f"  Score: {hit.getScore()}")
    print(f"  Charge: {hit.getCharge()}")
```

### 肽命中

单个肽与谱图匹配。

```python
# Access hit
hit = peptide_id.getHits()[0]

# Sequence information
sequence = hit.getSequence()
print(f"Sequence: {sequence.toString()}")
print(f"Mass: {sequence.getMonoWeight():.4f}")

# Score and rank
print(f"Score: {hit.getScore()}")
print(f"Rank: {hit.getRank()}")

# Charge state
print(f"Charge: {hit.getCharge()}")

# Protein accessions
accessions = hit.extractProteinAccessionsSet()
for acc in accessions:
    print(f"Protein: {acc.decode()}")

# Meta values (additional scores, errors)
if hit.metaValueExists("MS:1002252"):  # mass error
    mass_error = hit.getMetaValue("MS:1002252")
    print(f"Mass error: {mass_error:.4f} ppm")
```

### 蛋白质鉴定

蛋白质水平的识别信息。

```python
# Access protein identification
protein_id = protein_ids[0]

# Search engine info
print(f"Search engine: {protein_id.getSearchEngine()}")
print(f"Search engine version: {protein_id.getSearchEngineVersion()}")

# Search parameters
search_params = protein_id.getSearchParameters()
print(f"Database: {search_params.db}")
print(f"Enzyme: {search_params.digestion_enzyme.getName()}")
print(f"Missed cleavages: {search_params.missed_cleavages}")
print(f"Precursor tolerance: {search_params.precursor_mass_tolerance}")

# Protein hits
hits = protein_id.getHits()
for hit in hits:
    print(f"Accession: {hit.getAccession()}")
    print(f"Score: {hit.getScore()}")
    print(f"Coverage: {hit.getCoverage():.1f}%")
```

### 蛋白质命中

个体蛋白质鉴定。

```python
# Access protein hit
protein_hit = protein_id.getHits()[0]

# Protein information
print(f"Accession: {protein_hit.getAccession()}")
print(f"Description: {protein_hit.getDescription()}")
print(f"Sequence: {protein_hit.getSequence()}")

# Scoring
print(f"Score: {protein_hit.getScore()}")
print(f"Coverage: {protein_hit.getCoverage():.1f}%")

# Rank
print(f"Rank: {protein_hit.getRank()}")
```

## 序列对象

### AA序列

经过修饰的氨基酸序列。

```python
# Create sequence from string
seq = ms.AASequence.fromString("PEPTIDE")

# Basic properties
print(f"Sequence: {seq.toString()}")
print(f"Length: {seq.size()}")
print(f"Monoisotopic mass: {seq.getMonoWeight():.4f}")
print(f"Average mass: {seq.getAverageWeight():.4f}")

# Individual residues
for i in range(seq.size()):
    residue = seq.getResidue(i)
    print(f"Position {i}: {residue.getOneLetterCode()}")
    print(f"  Mass: {residue.getMonoWeight():.4f}")
    print(f"  Formula: {residue.getFormula().toString()}")

# Modified sequence
mod_seq = ms.AASequence.fromString("PEPTIDEM(Oxidation)K")
print(f"Modified: {mod_seq.isModified()}")

# Check modifications
for i in range(mod_seq.size()):
    residue = mod_seq.getResidue(i)
    if residue.isModified():
        print(f"Modification at {i}: {residue.getModificationName()}")

# N-terminal and C-terminal modifications
term_mod_seq = ms.AASequence.fromString("(Acetyl)PEPTIDE(Amidated)")
```

### 经验公式

分子式表示。

```python
# Create formula
formula = ms.EmpiricalFormula("C6H12O6")  # Glucose

# Properties
print(f"Formula: {formula.toString()}")
print(f"Monoisotopic mass: {formula.getMonoWeight():.4f}")
print(f"Average mass: {formula.getAverageWeight():.4f}")

# Element composition
print(f"Carbon atoms: {formula.getNumberOf(b'C')}")
print(f"Hydrogen atoms: {formula.getNumberOf(b'H')}")
print(f"Oxygen atoms: {formula.getNumberOf(b'O')}")

# Arithmetic operations
formula2 = ms.EmpiricalFormula("H2O")
combined = formula + formula2  # Add water
print(f"Combined: {combined.toString()}")
```

## 参数对象

### 参数

算法使用的通用参数容器。

```python
# Get algorithm parameters
algo = ms.GaussFilter()
params = algo.getParameters()

# List all parameters
for key in params.keys():
    value = params.getValue(key)
    print(f"{key}: {value}")

# Get specific parameter
gaussian_width = params.getValue("gaussian_width")
print(f"Gaussian width: {gaussian_width}")

# Set parameter
params.setValue("gaussian_width", 0.2)

# Apply modified parameters
algo.setParameters(params)

# Copy parameters
params_copy = ms.Param(params)
```

## 最佳实践

### 内存管理

```python
# For large files, use indexed access instead of full loading
indexed_mzml = ms.IndexedMzMLFileLoader()
indexed_mzml.load("large_file.mzML")

# Access specific spectrum without loading entire file
spec = indexed_mzml.getSpectrumById(100)
```

### 类型转换

```python
# Convert peak arrays to numpy
import numpy as np

mz, intensity = spec.get_peaks()
# These are already numpy arrays

# Can perform numpy operations
filtered_mz = mz[intensity > 1000]
```

### 对象复制

```python
# Create deep copy
exp_copy = ms.MSExperiment(exp)

# Modifications to copy don't affect original
```