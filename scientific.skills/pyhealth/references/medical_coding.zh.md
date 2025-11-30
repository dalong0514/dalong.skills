<!-- 此文件由机器翻译自 medical_coding.md -->

# PyHealth 医疗代码翻译

## 概述

医疗保健数据使用多种编码系统和标准。 PyHealth 的 MedCode 模块通过本体查找和跨系统映射实现医疗编码系统之间的翻译和映射。

## 核心课程

### 内部映射
处理系统内本体查找和分层导航。

**关键能力：**
- 使用属性（名称、描述）进行代码查找
- 祖先/后代层次结构遍历
- 代码标准化和转换
- 亲子关系导航

### 交叉映射
管理不同编码标准之间的跨系统映射。

**关键能力：**
- 编码系统之间的翻译
- 多对多关系处理
- 层级规范（针对药物）
- 双向映射支持

## 支持的编码系统

### 诊断代码

**ICD-9-CM（国际疾病分类，第 9 版，临床修改）**
- 传统诊断编码系统
- 具有 3-5 位代码的层次结构
- 2015 年之前用于美国医疗保健
- 用法：`from pyhealth.medcode import InnerMap`
  - `icd9_map = InnerMap.load("ICD9CM")`

**ICD-10-CM（国际疾病分类，第十次修订版，临床修改）**
- 现行诊断编码标准
- 字母数字代码（3-7 个字符）
- 比 ICD-9 更细化
- 用法：`from pyhealth.medcode import InnerMap`
  - `icd10_map = InnerMap.load("ICD10CM")`

**CCSCM（ICD-CM 临床分类软件）**
- 将 ICD 代码分为有临床意义的类别
- 降低分析维度
- 单级和多级层次结构
- 用法：`from pyhealth.medcode import CrossMap`
  - `icd_to_ccs = CrossMap.load("ICD9CM", "CCSCM")`

### 过程代码

**ICD-9-PROC（ICD-9 程序代码）**
- 住院手术分类
- 3-4 位数字代码
- 旧系统（2015 年之前）
- 用法：`from pyhealth.medcode import InnerMap`
  - `icd9proc_map = InnerMap.load("ICD9PROC")`

**ICD-10-PROC（ICD-10 程序编码系统）**
- 当前的程序编码标准
- 7 个字符的字母数字代码
- 比 ICD-9-PROC 更详细
- 用法：`from pyhealth.medcode import InnerMap`
  - `icd10proc_map = InnerMap.load("ICD10PROC")`

**CCSPROC（程序临床分类软件）**
- 将过程代码分组
- 简化程序分析
- 用法：`from pyhealth.medcode import CrossMap`
  - `proc_to_ccs = CrossMap.load("ICD9PROC", "CCSPROC")`

### 药物代码

**NDC（国家药品代码）**
- 美国FDA药品识别系统
- 10 或 11 位代码
- 产品级特异性（制造商、强度、包装）
- 用法：`from pyhealth.medcode import InnerMap`
  - `ndc_map = InnerMap.load("NDC")`

**接收标准**
- 标准化药物术语
- 标准化药物名称和关系
- 链接多个药物词汇
- 用法：`from pyhealth.medcode import CrossMap`
  - `ndc_to_rxnorm = CrossMap.load("NDC", "RXNORM")`

**ATC（解剖治疗化学分类）**
- WHO药物分类系统
- 5 级层次结构：
  - **1 级**：解剖学主要组（1 个字母）
  - **2级**：治疗亚组（2位数字）
  - **3级**：药理学亚组（1个字母）
  - **4 级**：化学分组（1 个字母）
  - **5级**：化学物质（2位数字）
- 示例：“C03CA01”=呋塞米
  - C = 心血管系统
  - C03 = 利尿剂
  - C03C = 高天花板利尿剂
  - C03CA = 磺胺类药物
  - C03CA01 = 呋塞米

**用途：**
```python
from pyhealth.medcode import CrossMap
ndc_to_atc = CrossMap.load("NDC", "ATC")
atc_codes = ndc_to_atc.map("00074-3799-13", level=3)  # Get ATC level 3
```

## 常用操作

### 内部映射操作

**1.代码查找**
<<<代码块_1>>>

**2.祖先穿越**
<<<代码块_2>>>

**3.后代遍历**
<<<代码块_3>>>

**4.代码标准化**
<<<代码块_4>>>

### CrossMap 操作

**1.直接翻译**
<<<代码块_5>>>

**2.分层药物图谱**
<<<代码块_6>>>

**3.双向映射**
```python
# Map in either direction
rxnorm_to_ndc = CrossMap.load("RXNORM", "NDC")
ndc_codes = rxnorm_to_ndc.map("197381")  # Get all NDC codes for RxNorm
```

## 工作流程示例

### 示例 1：标准化和分组诊断
```python
from pyhealth.medcode import InnerMap, CrossMap

# Load maps
icd9_map = InnerMap.load("ICD9CM")
icd_to_ccs = CrossMap.load("ICD9CM", "CCSCM")

# Process diagnosis codes
raw_codes = ["4280", "428.0", "42800"]

standardized = [icd9_map.standardize(code) for code in raw_codes]
# All become "428.0"

ccs_categories = [icd_to_ccs.map(code)[0] for code in standardized]
# All map to CCS category "108" (Heart failure)
```

### 示例2：药品分类分析
```python
from pyhealth.medcode import CrossMap

# Map NDC to ATC for drug class analysis
ndc_to_atc = CrossMap.load("NDC", "ATC")

patient_drugs = ["00074-3799-13", "00074-7286-01", "00456-0765-01"]

# Get therapeutic subgroups (ATC level 2)
drug_classes = []
for ndc in patient_drugs:
    atc_codes = ndc_to_atc.map(ndc, level=2)
    if atc_codes:
        drug_classes.append(atc_codes[0])

# Analyze drug class distribution
```

### 示例 3：ICD-9 到 ICD-10 的迁移
```python
from pyhealth.medcode import CrossMap

# Load ICD-9 to ICD-10 mapping
icd9_to_icd10 = CrossMap.load("ICD9CM", "ICD10CM")

# Convert historical ICD-9 codes
icd9_code = "428.0"
icd10_codes = icd9_to_icd10.map(icd9_code)
# Returns: ["I50.9", "I50.1", ...]  # Multiple possible ICD-10 codes

# Handle one-to-many mappings
for icd10_code in icd10_codes:
    print(f"ICD-9 {icd9_code} -> ICD-10 {icd10_code}")
```

## 与数据集集成

医疗代码翻译与 PyHealth 数据集无缝集成：

```python
from pyhealth.datasets import MIMIC4Dataset
from pyhealth.medcode import CrossMap

# Load dataset
dataset = MIMIC4Dataset(root="/path/to/data")

# Load code mapping
icd_to_ccs = CrossMap.load("ICD10CM", "CCSCM")

# Process patient diagnoses
for patient in dataset.iter_patients():
    for visit in patient.visits:
        diagnosis_events = [e for e in visit.events if e.vocabulary == "ICD10CM"]

        for event in diagnosis_events:
            ccs_codes = icd_to_ccs.map(event.code)
            print(f"Diagnosis {event.code} -> CCS {ccs_codes}")
```

## 用例

### 临床研究
- 跨不同编码系统标准化诊断
- 群组识别的群组相关条件
- 协调不同标准的多地点研究

### 药物安全性分析
- 按治疗类别对药物进行分类
- 在班级层面识别药物之间的相互作用
- 分析多药用药模式

### 医疗保健分析
- 减少诊断/程序维度
- 创建有意义的临床类别
- 支持跨编码系统变更的纵向分析

### 机器学习
- 创建一致的特征表示
- 处理训练/测试数据中的词汇不匹配
- 生成分层嵌入

## 最佳实践

1. **在映射之前始终标准化代码**，以确保格式一致
2. **适当处理一对多映射**（某些代码映射到多个目标）
3. **在映射药物时明确指定 ATC 级别**以避免歧义
4. **使用CCS类别**来减少诊断/程序维度
5. **验证映射**，因为某些代码可能没有直接翻译
6. **文档代码版本**（ICD-9 与 ICD-10）以维护数据来源