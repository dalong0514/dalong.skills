<!-- 此文件由机器翻译自 datasets.md -->

# PyHealth 数据集和数据结构

## 核心数据结构

### 活动
个别医疗事件的特征包括：
- **代码**：医疗代码（诊断、药物、手术、实验室测试）
- **词汇**：编码系统（ICD-9-CM、NDC、LOINC 等）
- **时间戳**：事件发生时间
- **值**：数值（用于实验室、生命体征）
- **单位**：测量单位

### 病人
跨访问期间按时间顺序组织的事件集合。每个患者包含：
- **患者_id**：唯一标识符
- **birth_datetime**：出生日期
- **性别**：患者性别
- **种族**：患者种族
- **访问**：访问对象列表

### 访问
医疗保健遭遇包含：
- **visit_id**：唯一标识符
- **encounter_time**：访问时间戳
- **discharge_time**：放电时间戳
- **visit_type**：遭遇类型（住院、门诊、急诊）
- **活动**：本次访问期间的活动列表

## 基础数据集类

**关键方法：**
- `get_patient(patient_id)`：检索单个患者记录
- `iter_patients()`：迭代所有患者
- `stats()`：获取数据集统计信息（患者、就诊、事件）
- `set_task(task_fn)`：定义预测任务

## 可用数据集

### 电子健康记录 (EHR) 数据集

**MIMIC-III 数据集** (`MIMIC3Dataset`)
- 来自贝斯以色列女执事医疗中心的重症监护室数据
- 40,000 多名重症监护患者
- 诊断、手术、药物、实验室结果
- 用法：`from pyhealth.datasets import MIMIC3Dataset`

**MIMIC-IV 数据集** (`MIMIC4Dataset`)
- 包含 70,000 多名患者的更新版本
- 提高数据质量和覆盖范围
- 增强人口统计和临床细节
- 用法：`from pyhealth.datasets import MIMIC4Dataset`

**eICU 数据集** (`eICUDataset`)
- 多中心重症监护数据库
- 超过 200,000 家医院的 200,000+ 次入院治疗
- 跨设施的标准化 ICU 数据
- 用法：`from pyhealth.datasets import eICUDataset`

**OMOP 数据集** (`OMOPDataset`)
- 观察性医疗结果合作伙伴格式
- 标准化通用数据模型
- 跨医疗保健系统的互操作性
- 用法：`from pyhealth.datasets import OMOPDataset`

**EHRShot 数据集** (`EHRShotDataset`)
- 少量学习的基准数据集
- 专门用于测试模型泛化能力
- 用法：`from pyhealth.datasets import EHRShotDataset`

### 生理信号数据集

**睡眠脑电图数据集：**
- `SleepEDFDataset`：用于睡眠分期的睡眠 EDF 数据库
- `SHHSDataset`：睡眠心脏健康研究数据
- `ISRUCDataset`：ISRUC-Sleep 数据库

**天普大学脑电图数据集：**
- `TUEVDataset`：异常脑电图事件检测
- `TUABDataset`：异常/正常脑电图分类
- `TUSZDataset`：癫痫检测

**所有信号数据集支持：**
- 多通道脑电图信号
- 标准化采样率
- 专家注释
- 睡眠阶段或异常标签

### 医学影像数据集

**COVID-19 CXR 数据集** (`COVID19CXRDataset`)
- 用于 COVID-19 分类的胸部 X 射线图像
- 多类别标签（COVID-19、肺炎、正常）
- 用法：`from pyhealth.datasets import COVID19CXRDataset`

### 基于文本的数据集

**医学转录数据集** (`MedicalTranscriptionsDataset`)
- 临床记录和转录
- 医学专业分类
- 基于文本的预测任务
- 用法：`from pyhealth.datasets import MedicalTranscriptionsDataset`

**心脏病学数据集** (`CardiologyDataset`)
- 心脏病患者记录
- 心血管疾病预测
- 用法：`from pyhealth.datasets import CardiologyDataset`

### 预处理数据集

**MIMIC 提取数据集** (`MIMICExtractDataset`)
- 预先提取的 MIMIC 特征
- 即用型基准测试数据
- 减少预处理要求
- 用法：`from pyhealth.datasets import MIMICExtractDataset`

## SampleDataset 类

将原始数据集转换为特定于任务的格式化样本。

**目的：** 将患者级数据转换为模型就绪的输入/输出对

**关键属性：**
- `input_schema`：定义输入数据结构
- `output_schema`：定义目标标签/预测
- `samples`：已处理样本列表

**使用模式：**
```python
# After setting task on BaseDataset
sample_dataset = dataset.set_task(task_fn)
```

## 数据分割函数

**患者级别分割** (`split_by_patient`)
- 确保没有患者出现在多次分割中
- 防止数据泄露
- 推荐用于临床预测任务

**访问级别拆分** (`split_by_visit`)
- 按个人访问划分
- 允许同一患者跨分区（谨慎使用）

**样本级分割** (`split_by_sample`)
- 随机样本分割
- 最灵活但可能导致泄漏

**参数：**
- `dataset`：要拆分的 SampleDataset
- `ratios`：分割比率元组（例如，[0.7, 0.1, 0.2]）
- `seed`：用于再现性的随机种子

## 通用工作流程

<<<代码块_1>>>

## 性能说明

- PyHealth 在医疗保健数据处理方面**比 pandas 快 3 倍**
- 针对大规模 EHR 数据集进行了优化
- 内存高效的患者迭代
- 用于特征提取的矢量化操作