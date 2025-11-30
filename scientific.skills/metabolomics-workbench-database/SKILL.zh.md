<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：代谢组学工作台数据库
描述：“通过 REST API 访问 NIH 代谢组学工作台（4,200 多项研究）。查询代谢物、RefMet 命名法、MS/NMR 数据、m/z 搜索、研究元数据，以进行代谢组学和生物标志物发现。”
---

# 代谢组学工作台数据库

## 概述

代谢组学工作台是一个由 NIH 共同基金资助的综合平台，托管在 UCSD，作为代谢组学研究数据的主要存储库。它提供了对 4,200 多项经过处理的研究（3,790 多项公开可用）的编程访问、通过 RefMet 标准化代谢物命名法以及跨多个分析平台（GC-MS、LC-MS、NMR）的强大搜索功能。

## 何时使用此技能

在查询代谢物结构、访问研究数据、标准化命名法、执行质谱搜索或通过代谢组学工作台 REST API 检索基因/蛋白质代谢物关联时，应使用此技能。

## 核心能力

### 1. 查询代谢物结构和数据

访问全面的代谢物信息，包括结构、标识符和外部数据库的交叉引用。

**关键操作：**
- 通过各种标识符检索化合物数据（PubChem CID、InChI Key、KEGG ID、HMDB ID 等）
- 将分子结构下载为 MOL 文件或 PNG 图像
- 访问标准化化合物分类
- 不同代谢物数据库之间的交叉引用

**查询示例：**
```python
import requests

# Get compound information by PubChem CID
response = requests.get('https://www.metabolomicsworkbench.org/rest/compound/pubchem_cid/5281365/all/json')

# Download molecular structure as PNG
response = requests.get('https://www.metabolomicsworkbench.org/rest/compound/regno/11/png')

# Get compound name by registry number
response = requests.get('https://www.metabolomicsworkbench.org/rest/compound/regno/11/name/json')
```

### 2. 访问研究元数据和实验结果

按各种标准查询代谢组学研究并检索完整的实验数据集。

**关键操作：**
- 按代谢物、机构、研究者或标题搜索研究
- 访问研究摘要、实验因素和分析详细信息
- 检索各种格式的完整实验数据
- 下载 mwTab 格式文件以获取完整的研究信息
- 查询非靶向代谢组学数据

**查询示例：**
<<<代码块_1>>>

### 3. 使用 RefMet 标准化代谢物命名法

使用 RefMet 数据库标准化代谢物名称并访问四个结构分辨率级别的系统分类。

**关键操作：**
- 将常见代谢物名称与标准化 RefMet 名称相匹配
- 通过化学式、精确质量或InChI键查询
- 访问层次分类（超类、主类、子类）
- 检索所有 RefMet 条目或按分类过滤

**查询示例：**
<<<代码块_2>>>

### 4. 执行质谱搜索

按质荷比 (m/z) 搜索具有指定离子加合物和耐受水平的化合物。

**关键操作：**
- 跨多个数据库搜索前体离子质量（Metabolomics Workbench、LIPIDS、RefMet）
- 指定离子加合物类型（M+H、M-H、M+Na、M+NH4、M+2H 等）
- 计算已知代谢物与特定加合物的精确质量
- 设置质量容差以实现灵活匹配

**查询示例：**
<<<代码块_3>>>

### 5. 通过分析和生物参数筛选研究

使用 MetStat 上下文查找匹配特定实验条件的研究。

**关键操作：**
- 通过分析方法（LCMS、GCMS、NMR）过滤
- 指定电离极性（正、负）
- 按色谱类型过滤（HILIC、RP、GC）
- 针对特定物种、样本来源或疾病
- 使用分号分隔格式组合多个过滤器

**查询示例：**
<<<代码块_4>>>

### 6. 访问基因和蛋白质信息

检索与代谢途径和代谢物代谢相关的基因和蛋白质数据。

**关键操作：**
- 通过符号、名称或 ID 查询基因
- 访问蛋白质序列和注释
- 基因 ID、RefSeq ID 和 UniProt ID 之间的交叉引用
- 检索基因代谢物关联

**查询示例：**
<<<代码块_5>>>

## 常见工作流程

### 工作流程 1：寻找特定代谢物的研究

要查找包含特定代谢物测量值的所有研究：

1. 首先使用RefMet标准化代谢物名称：
   <<<代码块_6>>>

2. 使用标准化名称搜索研究：
   ```python
   response = requests.get('https://www.metabolomicsworkbench.org/rest/study/refmet_name/Glucose/summary/json')
   ```

3. 从具体研究中检索实验数据：
   ```python
   response = requests.get('https://www.metabolomicsworkbench.org/rest/study/study_id/ST000001/data/json')
   ```

### 工作流程 2：从 MS 数据中识别化合物

要从质谱 m/z 值中识别潜在的化合物：

1. 使用适当的加合物和容差执行 m/z 搜索：
   ```python
   response = requests.get('https://www.metabolomicsworkbench.org/rest/moverz/MB/180.06/M+H/0.5/json')
   ```

2. 从结果中审查候选化合物
3. 检索候选化合物的详细信息：
   ```python
   response = requests.get('https://www.metabolomicsworkbench.org/rest/compound/regno/{regno}/all/json')
   ```

4. 下载结构进行确认：
   ```python
   response = requests.get('https://www.metabolomicsworkbench.org/rest/compound/regno/{regno}/png')
   ```

### 工作流程 3：探索疾病特异性代谢组学

要查找针对特定疾病和分析平台的代谢组学研究：

1. 使用 MetStat 筛选研究：
   ```python
   response = requests.get('https://www.metabolomicsworkbench.org/rest/metstat/LCMS;POSITIVE;;Human;;Cancer/json')
   ```

2. 检查结果中的研究 ID

3. 获取详细学习信息：
   ```python
   response = requests.get('https://www.metabolomicsworkbench.org/rest/study/study_id/ST{ID}/summary/json')
   ```

4. 检索完整的实验数据：
   ```python
   response = requests.get('https://www.metabolomicsworkbench.org/rest/study/study_id/ST{ID}/data/json')
   ```

## 输出格式

API 支持两种主要输出格式：
- **JSON**（默认）：机器可读格式，非常适合编程访问
- **TXT**：人类可读的制表符分隔文本格式

通过将 `/json` 或 `/txt` 附加到 API URL 来指定格式。省略format时，默认返回JSON。

## 最佳实践

1. **使用 RefMet 进行标准化**：在搜索研究之前始终通过 RefMet 标准化代谢物名称，以确保命名一致

2. **指定适当的加合物**：执行 m/z 搜索时，请使用适合您的分析方法的正确离子加合物类型（例如，M+H 用于正模式 ESI）

3. **设置合理的公差**：使用适当的质量公差值（通常低分辨率 MS 为 0.5 Da，高分辨率 MS 为 0.01 Da）

4. **缓存参考数据**：考虑缓存常用的参考数据（RefMet数据库、化合物信息）以最大程度地减少API调用

5. **处理分页**：对于大型结果集，准备好处理响应中的多个数据结构

6. **验证标识符**：尽可能跨多个数据库交叉引用代谢物标识符，以确保正确的化合物识别

## 资源

###参考资料/

详细的 API 参考文档可在 `references/api_reference.md` 中找到，包括：
- 完整的REST API端点规范
- 所有可用的上下文（化合物、研究、refmet、metstat、基因、蛋白质、moverz）
- 输入/输出参数详细信息
- 用于质谱分析的离子加合物类型
- 附加查询示例

当需要详细的 API 规范或使用不太常见的端点时，加载此参考文件。