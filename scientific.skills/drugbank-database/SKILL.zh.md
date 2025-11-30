<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：药物银行数据库
描述：访问并分析 DrugBank 数据库中的全面药物信息，包括药物特性、相互作用、靶标、途径、化学结构和药理学数据。在处理药物数据、药物发现研究、药理学研究、药物相互作用分析、靶点识别、化学相似性搜索、ADMET 预测或任何需要 DrugBank 中详细药物和药物靶点信息的任务时，应使用此技能。
---

# DrugBank 数据库

## 概述

DrugBank 是一个综合性生物信息学和化学信息学数据库，包含药物和药物靶点的详细信息。该技能支持以编程方式访问 DrugBank 数据，包括约 9,591 个药物条目（2,037 个 FDA 批准的小分子、241 种生物技术药物、96 种营养保健品和 6,000 多种实验化合物），每个条目有 200 多个数据字段。

## 核心能力

### 1. 数据访问和身份验证

使用 Python 通过适当的身份验证下载和访问 DrugBank 数据。该技能提供以下方面的指导：

- 安装和配置`drugbank-downloader`包
- 通过环境变量或配置文件安全地管理凭据
- 下载特定或最新的数据库版本
- 高效地打开和解析XML数据
- 使用缓存数据来优化性能

**何时使用**：设置 DrugBank 访问、下载数据库更新、初始项目配置。

**参考**：请参阅 `references/data-access.md` 了解详细的身份验证、下载过程、API 访问、缓存策略和故障排除。

### 2.药品信息查询

从数据库中提取全面的药物信息，包括标识符、化学性质、药理学、临床数据以及与外部数据库的交叉引用。

**查询能力**：
- 按 DrugBank ID、名称、CAS 编号或关键字搜索
- 提取药品基本信息（名称、类型、描述、适应症）
- 检索化学性质（SMILES、InChI、分子式）
- 获取药理学数据（作用机制、药效学、ADME）
- 访问外部标识符（PubChem、ChEMBL、UniProt、KEGG）
- 构建可搜索的药物数据集并导出到 DataFrames
- 按类型过滤药物（小分子、生物技术、营养保健品）

**何时使用**：检索特定药物信息、建立药物数据库、药理学研究、文献综述、药物分析。

**参考**：有关 XML 导航、查询函数、数据提取方法和性能优化，请参阅 `references/drug-queries.md`。

### 3. 药物-药物相互作用分析

分析药物相互作用 (DDI)，包括药物警戒和临床决策支持的机制、临床意义和相互作用网络。

**分析能力**：
- 提取特定药物的所有相互作用
- 建立双向互动网络
- 按严重性和机制对交互进行分类
- 检查药物对之间的相互作用
- 识别相互作用最多的药物
- 分析多药治疗方案的安全性
- 创建交互矩阵和网络图
- 在交互网络中执行社区检测
- 计算交互风险评分

**何时使用**：多药安全分析、临床决策支持、药物相互作用预测、药物警戒研究、识别禁忌症。

**参考**：交互提取、分类方法、网络分析和临床应用参见`references/interactions.md`。

### 4. 药物靶点和途径

获取有关药物-蛋白质相互作用的详细信息，包括靶标、酶、转运蛋白、载体和生物途径。

**目标分析能力**：
- 提取具有作用的药物靶点（抑制剂、激动剂、拮抗剂）
- 识别代谢酶（CYP450，II 期酶）
- 分析转运蛋白（摄取、流出）以进行 ADME 研究
- 将药物映射到生物途径（SMPDB）
- 寻找针对特定蛋白质的药物
- 识别具有共同目标的药物以进行重新利用
- 分析多药理学和脱靶效应
- 提取目标的基因本体（GO）术语
- 与 UniProt 交叉引用蛋白质数据

**何时使用**：作用机制研究、药物再利用研究、靶标识别、通路分析、预测脱靶效应、了解药物代谢。

**参考**：有关目标提取、通路分析、重新利用策略、CYP450 分析和转运蛋白分析，请参阅`references/targets-pathways.md`。
### 5. 化学性质和相似性

执行基于结构的分析，包括分子相似性搜索、属性计算、子结构搜索和 ADMET 预测。

**化学分析能力**：
- 提取化学结构（SMILES、InChI、分子式）
- 计算物理化学性质（MW、logP、PSA、氢键）
- 应用 Lipinski 的五规则和 Veber 的规则
- 计算分子之间的 Tanimoto 相似度
- 生成分子指纹（摩根、MACCS、拓扑）
- 使用 SMARTS 模式执行子结构搜索
- 寻找结构相似的药物进行重新利用
- 创建药物聚类的相似性矩阵
- 预测口服吸收和血脑屏障渗透性
- 使用 PCA 和聚类分析化学空间
- 导出化学性质数据库

**何时使用**：构效关系 (SAR) 研究、药物相似性搜索、QSAR 建模、药物相似性评估、ADMET 预测、化学空间探索。

**参考**：请参阅`references/chemical-analysis.md`了解结构提取、相似性计算、指纹生成、ADMET 预测和化学空间分析。

## 典型工作流程

### 药物发现工作流程
1. 使用`data-access.md`下载并访问最新的DrugBank数据
2. 使用`drug-queries.md`建立可检索的药品数据库
3. 使用`chemical-analysis.md`查找相似化合物
4. 使用`targets-pathways.md`来识别共享目标
5. 使用`interactions.md`检查候选组合的安全性

### 复方用药安全性分析
1. 使用`drug-queries.md`查找患者药物
2. 使用 `interactions.md` 检查所有成对交互
3. 使用`interactions.md`对交互严重性进行分类
4. 使用`interactions.md`计算总体风险评分
5.使用`targets-pathways.md`来理解交互机制

### 药物再利用研究
1. 使用`targets-pathways.md`查找具有共同目标的药物
2. 使用`chemical-analysis.md`查找结构相似的药物
3. 使用`drug-queries.md`提取适应症和药理学数据
4. 使用 `interactions.md` 评估潜在的联合疗法

### 药理学研究
1. 使用`drug-queries.md`提取感兴趣的药物
2. 使用 `targets-pathways.md` 识别所有蛋白质相互作用
3. 使用`targets-pathways.md`映射到生物途径
4. 使用 `chemical-analysis.md` 预测 ADMET 属性
5. 使用`interactions.md`来识别潜在的禁忌症

## 安装要求

### Python 包
```bash
uv pip install drugbank-downloader  # Core access
uv pip install bioversions          # Latest version detection
uv pip install lxml                 # XML parsing optimization
uv pip install pandas               # Data manipulation
uv pip install rdkit                # Chemical informatics (for similarity)
uv pip install networkx             # Network analysis (for interactions)
uv pip install scikit-learn         # ML/clustering (for chemical space)
```

### 帐户设置
1. 在 go.drugbank.com 创建免费帐户
2.接受许可协议（免费用于学术用途）
3. 获取用户名和密码凭证
4. 按照 `references/data-access.md` 中所述配置凭据

## 数据版本和再现性

始终指定 DrugBank 版本以进行可重复的研究：

<<<代码块_1>>>

记录出版物和分析脚本中使用的版本。

## 最佳实践

1. **凭据**：使用环境变量或配置文件，切勿硬编码
2. **版本控制**：指定准确的数据库版本以实现可重复性
3. **缓存**：缓存解析后的数据，避免重新下载和重新解析
4. **命名空间**：解析时正确处理XML命名空间
5. **验证**：使用前使用RDKit验证化学结构
6. **交叉引用**：使用外部标识符（UniProt、PubChem）进行集成
7. **临床背景**：在解释交互数据时始终考虑临床背景
8. **许可证合规性**：确保您的用例获得适当的许可

## 参考文档

所有详细的实施指南都组织在模块化参考文件中：

- **references/data-access.md**：身份验证、下载、解析、API 访问、缓存
- **references/drug-queries.md**：XML 导航、查询方法、数据提取、索引
- **references/interactions.md**：DDI提取、分类、网络分析、安全评分
- **references/targets-pathways.md**：目标/酶/转运蛋白提取、路径图谱、重新利用
- **references/chemical-analysis.md**：结构提取、相似性、指纹、ADMET 预测

根据您的具体分析要求按需加载这些参考文献。