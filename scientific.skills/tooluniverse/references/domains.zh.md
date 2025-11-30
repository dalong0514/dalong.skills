<!-- 此文件由机器翻译自 domains.md -->

# ToolUniverse 工具域和类别

## 概述

ToolUniverse 集成了多个研究领域的 600 多种科学工具。本文档按科学学科和用例对工具进行分类。

## 主要科学领域

### 生物信息学

**序列分析：**
- 序列比对和比较
- 多序列比对（MSA）
- BLAST 和同源性搜索
- 基序查找和模式匹配

**基因组学：**
- 基因表达分析
- RNA-seq数据处理
- 变体调用和注释
- 基因组组装和注释
- 拷贝数变异分析

**功能分析：**
- 基因本体（GO）丰富
- 通路分析（KEGG、Reactome）
- 基因集富集分析（GSEA）
- 蛋白质结构域分析

**示例工具：**
- GEO数据下载和分析
- DESeq2差异表达
- KEGG通路富集
- UniProt 序列检索
- VEP变异注释

### 化学信息学

**分子描述符：**
- 化学性质计算
- 分子指纹
- SMILES/InChI 转换
- 3D构象异构体生成

**药物发现：**
- 虚拟筛选
- 分子对接
- ADMET预测
- 药物相似性评估（Lipinski 的五法则）
- 毒性预测

**化学数据库：**
- PubChem 化合物搜索
- ChEMBL生物活性数据
- ZINC化合物库
- DrugBank 药品信息

**示例工具：**
- RDKit 分子描述符
- AutoDock分子对接
- ZINC文库筛选
- ChEMBL 目标化合物关联

### 结构生物学

**蛋白质结构：**
- AlphaFold结构预测
- PDB结构检索
- 结构对齐和比较
- 结合位点预测
- 蛋白质-蛋白质相互作用预测

**结构分析：**
- 二级结构预测
- 溶剂可及性计算
- 结构质量评估
- Ramachandran情节分析

**示例工具：**
- AlphaFold结构预测
- PDB结构下载
- Fpocket结合位点检测
- DSSP二级结构分配

### 蛋白质组学

**蛋白质分析：**
- 质谱数据分析
- 蛋白质鉴定
- 翻译后修饰分析
- 蛋白质定量

**蛋白质数据库：**
- UniProt 蛋白质信息
- STRING 蛋白质相互作用
- IntAct 交互数据库

**示例工具：**
- UniProt数据检索
- STRING互动网络
- 质谱峰分析

### 机器学习

**型号类型：**
- 分类模型
- 回归模型
- 聚类算法
- 神经网络
- 深度学习模型

**应用：**
- 预测建模
- 特征选择
- 降维
- 模式识别
- 生物标志物发现

**示例工具：**
- Scikit学习模型
- TensorFlow/PyTorch 模型
- XGBoost 预测器
- 随机森林分类器

### 医疗/临床

**疾病数据库：**
- OpenTargets 疾病-靶标关联
- OMIM 遗传性疾病
- ClinVar 致病变异
- DisGeNET 疾病基因关联

**临床数据：**
- 电子健康记录分析
- 临床试验数据
- 诊断工具
- 治疗建议

**示例工具：**
- OpenTargets疾病查询
- ClinVar变异分类
- OMIM疾病查询
- FDA药品批准数据

### 神经科学

**脑成像：**
- 功能磁共振成像数据分析
- 脑图谱绘制
- 连接性分析
- 神经影像管道

**神经数据：**
- 电生理学分析
- 尖峰列车分析
- 神经网络模拟

### 图像处理

**生物医学成像：**
- 显微镜图像分析
- 细胞分割
- 物体检测
- 图像增强
- 特征提取

**图像分析：**
- ImageJ/斐济工具
- CellProfiler 管道
- 深度学习分割

### 系统生物学

**网络分析：**
- 生物网络建设
- 网络拓扑分析
- 模块识别
- Hub基因鉴定

**建模：**
- 系统生物学模型
- 代谢网络建模
- 信号通路模拟

## 按用例划分的工具类别

###文学与知识

**文献检索：**
- PubMed 文章搜索
- 文章摘要
- 引文分析
- 知识提取

**知识库：**
- 本体查询（GO、DO、HPO）
- 数据库交叉引用
- 实体识别

### 数据访问

**公共存储库：**
- GEO（基因表达综合）
- SRA（序列读取存档）
- PDB（蛋白质数据库）
- ChEMBL（生物活性数据库）

**API访问：**
- RESTful API 客户端
- 数据库查询工具
- 批量数据检索
### 可视化

**情节生成：**
- 热图
- 火山地块
- 曼哈顿地块
- 网络图
- 分子结构

### 实用程序

**数据处理：**
- 格式转换
- 数据标准化
- 统计分析
- 质量控制

**工作流程管理：**
- 管道建设
- 任务编排
- 结果汇总

## 按领域查找工具

将特定于域的关键字与 Tool_Finder 结合使用：

```python
# Bioinformatics
tools = tu.run({
    "name": "Tool_Finder_Keyword",
    "arguments": {"description": "RNA-seq genomics", "limit": 10}
})

# Cheminformatics
tools = tu.run({
    "name": "Tool_Finder_Keyword",
    "arguments": {"description": "molecular docking SMILES", "limit": 10}
})

# Structural biology
tools = tu.run({
    "name": "Tool_Finder_Keyword",
    "arguments": {"description": "protein structure PDB", "limit": 10}
})

# Clinical
tools = tu.run({
    "name": "Tool_Finder_Keyword",
    "arguments": {"description": "disease clinical variants", "limit": 10}
})
```

## 跨域应用

许多科学问题需要来自多个领域的工具：

- **精准医学**：基因组学+临床+蛋白质组学
- **药物发现**：化学信息学 + 结构生物学 + 机器学习
- **癌症研究**：基因组学 + 通路 + 文献
- **神经退行性疾病**：基因组学 + 蛋白质组学 + 成像