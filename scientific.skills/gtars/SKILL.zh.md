<!-- 此文件由机器翻译自 SKILL.md -->

---
姓名：gtars
描述：使用 Rust 与 Python 绑定进行基因组区间分析的高性能工具包。在处理基因组区域、BED 文件、覆盖轨迹、重叠检测、ML 模型标记化或计算基因组学和机器学习应用程序中的片段分析时使用。
---

# Gtars：Rust 中的基因组工具和算法

## 概述

Gtars 是一个高性能 Rust 工具包，用于操作、分析和处理基因组区间数据。它提供了用于重叠检测、覆盖分析、机器学习标记化和参考序列管理的专用工具。

在处理以下情况时使用此技能：
- 基因组间隔文件（BED格式）
- 基因组区域之间的重叠检测
- 覆盖轨道生成（WIG、BigWig）
- 基因组机器学习预处理和标记化
- 单细胞基因组学中的片段分析
- 参考序列检索和验证

## 安装

### Python安装

安装 gtars Python 绑定：

```bash
uv uv pip install gtars
```

### CLI 安装

安装命令行工具（需要 Rust/Cargo）：

<<<代码块_1>>>

### Rust 库

添加到 Rust 项目的 Cargo.toml 中：

<<<代码块_2>>>

## 核心能力

Gtars 被组织成专门的模块，每个模块都专注于特定的基因组分析任务：

### 1. 重叠检测和 IGD 索引

使用集成基因组数据库 (IGD) 数据结构有效检测基因组间隔之间的重叠。

**何时使用：**
- 寻找重叠的监管要素
- 变异注释
- 比较 ChIP-seq 峰
- 识别共享的基因组特征

**简单示例：**
<<<代码块_3>>>

请参阅 `references/overlap.md` 以获取全面的重叠检测文档。

### 2. 覆盖轨道生成

使用 uniwig 模块从测序数据生成覆盖轨迹。

**何时使用：**
- ATAC-seq 可访问性配置文件
- ChIP-seq 覆盖可视化
- RNA-seq 读取覆盖率
- 差异覆盖分析

**简单示例：**
<<<代码块_4>>>

有关详细的覆盖率分析工作流程，请参阅 `references/coverage.md`。

### 3. 基因组标记化

将基因组区域转换为机器学习应用程序的离散标记，特别是基因组数据的深度学习模型。

**何时使用：**
- 基因组机器学习模型的预处理
- 与 geniml 库集成
- 创建位置编码
- 在基因组序列上训练 Transformer 模型

**简单示例：**
<<<代码块_5>>>

请参阅 `references/tokenizers.md` 以获取标记化文档。

### 4.参考序列管理

按照 GA4GH refget 协议处理参考基因组序列并计算摘要。

**何时使用：**
- 验证参考基因组的完整性
- 提取特定的基因组序列
- 计算序列摘要
- 交叉参考比较

**简单示例：**
<<<代码块_6>>>

有关参考序列操作，请参阅`references/refget.md`。

### 5. 片段处理

分割和分析片段文件，对于单细胞基因组数据特别有用。

**何时使用：**
- 处理单细胞 ATAC-seq 数据
- 通过细胞条形码分割片段
- 基于聚类的片段分析
- 片段质量控制

**简单示例：**
```bash
# Split fragments by clusters
gtars fragsplit cluster-split --input fragments.tsv --clusters clusters.txt --output-dir ./by_cluster/
```

有关片段处理命令，请参阅`references/cli.md`。

### 6. 片段评分

分数片段与参考数据集重叠。

**何时使用：**
- 评估片段富集
- 将实验数据与参考数据进行比较
- 质量指标计算
- 跨样本批量评分

**简单示例：**
```bash
# Score fragments against reference
gtars scoring score --fragments fragments.bed --reference reference.bed --output scores.txt
```

## 常见工作流程

### 工作流程 1：峰重叠分析

识别重叠的基因组特征：

```python
import gtars

# Load two region sets
peaks = gtars.RegionSet.from_bed("chip_peaks.bed")
promoters = gtars.RegionSet.from_bed("promoters.bed")

# Find overlaps
overlapping_peaks = peaks.filter_overlapping(promoters)

# Export results
overlapping_peaks.to_bed("peaks_in_promoters.bed")
```

### 工作流程 2：覆盖跟踪管道

生成可视化覆盖轨迹：

```bash
# Step 1: Generate coverage
gtars uniwig generate --input atac_fragments.bed --output coverage.wig --resolution 10

# Step 2: Convert to BigWig for genome browsers
gtars uniwig generate --input atac_fragments.bed --output coverage.bw --format bigwig
```

### 工作流程 3：ML 预处理

为机器学习准备基因组数据：

```python
from gtars.tokenizers import TreeTokenizer
import gtars

# Step 1: Load training regions
regions = gtars.RegionSet.from_bed("training_peaks.bed")

# Step 2: Create tokenizer
tokenizer = TreeTokenizer.from_bed_file("training_peaks.bed")

# Step 3: Tokenize regions
tokens = [tokenizer.tokenize(r.chromosome, r.start, r.end) for r in regions]

# Step 4: Use tokens in ML pipeline
# (integrate with geniml or custom models)
```

## Python 与 CLI 用法

**在以下情况下使用 Python API：**
- 与分析管道集成
- 需要程序控制
- 使用 NumPy/Pandas
- 构建自定义工作流程

**在以下情况下使用 CLI：**
- 快速一次性分析
- 外壳脚本
- 批处理文件
- 原型制作工作流程

## 参考文档

全面的模块文档：

- **`references/python-api.md`** - 包含 RegionSet 操作、NumPy 集成和数据导出的完整 Python API 参考
- **`references/overlap.md`>** - IGD 索引、重叠检测和集合操作
- **`references/coverage.md`** - 使用 uniwig 生成覆盖轨道
- **`references/tokenizers.md`>** - ML 应用程序的基因组标记化
- **`references/refget.md`** - 参考序列管理和摘要
- **`references/cli.md`** - 命令行界面完整参考

## 与 geniml 集成

Gtars 作为 geniml Python 包的基础，为机器学习工作流程提供核心基因组区间操作。在处理 geniml 相关任务时，使用 gtar 进行数据预处理和标记化。

## 性能特点

- **原生 Rust 性能**：快速执行，内存开销低
- **并行处理**：大型数据集的多线程操作
- **内存效率**：流和内存映射文件支持
- **零复制操作**：NumPy 集成与最少的数据复制

## 数据格式

Gtars 使用标准基因组格式：

- **BED**：基因组间隔（3 柱或扩展）
- **WIG/BigWig**：覆盖曲目
- **FASTA**：参考序列
- **片段 TSV**：带有条形码的单细胞片段文件

## 错误处理和调试

启用详细日志记录以进行故障排除：

```python
import gtars

# Enable debug logging
gtars.set_log_level("DEBUG")
```

```bash
# CLI verbose mode
gtars --verbose <command>
```