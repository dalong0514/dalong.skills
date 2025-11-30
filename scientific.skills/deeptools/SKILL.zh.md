<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：深度工具
描述：“NGS 分析工具包。BAM 到 bigWig 的转换、QC（相关性、PCA、指纹）、热图/概况（TSS、峰），用于 ChIP-seq、RNA-seq、ATAC-seq 可视化。”
---

# deepTools：NGS 数据分析工具包

## 概述

deepTools 是一套全面的 Python 命令行工具，旨在处理和分析高通量测序数据。使用 deepTools 为 ChIP-seq、RNA-seq、ATAC-seq、MNase-seq 和其他 NGS 实验执行质量控制、标准化数据、比较样本并生成出版质量的可视化结果。

**核心能力：**
- 将 BAM 对齐转换为标准化覆盖轨迹（bigWig/bedGraph）
- 质量控制评估（指纹、相关性、覆盖率）
- 样本比较及相关性分析
- 围绕基因组特征生成热图和剖面图
- 富集分析和峰区域可视化

## 何时使用此技能

该技能应该在以下情况下使用：

- **文件转换**：“将 BAM 转换为 bigWig”、“生成覆盖轨迹”、“标准化 ChIP-seq 数据”
- **质量控制**：“检查 ChIP 质量”、“比较重复”、“评估测序深度”、“QC 分析”
- **可视化**：“围绕 TSS 创建热图”、“绘制 ChIP 信号”、“可视化富集”、“生成剖面图”
- **样本比较**：“比较治疗与对照”、“关联样本”、“PCA 分析”
- **分析工作流程**：“分析 ChIP-seq 数据”、“RNA-seq 覆盖率”、“ATAC-seq 分析”、“完整工作流程”
- **使用特定文件类型**：基因组学上下文中的 BAM 文件、bigWig 文件、BED 区域文件

## 快速入门

对于 deepTools 的新用户，请从文件验证和常见工作流程开始：

### 1. 验证输入文件

在运行任何分析之前，使用验证脚本验证 BAM、bigWig 和 BED 文件：

```bash
python scripts/validate_files.py --bam sample1.bam sample2.bam --bed regions.bed
```

这会检查文件是否存在、BAM 索引和格式正确性。

### 2.生成工作流程模板

对于标准分析，请使用工作流程生成器创建自定义脚本：

<<<代码块_1>>>

### 3. 最常见的操作

有关常用命令和参数，请参阅 `assets/quick_reference.md`。

## 安装

<<<代码块_2>>>

## 核心工作流程

deepTools 工作流程通常遵循以下模式：**QC → 标准化 → 比较/可视化**

### ChIP-seq 质量控制工作流程

当用户请求 ChIP-seq QC 或质量评估时：

1. **使用 `scripts/workflow_generator.py chipseq_qc` 生成工作流程脚本**
2. **关键质量控制步骤**：
   - 样本相关性（multiBamSummary+plotCorrelation）
   - 主成分分析（plotPCA）
   - 覆盖率评估（plotCoverage）
   - 片段大小验证（bamPEFragmentSize）
   - ChIP 富集强度 (plotFingerprint)

**解释结果：**
- **相关性**：重复应以高相关性聚集在一起（>0.9）
- **指纹**：强 ChIP 显示急剧上升；平坦的对角线表示富集不良
- **覆盖率**：评估测序深度是否足以进行分析

完整工作流程详细信息位于 `references/workflows.md` →“ChIP-seq 质量控制工作流程”

### ChIP-seq 完整分析工作流程

对于从 BAM 到可视化的完整 ChIP-seq 分析：

1. **生成覆盖轨道**并进行归一化（bamCoverage）
2. **创建比较轨迹**（bamCompare for log2 比率）
3. **围绕特征计算信号矩阵**（computeMatrix）
4. **生成可视化**（plotHeatmap、plotProfile）
5. **峰处的富集分析** (plotEnrichment)

使用`scripts/workflow_generator.py chipseq_analysis`生成模板。

`references/workflows.md` 中的完整命令序列 →“ChIP-seq 分析工作流程”

### RNA-seq 覆盖工作流程

对于链特异性 RNA-seq 覆盖轨迹：

使用 bamCoverage 和 `--filterRNAstrand` 来分隔正向和反向链。

**重要提示：** 切勿使用 `--extendReads` 进行 RNA-seq（会延伸到剪接点）。

使用归一化：CPM 用于固定箱，RPKM 用于基因水平分析。

可用模板：`scripts/workflow_generator.py rnaseq_coverage`

详细信息请参见`references/workflows.md` →“RNA-seq 覆盖工作流程”

### ATAC-seq 分析工作流程

ATAC-seq 需要 Tn5 偏移校正：

1. **使用alignmentSieve 和 `--ATACshift` 进行移位读取**
2. **使用 bamCoverage 生成覆盖率**
3. **分析片段大小**（预计核小体阶梯图案）
4. **可视化峰值**（如果有）

模板：`scripts/workflow_generator.py atacseq`

`references/workflows.md` 中的完整工作流程 →“ATAC-seq 工作流程”

## 工具类别和常见任务

### BAM/bigWig 处理

**将 BAM 转换为标准化覆盖率：**
<<<代码块_3>>>

**比较两个样本（log2 比率）：**
<<<代码块_4>>>

**关键工具：** bamCoverage、bamCompare、multiBamSummary、multiBigwigSummary、 CorrectGCBias、alignmentSieve

完整参考：`references/tools_reference.md` →“BAM 和 bigWig 文件处理工具”

### 质量控制

**检查 ChIP 富集：**
<<<代码块_5>>>

**相关性示例：**
<<<代码块_6>>>

**关键工具：**plotFingerprint、plotCoverage、plotCorrelation、plotPCA、bamPEFragmentSize

完整参考：`references/tools_reference.md` →“质量控制工具”

### 可视化

**围绕 TSS 创建热图：**
```bash
# Compute matrix
computeMatrix reference-point -S signal.bw -R genes.bed \
    -b 3000 -a 3000 --referencePoint TSS -o matrix.gz

# Generate heatmap
plotHeatmap -m matrix.gz -o heatmap.png \
    --colorMap RdBu --kmeans 3
```

**创建剖面图：**
```bash
plotProfile -m matrix.gz -o profile.png \
    --plotType lines --colors blue red
```

**关键工具：**computeMatrix、plotHeatmap、plotProfile、plotEnrichment

完整参考：`references/tools_reference.md` →“可视化工具”

## 标准化方法

选择正确的标准化对于有效比较至关重要。请参阅`references/normalization_methods.md`以获得全面的指导。

**快速选型指南：**

- **ChIP-seq 覆盖率**：使用 RPGC 或 CPM
- **ChIP-seq 比较**：将 bamCompare 与 log2 和 readCount 结合使用
- **RNA-seq bins**：使用 CPM
- **RNA-seq 基因**：使用 RPKM（考虑基因长度）
- **ATAC-seq**：使用 RPGC 或 CPM

**标准化方法：**
- **RPGC**：1×基因组覆盖率（需要-- effectiveGenomeSize）
- **CPM**：每百万映射读取的计数
- **RPKM**：每百万每 kb 的读取次数（考虑区域长度）
- **BPM**：每百万个垃圾箱
- **无**：原始计数（不建议用于比较）

完整说明：`references/normalization_methods.md`

## 有效基因组大小

RPGC 标准化需要有效的基因组大小。共同价值观：

|有机体 |组装|尺寸|用途 |
|----------|----------|------|--------|
|人类 | GRCh38/hg38 | 2,913,022,398 | 2,913,022,398 `--effectiveGenomeSize 2913022398` |
|鼠标| GRCm38/mm10 | 2,652,783,500 | `--effectiveGenomeSize 2652783500` |
|斑马鱼 | GRCz11 | 1,368,780,147 | 1,368,780,147 `--effectiveGenomeSize 1368780147` |
| *果蝇* | DM6 | 142,573,017 | 142,573,017 `--effectiveGenomeSize 142573017` |
| *C。线虫* | ce10/ce11 | 100,286,401 | `--effectiveGenomeSize 100286401` |

具有特定读取长度值的完整表：`references/effective_genome_sizes.md`

## 跨工具的通用参数

许多 deepTools 命令共享这些选项：

**性能：**
- `--numberOfProcessors, -p`：启用并行处理（始终使用可用内核）
- `--region`：处理特定区域进行测试（例如，`chr1:1-1000000`）

**读取过滤：**
- `--ignoreDuplicates`：删除 PCR 重复项（建议用于大多数分析）
- `--minMappingQuality`：按对齐质量过滤（例如，`--minMappingQuality 10`）
- `--minFragmentLength` / `--maxFragmentLength`：片段长度范围
- `--samFlagInclude` / `--samFlagExclude`：SAM 标志过滤

**读取处理：**
- `--extendReads`：扩展到片段长度（ChIP-seq：是，RNA-seq：否）
- `--centerReads`：以片段中点为中心以获得更清晰的信号

## 最佳实践

### 文件验证
**始终首先验证文件**使用 `scripts/validate_files.py` 检查：
- 文件存在性和可读性
- 存在 BAM 索引（.bai 文件）
- BED格式正确性
- 文件大小合理

### 分析策略

1. **从 QC 开始**：在继续之前运行相关性、覆盖率和指纹分析
2. **小区域测试**：使用`--region chr1:1-10000000`进行参数测试
3. **文档命令**：保存完整的命令行以实现可重复性
4. **使用一致的归一化**：在比较样本中应用相同的方法
5. **验证基因组组装**：确保 BAM 和 BED 文件使用匹配的基因组构建

### ChIP-seq 特定

- **始终延长 ChIP-seq 的读取**：`--extendReads 200`
- **删除重复项**：在大多数情况下使用 `--ignoreDuplicates`
- **首先检查富集**：在详细分析之前运行plotFingerprint
- **GC 校正**：仅在检测到显着偏差时适用； GC 校正后切勿使用 `--ignoreDuplicates`

### RNA-seq 特异性

- **切勿延长 RNA-seq 的读取**（会跨越剪接点）
- **链特定**：对链库使用 `--filterRNAstrand forward/reverse`
- **标准化**：箱的 CPM，基因的 RPKM

### ATAC-seq 特定

- **应用Tn5校正**：使用alignmentSieve和`--ATACshift`
- **片段过滤**：设置适当的最小/最大片段长度
- **检查核小体模式**：片段大小图应显示阶梯模式

### 性能优化

1. **使用多个处理器**：`--numberOfProcessors 8`（或可用内核）
2. **增加 bin 大小**以加快处理速度并缩小文件
3. **针对内存有限的系统单独处理染色体**
4. **使用alignmentSieve预过滤BAM文件**以创建可重复使用的过滤文件
5. **使用 bigWig 而不是 bedGraph**：压缩且处理速度更快

## 故障排除

### 常见问题

**BAM 索引缺失：**
```bash
samtools index input.bam
```

**内存不足：**
使用 `--region` 单独处理染色体：
```bash
bamCoverage --bam input.bam -o chr1.bw --region chr1
```

**处理速度慢：**
增加 `--numberOfProcessors` 和/或增加 `--binSize`

**bigWig 文件太大：**
增加 bin 大小：`--binSize 50` 或更大

### 验证错误

运行验证脚本来识别问题：
```bash
python scripts/validate_files.py --bam *.bam --bed regions.bed
```

脚本输出中解释了常见错误和解决方案。

## 参考文档

该技能包括全面的参考文档：

### 参考文献/tools_reference.md
按类别组织的所有 deepTools 命令的完整文档：
- BAM 和 bigWig 处理工具（9 个工具）
- 质量控制工具（6 个工具）
- 可视化工具（3个工具）
- 杂项工具（2个工具）

每个工具包括：
- 目的和概述
- 关键参数及解释
- 使用示例
- 重要说明和最佳实践

**在以下情况下使用此参考：** 用户询问特定工具、参数或详细用法。

### 参考文献/workflows.md
常见分析的完整工作流程示例：
- ChIP-seq 质量控制工作流程
- ChIP-seq 完整的分析工作流程
- RNA-seq覆盖工作流程
- ATAC-seq 分析工作流程
- 多样本比较工作流程
- 峰区域分析工作流程
- 故障排除和性能提示

**在以下情况下使用此参考：** 用户需要完整的分析管道或工作流程示例。

### 参考文献/normalization_methods.md
标准化方法综合指南：
- 每种方法的详细解释（RPGC、CPM、RPKM、BPM等）
- 何时使用每种方法
- 公式和解释
- 按实验类型选择指南
- 常见的陷阱和解决方案
- 快速参考表

**在以下情况下使用此参考：** 用户询问标准化、比较样本或使用哪种方法。

### 参考文献/ effective_genome_sizes.md
有效基因组大小值和用途：
- 常见的生物体值（人类、小鼠、苍蝇、蠕虫、斑马鱼）
- 读取长度特定值
- 计算方法
- 何时以及如何在命令中使用
- 自定义基因组计算说明

**在以下情况下使用此参考：** 用户需要基因组大小来进行 RPGC 标准化或 GC 偏差校正。

## 帮助脚本

### 脚本/validate_files.py

验证 BAM、bigWig 和 BED 文件以进行 deepTools 分析。检查文件是否存在、索引和格式。

**用途：**
```bash
python scripts/validate_files.py --bam sample1.bam sample2.bam \
    --bed peaks.bed --bigwig signal.bw
```

**何时使用：** 在开始任何分析之前，或在排除错误时。

### 脚本/workflow_generator.py

为常见的 deepTools 工作流程生成可定制的 bash 脚本模板。

**可用的工作流程：**
- `chipseq_qc`：ChIP-seq 质量控制
- `chipseq_analysis`：完整的 ChIP-seq 分析
- `rnaseq_coverage`：链特异性 RNA-seq 覆盖率
- `atacseq`：带有 Tn5 校正的 ATAC-seq

**用途：**
```bash
# List workflows
python scripts/workflow_generator.py --list

# Generate workflow
python scripts/workflow_generator.py chipseq_qc -o qc.sh \
    --input-bam Input.bam --chip-bams "ChIP1.bam ChIP2.bam" \
    --genome-size 2913022398 --threads 8

# Run generated workflow
chmod +x qc.sh
./qc.sh
```

**何时使用：** 用户请求标准工作流程或需要模板脚本进行自定义。

## 资产

### 资产/quick_reference.md

快速参考卡，其中包含最常见的命令、有效的基因组大小和典型的工作流程模式。

**何时使用：** 用户需要快速的命令示例，而无需详细的文档。

## 处理用户请求

### 对于新用户

1.从安装验证开始
2. 使用 `scripts/validate_files.py` 验证输入文件
3. 根据实验类型推荐合适的工作流程
4.使用`scripts/workflow_generator.py`生成工作流程模板
5. 定制和执行指导

### 对于有经验的用户

1. 为请求的操作提供特定的工具命令
2. 参考`references/tools_reference.md`中的相应部分
3. 提出优化和最佳实践建议
4. 提供问题排查

### 对于特定任务

**“将 BAM 转换为 bigWig”：**
- 使用 bamCoverage 进行适当的标准化
- 根据用例推荐RPGC或CPM
- 为生物体提供有效的基因组大小
- 建议相关参数（extendReads、ignoreDuplicates、binSize）

**“检查 ChIP 质量”：**
- 运行完整的QC工作流程或专门使用plotFingerprint
- 解释结果的解释
- 根据结果建议后续行动

**“创建热图”：**
- 指导完成两步过程：computeMatrix→plotHeatmap
- 帮助选择适当的矩阵模式（参考点与比例区域）
- 建议可视化参数和聚类选项

**“比较样本”：**
- 推荐 bamCompare 进行两个样本比较
- 建议多个样本使用multiBamSummary +plotCorrelation
- 指导归一化方法选择

### 参考文档

当用户需要详细信息时：
- **工具详细信息**：直接到`references/tools_reference.md`中的特定部分
- **工作流程**：使用 `references/workflows.md` 进行完整的分析管道
- **标准化**：参考`references/normalization_methods.md`进行方法选择
- **基因组大小**：参考`references/effective_genome_sizes.md`

使用 grep 模式搜索参考文献：
```bash
# Find tool documentation
grep -A 20 "^### toolname" references/tools_reference.md

# Find workflow
grep -A 50 "^## Workflow Name" references/workflows.md

# Find normalization method
grep -A 15 "^### Method Name" references/normalization_methods.md
```

## 交互示例

**用户：“我需要分析我的 ChIP-seq 数据”**

应对办法：
1. 询问可用的文件（BAM 文件、峰、基因）
2. 使用验证脚本验证文件
3.生成chipseq_analysis工作流程模板
4.针对他们的特定文件和组织进行定制
5. 解释脚本运行时的每个步骤

**用户：“我应该使用哪种标准化？”**

应对办法：
1.询问实验类型（ChIP-seq、RNA-seq等）
2.询问比较目标（样本内或样本间）
3. 查阅`references/normalization_methods.md`选型指南
4. 有理由推荐适当的方法
5. 提供带有参数的命令示例

**用户：“围绕 TSS 创建热图”**

应对办法：
1.验证bigWig和gene BED文件是否可用
2. 在 TSS 处使用带有参考点模式的computeMatrix
3. 生成具有适当可视化参数的plotHeatmap
4. 如果数据集很大，建议聚类
5. 提供剖面图作为补充

## 重要提醒

- **首先进行文件验证**：在分析之前始终验证输入文件
- **标准化很重要**：为比较类型选择适当的方法
- **仔细扩展读取**：对于 ChIP-seq 是，对于 RNA-seq 否
- **使用所有核心**：将 `--numberOfProcessors` 设置为可用核心
- **区域测试**：使用`--region`进行参数测试
- **首先检查 QC**：在详细分析之前进行质量控制
- **记录一切**：保存命令以实现可重复性
- **参考文档**：使用全面的参考资料获取详细指导