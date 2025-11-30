<!-- 此文件由机器翻译自 research_pipeline.md -->

# 研究管道 API 参考

## 核心课程

### 迪纳里奥

用于编排研究工作流程的主类。

#### 初始化

```python
from denario import Denario

den = Denario(project_dir="path/to/project")
```

**参数：**
- `project_dir` (str)：存储所有输出的研究项目目录的路径

#### 方法

##### set_data_description()

通过描述可用数据和分析工具来定义研究背景。

<<<代码块_1>>>

**参数：**
- `description` (str)：描述数据集、可用工具、研究领域和任何相关上下文的文本

**示例：**
<<<代码块_2>>>

**目的：** 通过提供有关哪些数据可用以及哪些分析可行的背景信息，为自动化创意生成奠定了基础。

##### get_idea()

根据数据描述生成研究假设。

<<<代码块_3>>>

**返回：** 研究想法/假设（内部存储在项目目录中）

**输出：** 创建一个包含生成的研究问题或假设的文件

**示例：**
<<<代码块_4>>>

##### set_idea()

手动指定一个研究想法，而不是生成一个。

<<<代码块_5>>>

**参数：**
- `idea` (str)：要调查的研究假设或问题

**示例：**
<<<代码块_6>>>

**用例：**当您有特定的研究方向并想要跳过自动生成想法时。

##### get_method()

根据想法和数据描述制定研究方法。

```python
den.get_method()
```

**返回：** 方法文档（内部存储在项目目录中）

**输出：** 创建结构化方法，包括：
- 分析方法
- 应用的统计方法
- 验证策略
- 预期产出

**示例：**
```python
den.get_method()
# Generates methodology: "Apply seasonal decomposition, compute correlation coefficients,
# perform statistical significance tests, generate visualization plots..."
```

##### set_method()

提供自定义方法而不是生成方法。

```python
den.set_method(method: str)
den.set_method(method: Path)  # Can also accept file paths
```

**参数：**
- `method`（str 或 Path）：方法描述或包含方法的 Markdown 文件路径

**示例：**
```python
# From string
den.set_method("""
1. Apply seasonal decomposition using STL
2. Compute Pearson correlation coefficients
3. Perform Mann-Kendall trend test
4. Generate time-series plots with confidence intervals
""")

# From file
den.set_method("methodology.md")
```

##### 获取结果()

执行方法、执行计算并生成结果。

```python
den.get_results()
```

**返回：** 带有分析输出的结果文档（内部存储在项目目录中）

**输出：** 创建结果，包括：
- 计算统计数据
- 生成的图形和可视化
- 数据表
- 分析结果

**示例：**
```python
den.get_results()
# Executes the methodology, runs analyses, creates plots, compiles findings
```

**注意：** 这是实际计算工作发生的地方。代理执行代码来执行方法中指定的分析。

##### set_results()

提供预先计算的结果而不是生成它们。

```python
den.set_results(results: str)
den.set_results(results: Path)  # Can also accept file paths
```

**参数：**
- `results`（str 或 Path）：结果描述或包含结果的 markdown 文件的路径

**示例：**
```python
# From string
den.set_results("""
Analysis Results:
- Correlation coefficient: 0.78 (p < 0.001)
- Seasonal amplitude: 5.2°C
- Long-term trend: +0.15°C per decade
- Figure 1: Seasonal decomposition (see attached)
""")

# From file
den.set_results("results.md")
```

**用例：** 在外部执行分析时或在不重新运行计算的情况下迭代论文写作时。

##### get_paper()

使用研究结果生成可供发表的 LaTeX 论文。

```python
den.get_paper(journal: Journal = None)
```

**参数：**
- `journal`（日志，可选）：用于格式化的目标日志。默认为通用格式。

**返回：** 具有正确格式的 LaTeX 纸张（存储在项目目录中）

**输出：** 创建：
- 完整的LaTeX源文件
- 编译的 PDF（如果 LaTeX 可用）
- 综合图表
- 格式正确的参考书目

**示例：**
```python
from denario import Journal

den.get_paper(journal=Journal.APS)
# Generates paper.tex and paper.pdf formatted for APS journals
```

### 期刊枚举

支持的日志格式的枚举。

```python
from denario import Journal
```

#### 可用期刊

- `Journal.APS` - 美国物理学会格式
  - 适用于物理评论、物理评论快报等。
  - 使用 RevTeX 文档类

可能还有其他期刊格式可用。查看最新的第纳尔文档以获取完整列表。

#### 用法

```python
from denario import Denario, Journal

den = Denario(project_dir="./research")
# ... complete workflow ...
den.get_paper(journal=Journal.APS)
```

## 工作流程模式

### 全自动管道

让第纳尔处理每个阶段：

```python
from denario import Denario, Journal

den = Denario(project_dir="./automated_research")

# Define context
den.set_data_description("""
Dataset: Sensor readings from IoT devices
Tools: pandas, numpy, sklearn, matplotlib
Goal: Anomaly detection in sensor networks
""")

# Automate entire pipeline
den.get_idea()        # Generate research idea
den.get_method()      # Develop methodology
den.get_results()     # Execute analysis
den.get_paper(journal=Journal.APS)  # Create paper
```

### 自定义创意，自动执行

提供您的研究问题，自动化其余部分：

```python
den = Denario(project_dir="./custom_idea")

den.set_data_description("Dataset: Financial time-series data...")

# Manual idea
den.set_idea("Investigate predictive models for stock market volatility using LSTM networks")

# Automated execution
den.get_method()
den.get_results()
den.get_paper(journal=Journal.APS)
```

### 完全手动生成模板

仅将第纳里奥用于纸张格式：

```python
den = Denario(project_dir="./manual_research")

# Provide everything manually
den.set_data_description("Pre-existing dataset description...")
den.set_idea("Pre-defined research hypothesis")
den.set_method("methodology.md")  # Load from file
den.set_results("results.md")      # Load from file

# Generate formatted paper
den.get_paper(journal=Journal.APS)
```

### 迭代细化

优化特定阶段而不重新运行所有内容：

```python
den = Denario(project_dir="./iterative")

# Initial run
den.set_data_description("Dataset description...")
den.get_idea()
den.get_method()
den.get_results()

# Refine methodology after reviewing results
den.set_method("""
Revised methodology:
- Use different statistical test
- Add sensitivity analysis
- Include cross-validation
""")

# Re-run only downstream stages
den.get_results()  # Re-execute with new method
den.get_paper(journal=Journal.APS)
```

## 项目目录结构

运行完整的工作流程后，项目目录包含：

```
project_dir/
├── data_description.txt    # Input: data context
├── idea.md                 # Generated or provided research idea
├── methodology.md          # Generated or provided methodology
├── results.md              # Generated or provided results
├── figures/                # Generated visualizations
│   ├── figure_1.png
│   ├── figure_2.png
│   └── ...
├── paper.tex               # Generated LaTeX source
├── paper.pdf               # Compiled PDF (if LaTeX available)
└── logs/                   # Agent execution logs
    └── ...
```

## 高级功能

### 多代理编排

Denario 使用 AG2 和 LangGraph 框架来协调多个专门代理：
- **创意代理**：根据数据描述生成研究假设
- **方法代理**：开发分析方法
- **执行代理**：运行计算并创建可视化
- **写作代理**：制作可供出版的手稿

这些代理自动协作，每个阶段都建立在先前的输出之上。

### 与科学工具集成

Denario 与常见的科学 Python 库集成：

- **pandas**：数据操作和分析
- **scikit-learn**：机器学习算法
- **scipy**：科学计算和统计
- **matplotlib/seaborn**：可视化
- **numpy**：数值运算

生成结果时，denario 可以使用这些库自动编写和执行代码。

### 再现性

所有阶段都会生成保存到项目目录的结构化输出：

- 版本控制友好（markdown 和 LaTeX）
- 可审计（代理决策和代码执行的日志）
- 可重复（保存的方法可以重新运行）

### 文献检索

Denario 包括文献检索功能，为研究思想和方法开发提供背景。请参阅`examples.md`了解文献检索工作流程。

## 错误处理

### 常见问题

**缺失数据说明：**
```python
den = Denario(project_dir="./project")
den.get_idea()  # Error: must call set_data_description() first
```

**解决方案：** 始终在产生想法之前设置数据描述。

**缺少先决条件阶段：**
```python
den = Denario(project_dir="./project")
den.get_results()  # Error: must have idea and method first
```

**解决方案：** 按照工作流程顺序或手动设置先决条件阶段。

**LaTeX 编译错误：**
```python
den.get_paper()  # May fail if LaTeX not installed
```

**解决方案：** 安装 LaTeX 发行版或使用预装 LaTeX 的 Docker 镜像。

## 最佳实践

### 数据描述质量

提供详细的背景信息以更好地产生创意：

```python
# Good: Detailed and specific
den.set_data_description("""
Dataset: 10 years of daily temperature readings from 50 weather stations
Format: CSV with columns [date, station_id, temperature, humidity]
Tools available: pandas, scipy, sklearn, matplotlib, seaborn
Domain: Climatology
Research interests: Climate change, seasonal patterns, regional variations
Known challenges: Missing data in 2015, station 23 has calibration issues
""")

# Bad: Too vague
den.set_data_description("Temperature data from weather stations")
```

### 方法验证

在执行之前检查生成的方法：

```python
den.get_method()
# Review the methodology.md file in project_dir
# If needed, refine with set_method()
```

### 增量开发

逐步建立研究渠道：

```python
# Stage 1: Validate idea generation
den.set_data_description("...")
den.get_idea()
# Review idea.md, adjust if needed

# Stage 2: Validate methodology
den.get_method()
# Review methodology.md, adjust if needed

# Stage 3: Execute and validate results
den.get_results()
# Review results.md and figures/

# Stage 4: Generate paper
den.get_paper(journal=Journal.APS)
```

### 版本控制集成

在项目目录中初始化git以进行跟踪：

```bash
cd project_dir
git init
git add .
git commit -m "Initial research workflow"
```

在每个阶段之后承诺跟踪研究的进展。