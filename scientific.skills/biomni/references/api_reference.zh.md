<!-- 此文件由机器翻译自 api_reference.md -->

# Biomni API 参考

biomni 框架的综合 API 文档。

## A1特工级

A1 类是与 biomni 交互的主要接口。

### 初始化

```python
from biomni.agent import A1

agent = A1(
    path: str,              # Path to data lake directory
    llm: str,               # LLM model identifier
    verbose: bool = True,   # Enable verbose logging
    mcp_config: str = None  # Path to MCP server configuration
)
```

**参数：**

- **`path`** （str，必需）- biomni 数据湖的目录路径（~11GB）。如果数据不存在，则首次使用时会自动下载。

- **`llm`** （str，必需）- LLM 模型标识符。选项包括：
  - `'claude-sonnet-4-20250514'` - 建议平衡性能
  - `'claude-opus-4-20250514'` - 最大容量
  - `'gpt-4'`、`'gpt-4-turbo'` - OpenAI 模型
  - `'gemini-2.0-flash-exp'` - 谷歌双子座
  - `'llama-3.3-70b-versatile'` - 通过 Groq
  - 通过提供者配置自定义模型端点

- **`verbose`** （布尔值，可选，默认=True） - 启用代理推理、工具使用和代码执行的详细日志记录。

- **`mcp_config`** （str，可选） - 用于外部工具集成的 MCP（模型上下文协议）服务器配置文件的路径。

**示例：**
<<<代码块_1>>>

### 核心方法

#### `go(query: str) -> str`

自主执行生物医学研究任务。

<<<代码块_2>>>

**参数：**
- **`query`** (str) - 要执行的生物医学任务的自然语言描述

**退货：**
- **`str`** - 代理的最终答案或分析结果

**行为：**
1. 将查询分解为可执行的子任务
2.从综合数据库中检索相关知识
3.生成并执行Python代码进行分析
4. 迭代结果直到任务完成
5. 返回最终综合答案

**示例：**
<<<代码块_3>>>

#### `save_conversation_history(output_path: str, format: str = 'pdf')`

保存完整的对话历史记录，包括任务、推理、代码和结果。

<<<代码块_4>>>

**参数：**
- **`output_path`** (str) - 保存报告的文件路径
- **`format`**（str，可选，默认='pdf'） - 输出格式：`'pdf'`、`'html'` 或 `'markdown'`

**示例：**
<<<代码块_5>>>

#### `reset()`

重置代理状态并清除对话历史记录。

<<<代码块_6>>>

在开始新的独立任务时使用以清除先前的上下文。

**示例：**
```python
# Task 1
agent.go("Analyze dataset A")
agent.save_conversation_history("task1.pdf")

# Reset for fresh context
agent.reset()

# Task 2 - independent of Task 1
agent.go("Analyze dataset B")
```

### 通过default_config进行配置

全局配置参数可通过 `biomni.config.default_config` 访问。

```python
from biomni.config import default_config

# LLM Configuration
default_config.llm = "claude-sonnet-4-20250514"
default_config.llm_temperature = 0.7

# Execution Parameters
default_config.timeout_seconds = 1200  # 20 minutes
default_config.max_iterations = 50     # Max reasoning loops
default_config.max_tokens = 4096       # Max tokens per LLM call

# Code Execution
default_config.enable_code_execution = True
default_config.sandbox_mode = False    # Enable for restricted execution

# Data and Caching
default_config.data_cache_dir = "./biomni_cache"
default_config.enable_caching = True
```

**关键参数：**

- **`timeout_seconds`** (int, default=1200) - 任务执行的最长时间。增加复杂分析。

- **`max_iterations`** (int, default=50) - 最大代理推理循环。防止无限循环。

- **`enable_code_execution`** (bool, default=True) - 允许代理执行生成的代码。仅针对代码生成禁用。

- **`sandbox_mode`** (bool, default=False) - 启用沙盒代码执行（需要额外设置）。

## BiomniEval1 评估框架

用于在生物医学任务上对代理性能进行基准测试的框架。

### 初始化

```python
from biomni.eval import BiomniEval1

evaluator = BiomniEval1(
    dataset_path: str = None,  # Path to evaluation dataset
    metrics: list = None        # Evaluation metrics to compute
)
```

**示例：**
```python
evaluator = BiomniEval1()
```

### 方法

#### `evaluate(task_type: str, instance_id: str, answer: str) -> float`

根据真实情况评估代理的答案。

```python
score = evaluator.evaluate(
    task_type: str,     # Task category
    instance_id: str,   # Specific task instance
    answer: str         # Agent-generated answer
)
```

**参数：**
- **`task_type`** (str) - 任务类别：`'crispr_design'`、`'scrna_analysis'`、`'gwas_interpretation'`、`'drug_admet'`、`'clinical_diagnosis'`
- **`instance_id`** (str) - 数据集中任务实例的唯一标识符
- **`answer`** (str) - 代理的评估答案

**退货：**
- **`float`** - 评估分数（0.0 到 1.0）

**示例：**
```python
# Generate answer
result = agent.go("Design CRISPR screen for autophagy genes")

# Evaluate
score = evaluator.evaluate(
    task_type='crispr_design',
    instance_id='autophagy_001',
    answer=result
)
print(f"Score: {score:.2f}")
```

#### `load_dataset() -> dict`

加载 Biomni-Eval1 基准数据集。

```python
dataset = evaluator.load_dataset()
```

**退货：**
- **`dict`** - 包含按任务类型组织的任务实例的字典

**示例：**
```python
dataset = evaluator.load_dataset()

for task_type, instances in dataset.items():
    print(f"{task_type}: {len(instances)} instances")
```

#### `run_benchmark(agent: A1, task_types: list = None) -> dict`

对代理运行完整的基准评估。

```python
results = evaluator.run_benchmark(
    agent: A1,
    task_types: list = None  # Specific task types or None for all
)
```

**退货：**
- **`dict`** - 结果包含每个任务的分数、时间和详细指标

**示例：**
```python
results = evaluator.run_benchmark(
    agent=agent,
    task_types=['crispr_design', 'scrna_analysis']
)

print(f"Overall accuracy: {results['mean_score']:.2f}")
print(f"Average time: {results['mean_time']:.1f}s")
```

## 数据湖API

以编程方式访问集成生物医学数据库。

### 基因数据库查询

```python
from biomni.data import GeneDB

gene_db = GeneDB(path='./biomni_data')

# Query gene information
gene_info = gene_db.get_gene('BRCA1')
# Returns: {'symbol': 'BRCA1', 'name': '...', 'function': '...', ...}

# Search genes by pathway
pathway_genes = gene_db.search_by_pathway('DNA repair')
# Returns: List of gene symbols in pathway

# Get gene interactions
interactions = gene_db.get_interactions('TP53')
# Returns: List of interacting genes with interaction types
```

### 蛋白质结构获取

```python
from biomni.data import ProteinDB

protein_db = ProteinDB(path='./biomni_data')

# Get AlphaFold structure
structure = protein_db.get_structure('P38398')  # BRCA1 UniProt ID
# Returns: Path to PDB file or structure object

# Search PDB database
pdb_entries = protein_db.search_pdb('kinase', resolution_max=2.5)
# Returns: List of PDB IDs matching criteria
```

### 临床数据访问

```python
from biomni.data import ClinicalDB

clinical_db = ClinicalDB(path='./biomni_data')

# Query ClinVar variants
variant_info = clinical_db.get_variant('rs429358')  # APOE4 variant
# Returns: {'significance': '...', 'disease': '...', 'frequency': ...}

# Search OMIM for disease
disease_info = clinical_db.search_omim('Alzheimer')
# Returns: List of OMIM entries with gene associations
```

### 文献检索

```python
from biomni.data import LiteratureDB

lit_db = LiteratureDB(path='./biomni_data')

# Search PubMed abstracts
papers = lit_db.search('CRISPR screening cancer', max_results=10)
# Returns: List of paper dictionaries with titles, abstracts, PMIDs

# Get citations for paper
citations = lit_db.get_citations('PMID:12345678')
# Returns: List of citing papers
```

## MCP 服务器集成

通过模型上下文协议使用外部工具扩展 biomni。

### 配置格式

创建`.biomni/mcp_config.json`：

```json
{
  "servers": {
    "fda-drugs": {
      "command": "python",
      "args": ["-m", "mcp_server_fda"],
      "env": {
        "FDA_API_KEY": "${FDA_API_KEY}"
      }
    },
    "web-search": {
      "command": "npx",
      "args": ["-y", "@modelcontextprotocol/server-brave-search"],
      "env": {
        "BRAVE_API_KEY": "${BRAVE_API_KEY}"
      }
    }
  }
}
```

### 在任务中使用 MCP 工具
```python
# Initialize with MCP config
agent = A1(
    path='./data',
    llm='claude-sonnet-4-20250514',
    mcp_config='./.biomni/mcp_config.json'
)

# Agent can now use MCP tools automatically
result = agent.go("""
Search for FDA-approved drugs targeting EGFR.
Get their approval dates and indications.
""")
# Agent uses fda-drugs MCP server automatically
```

## 错误处理

常见异常情况及处理策略：

```python
from biomni.exceptions import (
    BiomniException,
    LLMError,
    CodeExecutionError,
    DataNotFoundError,
    TimeoutError
)

try:
    result = agent.go("Complex biomedical task")
except TimeoutError:
    # Task exceeded timeout_seconds
    print("Task timed out. Consider increasing timeout.")
    default_config.timeout_seconds = 3600
except CodeExecutionError as e:
    # Generated code failed to execute
    print(f"Code execution error: {e}")
    # Review generated code in conversation history
except DataNotFoundError:
    # Required data not in data lake
    print("Data not found. Ensure data lake is downloaded.")
except LLMError as e:
    # LLM API error
    print(f"LLM error: {e}")
    # Check API keys and rate limits
```

## 最佳实践

### 高效的 API 使用

1. **对相关任务重用代理实例**以维护上下文
2. **根据任务复杂程度设置适当的超时**
3. **使用缓存**避免冗余数据下载
4. **监控迭代**以尽早检测推理循环

### 生产部署

```python
from biomni.agent import A1
from biomni.config import default_config
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)

# Production settings
default_config.timeout_seconds = 3600
default_config.max_iterations = 100
default_config.sandbox_mode = True  # Enable sandboxing

# Initialize with error handling
try:
    agent = A1(path='/data/biomni', llm='claude-sonnet-4-20250514')
    result = agent.go(task_query)
    agent.save_conversation_history(f'reports/{task_id}.pdf')
except Exception as e:
    logging.error(f"Task {task_id} failed: {e}")
    # Handle failure appropriately
```

### 内存管理

对于大规模分析：

```python
# Process datasets in chunks
chunk_results = []
for chunk in dataset_chunks:
    agent.reset()  # Clear memory between chunks
    result = agent.go(f"Analyze chunk: {chunk}")
    chunk_results.append(result)

# Combine results
final_result = combine_results(chunk_results)
```