<!-- 此文件由机器翻译自 basic_inference.md -->

# 使用 Arboreto 进行基本 GRN 推理

## 输入数据要求

Arboreto 需要以下两种格式之一的基因表达数据：

### Pandas DataFrame（推荐）
- **行**：观察结果（细胞、样本、条件）
- **列**：基因（以基因名称作为列标题）
- **格式**：数字表达式值

示例：
```python
import pandas as pd

# Load expression matrix with genes as columns
expression_matrix = pd.read_csv('expression_data.tsv', sep='\t')
# Columns: ['gene1', 'gene2', 'gene3', ...]
# Rows: observation data
```

### NumPy 数组
- **形状**：（观察、基因）
- **要求**：单独提供与列顺序匹配的基因名称列表

示例：
<<<代码块_1>>>

## 转录因子 (TF)

（可选）提供转录因子名称列表以限制监管推断：

<<<代码块_2>>>

如果未提供，则所有基因都被视为潜在的调节因子。

## 基本推理工作流程

### 使用 Pandas DataFrame

<<<代码块_3>>>

**关键**：需要 `if __name__ == '__main__':` 防护，因为 Dask 会在内部生成新进程。

### 使用 NumPy 数组

<<<代码块_4>>>

## 输出格式

Arboreto 返回一个包含三列的 Pandas DataFrame：

|专栏 |描述 |
|--------|-------------|
| `TF` |转录因子（调节因子）基因名称 |
| `target` |目标基因名称 |
| `importance` |监管重要性得分（越高=监管越强）|

输出示例：
<<<代码块_5>>>

## 设置随机种子

为了获得可重现的结果，请提供种子参数：

<<<代码块_6>>>

## 算法选择

大多数情况下使用 `grnboost2()`（速度更快，处理大型数据集）：
```python
from arboreto.algo import grnboost2
network = grnboost2(expression_data=expression_matrix)
```

使用 `genie3()` 进行比较或特定要求：
```python
from arboreto.algo import genie3
network = genie3(expression_data=expression_matrix)
```

详细算法比较请参见`references/algorithms.md`。