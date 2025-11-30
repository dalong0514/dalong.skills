<!-- 此文件由机器翻译自 SKILL.md -->

---
名称： 极地
描述：“快速 DataFrame 库 (Apache Arrow)。选择、过滤、group_by、联接、惰性求值、CSV/Parquet I/O、表达式 API，用于高性能数据分析工作流程。”
---

# 北极星

## 概述

Polars 是一个基于 Apache Arrow 构建的用于 Python 和 Rust 的闪电般快速的 DataFrame 库。使用 Polars 基于表达式的 API、惰性求值框架和高性能数据操作功能，实现高效数据处理、pandas 迁移和数据管道优化。

## 快速入门

### 安装和基本使用

安装极地：
```python
uv pip install polars
```

基本 DataFrame 创建和操作：
<<<代码块_1>>>

## 核心概念

### 表达式

表达式是 Polars 操作的基本构建块。它们描述数据的转换，并且可以组合、重用和优化。

**关键原则：**
- 使用 `pl.col("column_name")` 引用列
- 构建复杂转换的链式方法
- 表达式是惰性的，仅在上下文中执行（select、with_columns、filter、group_by）

**示例：**
<<<代码块_2>>>

### 懒惰与急切评估

**Eager (DataFrame):** 操作立即执行
<<<代码块_3>>>

**惰性（LazyFrame）：**操作构建查询计划，在执行前进行优化
<<<代码块_4>>>

**何时使用惰性：**
- 处理大型数据集
- 复杂的查询管道
- 当只需要某些列/行时
- 性能至关重要

**惰性求值的好处：**
- 自动查询优化
- 谓词下推
- 投影下推
- 并行执行

有关详细概念，请加载`references/core_concepts.md`。

## 常用操作

### 选择
选择和操作列：
<<<代码块_5>>>

### 过滤器
按条件过滤行：
<<<代码块_6>>>

### 带列
添加或修改列，同时保留现有列：
```python
# Add new columns
df.with_columns(
    age_plus_10=pl.col("age") + 10,
    name_upper=pl.col("name").str.to_uppercase()
)

# Parallel computation (all columns computed in parallel)
df.with_columns(
    pl.col("value") * 10,
    pl.col("value") * 100,
)
```

### 分组依据和聚合
分组数据和计算聚合：
```python
# Basic grouping
df.group_by("city").agg(
    pl.col("age").mean().alias("avg_age"),
    pl.len().alias("count")
)

# Multiple group keys
df.group_by("city", "department").agg(
    pl.col("salary").sum()
)

# Conditional aggregations
df.group_by("city").agg(
    (pl.col("age") > 30).sum().alias("over_30")
)
```

有关详细操作模式，请加载`references/operations.md`。

## 聚合和窗口函数

### 聚合函数
`group_by` 上下文中的常见聚合：
- `pl.len()` - 计算行数
- `pl.col("x").sum()` - 总和值
- `pl.col("x").mean()` - 平均值
- `pl.col("x").min()` / `pl.col("x").max()` - 极端
- `pl.first()` / `pl.last()` - 第一个/最后一个值

### 带有 `over()` 的窗口函数
在保留行数的同时应用聚合：
```python
# Add group statistics to each row
df.with_columns(
    avg_age_by_city=pl.col("age").mean().over("city"),
    rank_in_city=pl.col("salary").rank().over("city")
)

# Multiple grouping columns
df.with_columns(
    group_avg=pl.col("value").mean().over("category", "region")
)
```

**映射策略：**
- `group_to_rows`（默认）：保留原始行顺序
- `explode`：更快，但将行分组在一起
- `join`：创建列表列

## 数据输入/输出

### 支持的格式
Polars 支持读取和写入：
- CSV、Parquet、JSON、Excel
- 数据库（通过连接器）
- 云存储（S3、Azure、GCS）
- 谷歌大查询
- 多个/分区文件

### 常见 I/O 操作

**CSV：**
```python
# Eager
df = pl.read_csv("file.csv")
df.write_csv("output.csv")

# Lazy (preferred for large files)
lf = pl.scan_csv("file.csv")
result = lf.filter(...).select(...).collect()
```

**镶木地板（建议用于性能）：**
```python
df = pl.read_parquet("file.parquet")
df.write_parquet("output.parquet")
```

**JSON：**
```python
df = pl.read_json("file.json")
df.write_json("output.json")
```

要获取全面的 I/O 文档，请加载 `references/io_guide.md`。

## 转换

### 加入
合并数据框：
```python
# Inner join
df1.join(df2, on="id", how="inner")

# Left join
df1.join(df2, on="id", how="left")

# Join on different column names
df1.join(df2, left_on="user_id", right_on="id")
```

### 连接
堆栈数据帧：
```python
# Vertical (stack rows)
pl.concat([df1, df2], how="vertical")

# Horizontal (add columns)
pl.concat([df1, df2], how="horizontal")

# Diagonal (union with different schemas)
pl.concat([df1, df2], how="diagonal")
```

### 旋转和取消旋转
重塑数据：
```python
# Pivot (wide format)
df.pivot(values="sales", index="date", columns="product")

# Unpivot (long format)
df.unpivot(index="id", on=["col1", "col2"])
```

有关详细转换示例，请加载 `references/transformations.md`。

## 大熊猫迁徙

Polars 通过更干净的 API 比 pandas 提供了显着的性能改进。主要区别：

### 概念差异
- **无索引**：Polars 仅使用整数位置
- **严格输入**：无静默类型转换
- **惰性评估**：通过 LazyFrame 提供
- **默认并行**：操作自动并行

### 常用操作映射

|运营|熊猫 |极地 |
|------------|--------|--------|
|选择列 | `df["col"]` | `df.select("col")` |
|过滤| `df[df["col"] > 10]` | `df.filter(pl.col("col") > 10)` |
|添加栏目| `df.assign(x=...)` | `df.with_columns(x=...)` |
|分组依据 | `df.groupby("col").agg(...)` | `df.group_by("col").agg(...)` |
|窗口| `df.groupby("col").transform(...)` | `df.with_columns(...).over("col")` |

### 关键语法模式

**熊猫顺序（慢）：**
```python
df.assign(
    col_a=lambda df_: df_.value * 10,
    col_b=lambda df_: df_.value * 100
)
```

**极地平行（快速）：**
```python
df.with_columns(
    col_a=pl.col("value") * 10,
    col_b=pl.col("value") * 100,
)
```

如需全面的迁移指南，请加载`references/pandas_migration.md`。

## 最佳实践

### 性能优化

1. **对大型数据集使用惰性求值：**
   ```python
   lf = pl.scan_csv("large.csv")  # Don't use read_csv
   result = lf.filter(...).select(...).collect()
   ```

2. **避免在热路径中使用Python函数：**
   - 留在表达式 API 内以实现并行化
- 仅在必要时使用 `.map_elements()`
   - 更喜欢原生 Polars 操作

3. **对非常大的数据使用流式传输：**
   ```python
   lf.collect(streaming=True)
   ```

4. **尽早仅选择需要的列：**
   ```python
   # Good: Select columns early
   lf.select("col1", "col2").filter(...)

   # Bad: Filter on all columns first
   lf.filter(...).select("col1", "col2")
   ```

5. **使用适当的数据类型：**
   - 低基数字符串的分类
   - 适当的整数大小（i32 与 i64）
   - 时态数据的日期类型

### 表达模式

**条件操作：**
```python
pl.when(condition).then(value).otherwise(other_value)
```

**跨多列的列操作：**
```python
df.select(pl.col("^.*_value$") * 2)  # Regex pattern
```

**空处理：**
```python
pl.col("x").fill_null(0)
pl.col("x").is_null()
pl.col("x").drop_nulls()
```

如需其他最佳实践和模式，请加载 `references/best_practices.md`。

## 资源

该技能包括全面的参考文档：

###参考资料/
- `core_concepts.md` - 表达式、惰性求值和类型系统的详细解释
- `operations.md` - 所有常见操作的综合指南和示例
- `pandas_migration.md` - 从 pandas 到 Polars 的完整迁移指南
- `io_guide.md` - 所有支持格式的数据 I/O 操作
- `transformations.md` - 连接、串联、旋转和重塑操作
- `best_practices.md` - 性能优化技巧和常见模式

当用户需要有关特定主题的详细信息时，根据需要加载这些参考资料。