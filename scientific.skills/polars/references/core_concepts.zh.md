<!-- 此文件由机器翻译自 core_concepts.md -->

# Polars 核心概念

## 表达式

表达式是 Polars API 的基础。它们是可组合的单元，描述数据转换而不立即执行它们。

### 什么是表达式？

表达式描述了数据的转换。它仅在特定上下文中具体化（执行）：
- `select()` - 选择并转换列
- `with_columns()` - 添加或修改列
- `filter()` - 过滤行
- `group_by().agg()` - 聚合数据

### 表达式语法

**基础栏目参考：**
```python
pl.col("column_name")
```

**计算表达式：**
<<<代码块_1>>>

### 表达式上下文

**选择上下文：**
<<<代码块_2>>>

**With_columns 上下文：**
<<<代码块_3>>>

**过滤上下文：**
<<<代码块_4>>>

**按上下文分组：**
<<<代码块_5>>>

### 表达式扩展

一次对多个列应用操作：

**所有栏目：**
<<<代码块_6>>>

**模式匹配：**
```python
# All columns ending with "_value"
df.select(pl.col("^.*_value$") * 100)

# All numeric columns
df.select(pl.col(pl.NUMERIC_DTYPES) + 1)
```

**排除模式：**
```python
df.select(pl.all().exclude("id", "name"))
```

### 表达式组成

表达式可以存储和重用：

```python
# Define reusable expressions
age_expression = pl.col("age") * 12
name_expression = pl.col("name").str.to_uppercase()

# Use in multiple contexts
df.select(age_expression, name_expression)
df.with_columns(age_months=age_expression)
```

## 数据类型

Polars 拥有基于 Apache Arrow 的严格类型系统。

### 核心数据类型

**数字：**
- `Int8`、`Int16`、`Int32`、`Int64` - 有符号整数
- `UInt8`、`UInt16`、`UInt32`、`UInt64` - 无符号整数
- `Float32`、`Float64` - 浮点数

**文字：**
- `Utf8` / `String` - UTF-8 编码字符串
- `Categorical` - 分类字符串（低基数）
- `Enum` - 固定字符串值集

**颞部：**
- `Date` - 日历日期（无时间）
- `Datetime` - 带有可选时区的日期和时间
- `Time` - 一天中的时间
- `Duration` - 持续时间/差异

**布尔值：**
- `Boolean` - 真/假值

**嵌套：**
- `List` - 可变长度列表
- `Array` - 固定长度数组
- `Struct` - 嵌套记录结构

**其他：**
- `Binary` - 二进制数据
- `Object` - Python 对象（避免在生产中使用）
- `Null` - 空类型

### 类型转换

显式地在类型之间进行转换：

```python
# Cast to different type
df.select(
    pl.col("age").cast(pl.Float64),
    pl.col("date_string").str.strptime(pl.Date, "%Y-%m-%d"),
    pl.col("id").cast(pl.Utf8)
)
```

### 空处理

Polars 在所有类型中使用一致的 null 处理：

**检查空值：**
```python
df.filter(pl.col("value").is_null())
df.filter(pl.col("value").is_not_null())
```

**填充空值：**
```python
pl.col("value").fill_null(0)
pl.col("value").fill_null(strategy="forward")
pl.col("value").fill_null(strategy="backward")
pl.col("value").fill_null(strategy="mean")
```

**删除空值：**
```python
df.drop_nulls()  # Drop any row with nulls
df.drop_nulls(subset=["col1", "col2"])  # Drop rows with nulls in specific columns
```

### 分类数据

对基数较低（重复值）的字符串列使用分类类型：

```python
# Cast to categorical
df.with_columns(
    pl.col("category").cast(pl.Categorical)
)

# Benefits:
# - Reduced memory usage
# - Faster grouping and joining
# - Maintains order information
```

## 懒惰与急切评估

Polars 支持两种执行模式：eager（DataFrame）和lazyFrame（LazyFrame）。

### 热切评估（数据帧）

操作立即执行：

```python
import polars as pl

# DataFrame operations execute right away
df = pl.read_csv("data.csv")  # Reads file immediately
result = df.filter(pl.col("age") > 25)  # Filters immediately
final = result.select("name", "age")  # Selects immediately
```

**何时使用 eager：**
- 适合内存的小数据集
- 笔记本中的交互式探索
- 简单的一次性操作
- 需要立即反馈

### 惰性评估 (LazyFrame)

操作构建一个查询计划，并在执行前进行优化：

```python
import polars as pl

# LazyFrame operations build a query plan
lf = pl.scan_csv("data.csv")  # Doesn't read yet
lf2 = lf.filter(pl.col("age") > 25)  # Adds to plan
lf3 = lf2.select("name", "age")  # Adds to plan
df = lf3.collect()  # NOW executes optimized plan
```

**何时使用惰性：**
- 大型数据集
- 复杂的查询管道
- 只需要数据的子集
- 性能至关重要
- 需要流媒体

### 查询优化

Polars 自动优化惰性查询：

**谓词下推：**
尽可能将筛选操作推送到数据源：
```python
# Only reads rows where age > 25 from CSV
lf = pl.scan_csv("data.csv")
result = lf.filter(pl.col("age") > 25).collect()
```

**投影下推：**
仅从数据源读取所需的列：
```python
# Only reads "name" and "age" columns from CSV
lf = pl.scan_csv("data.csv")
result = lf.select("name", "age").collect()
```

**查询计划检查：**
```python
# View the optimized query plan
lf = pl.scan_csv("data.csv")
result = lf.filter(pl.col("age") > 25).select("name", "age")
print(result.explain())  # Shows optimized plan
```

### 流媒体模式

处理大于内存的数据：

```python
# Enable streaming for very large datasets
lf = pl.scan_csv("very_large.csv")
result = lf.filter(pl.col("age") > 25).collect(streaming=True)
```

**流媒体的好处：**
- 处理大于RAM的数据
- 降低峰值内存使用量
- 基于块的处理
- 自动内存管理

**流媒体限制：**
- 并非所有操作都支持流式传输
- 对于小数据可能会更慢
- 某些操作需要具体化整个数据集

### 在渴望和懒惰之间转换

**渴望到懒惰：**
```python
df = pl.read_csv("data.csv")
lf = df.lazy()  # Convert to LazyFrame
```

**懒惰到渴望：**
```python
lf = pl.scan_csv("data.csv")
df = lf.collect()  # Execute and return DataFrame
```

## 内存格式

Polars 使用 Apache Arrow 列式内存格式：

**好处：**
- 与其他 Arrow 库共享零拷贝数据
- 高效的柱状作业
- SIMD矢量化
- 减少内存开销
- 快速序列化

**影响：**
- 数据按列存储，而不是按行存储
- 列操作速度非常快
- 随机行访问比 pandas 慢
- 最适合分析工作负载

## 并行化
Polars 使用 Rust 的并发性自动并行化操作：

**什么被并行化：**
- 组内的聚合
- 窗口功能
- 大多数表达评估
- 文件读取（多个文件）
- 加盟运营

**并行化要避免什么：**
- Python 用户定义函数 (UDF)
- `.map_elements()` 中的 Lambda 函数
- 顺序`.pipe()`链

**最佳实践：**
```python
# Good: Stays in expression API (parallelized)
df.with_columns(
    pl.col("value") * 10,
    pl.col("value").log(),
    pl.col("value").sqrt()
)

# Bad: Uses Python function (sequential)
df.with_columns(
    pl.col("value").map_elements(lambda x: x * 10)
)
```

## 严格类型系统

Polars 强制执行严格的输入：

**无静默转换：**
```python
# This will error - can't mix types
# df.with_columns(pl.col("int_col") + "string")

# Must cast explicitly
df.with_columns(
    pl.col("int_col").cast(pl.Utf8) + "_suffix"
)
```

**好处：**
- 防止无声错误
- 可预测的行为
- 更好的性能
- 更清晰的代码意图

**整数空值：**
与 pandas 不同，整数列可以有空值而无需转换为浮点数：
```python
# In pandas: Int column with null becomes Float
# In polars: Int column with null stays Int (with null values)
df = pl.DataFrame({"int_col": [1, 2, None, 4]})
# dtype: Int64 (not Float64)
```