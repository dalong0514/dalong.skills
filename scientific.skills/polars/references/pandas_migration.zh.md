<!-- 此文件由机器翻译自 pandas_migration.md -->

# 大熊猫到极地迁徙指南

本指南通过全面的操作映射和关键差异帮助您从 pandas 迁移到 Polars。

## 核心概念差异

### 1.没有索引系统

**Pandas：** 使用基于行的索引系统
```python
df.loc[0, "column"]
df.iloc[0:5]
df.set_index("id")
```

**极坐标：** 仅使用整数位置
<<<代码块_1>>>

### 2. 内存格式

**Pandas：** 面向行的 NumPy 数组
**极坐标：** 柱状阿帕奇箭头格式

**影响：**
- Polars 的色谱柱操作速度更快
- Polars 使用更少的内存
- Polars具有更好的数据共享能力

### 3. 并行化

**Pandas：** 主要是单线程（需要 Dask 来实现并行性）
**Polars：** 默认情况下使用 Rust 的并发性并行

### 4. 惰性评估

**Pandas：** 只有急切的评估
**Polars：** 具有查询优化的急切 (DataFrame) 和惰性 (LazyFrame)

### 5.类型严格

**Pandas：** 允许静默类型转换
**Polars：** 严格输入，需要显式转换

**示例：**
<<<代码块_2>>>

## 操作映射

### 数据选择

|运营|熊猫 |极地 |
|------------|--------|--------|
|选择列 | `df["col"]` 或 `df.col` | `df.select("col")` 或 `df["col"]` |
|选择多个 | `df[["a", "b"]]` | `df.select("a", "b")` |
|按位置选择 | `df.iloc[:, 0:3]` | `df.select(pl.col(df.columns[0:3]))` |
|按条件选择| `df[df["age"] > 25]` | `df.filter(pl.col("age") > 25)` |

### 数据过滤

|运营|熊猫 |极地 |
|------------|--------|--------|
|单一条件 | `df[df["age"] > 25]` | `df.filter(pl.col("age") > 25)` |
|多种条件 | `df[(df["age"] > 25) & (df["city"] == "NY")]` | `df.filter(pl.col("age") > 25, pl.col("city") == "NY")` |
|查询方式 | `df.query("age > 25")` | `df.filter(pl.col("age") > 25)` |
|伊辛 | `df[df["city"].isin(["NY", "LA"])]` | `df.filter(pl.col("city").is_in(["NY", "LA"]))` |
|伊斯纳 | `df[df["value"].isna()]` | `df.filter(pl.col("value").is_null())` |
|不娜| `df[df["value"].notna()]` | `df.filter(pl.col("value").is_not_null())` |

### 添加/修改列

|运营|熊猫 |极地 |
|------------|--------|--------|
|添加栏目| `df["new"] = df["old"] * 2` | `df.with_columns(new=pl.col("old") * 2)` |
|多栏 | `df.assign(a=..., b=...)` | `df.with_columns(a=..., b=...)` |
|条件栏 | `np.where(condition, a, b)` | `pl.when(condition).then(a).otherwise(b)` |

**重要区别 - 并行执行：**

<<<代码块_3>>>

### 分组和聚合

|运营|熊猫 |极地 |
|------------|--------|--------|
|分组依据 | `df.groupby("col")` | `df.group_by("col")` |
|聚合单| `df.groupby("col")["val"].mean()` | `df.group_by("col").agg(pl.col("val").mean())` |
|聚合多个 | `df.groupby("col").agg({"val": ["mean", "sum"]})` | `df.group_by("col").agg(pl.col("val").mean(), pl.col("val").sum())` |
|尺寸| `df.groupby("col").size()` | `df.group_by("col").agg(pl.len())` |
|计数 | `df.groupby("col").count()` | `df.group_by("col").agg(pl.col("*").count())` |

### 窗口函数

|运营|熊猫 |极地 |
|------------|--------|--------|
|转变| `df.groupby("col").transform("mean")` | `df.with_columns(pl.col("val").mean().over("col"))` |
|排名| `df.groupby("col")["val"].rank()` | `df.with_columns(pl.col("val").rank().over("col"))` |
|班次| `df.groupby("col")["val"].shift(1)` | `df.with_columns(pl.col("val").shift(1).over("col"))` |
|累积 | `df.groupby("col")["val"].cumsum()` | `df.with_columns(pl.col("val").cum_sum().over("col"))` |

### 加入

|运营|熊猫 |极地 |
|------------|--------|--------|
|内连接| `df1.merge(df2, on="id")` | `df1.join(df2, on="id", how="inner")` |
|左连接| `df1.merge(df2, on="id", how="left")` | `df1.join(df2, on="id", how="left")` |
|不同的键| `df1.merge(df2, left_on="a", right_on="b")` | `df1.join(df2, left_on="a", right_on="b")` |

### 连接

|运营|熊猫 |极地 |
|------------|--------|--------|
|垂直| `pd.concat([df1, df2], axis=0)` | `pl.concat([df1, df2], how="vertical")` |
|卧式| `pd.concat([df1, df2], axis=1)` | `pl.concat([df1, df2], how="horizontal")` |

### 排序

|运营|熊猫 |极地 |
|------------|--------|--------|
|按列排序| `df.sort_values("col")` | `df.sort("col")` |
|降序| `df.sort_values("col", ascending=False)` | `df.sort("col", descending=True)` |
|多栏| `df.sort_values(["a", "b"])` | `df.sort("a", "b")` |

### 重塑

|运营|熊猫 |极地 |
|------------|--------|--------|
|枢轴| `df.pivot(index="a", columns="b", values="c")` | `df.pivot(values="c", index="a", columns="b")` |
|融化 | `df.melt(id_vars="id")` | `df.unpivot(index="id")` |

### I/O 操作

|运营|熊猫 |极地 |
|------------|--------|--------|
|读取 CSV | `pd.read_csv("file.csv")` | `pl.read_csv("file.csv")` 或 `pl.scan_csv()` |
|写入 CSV | `df.to_csv("file.csv")` | `df.write_csv("file.csv")` |
|阅读镶木地板 | `pd.read_parquet("file.parquet")` | `pl.read_parquet("file.parquet")` |
|写实木复合地板| `df.to_parquet("file.parquet")` | `df.write_parquet("file.parquet")` |
|阅读 Excel | `pd.read_excel("file.xlsx")` | `pl.read_excel("file.xlsx")` |

### 字符串操作

|运营|熊猫 |极地 |
|------------|--------|--------|
|上层| `df["col"].str.upper()` | `df.select(pl.col("col").str.to_uppercase())` |
|降低| `df["col"].str.lower()` | `df.select(pl.col("col").str.to_lowercase())` |
|包含 | `df["col"].str.contains("pattern")` | `df.filter(pl.col("col").str.contains("pattern"))` |
|更换| `df["col"].str.replace("old", "new")` | `df.select(pl.col("col").str.replace("old", "new"))` |
|分裂| `df["col"].str.split(" ")` | `df.select(pl.col("col").str.split(" "))` |

### 日期时间操作

|运营|熊猫 |极地 |
|------------|--------|--------|
|解析日期 | `pd.to_datetime(df["col"])` | `df.select(pl.col("col").str.strptime(pl.Date, "%Y-%m-%d"))` |
|年份| `df["date"].dt.year` | `df.select(pl.col("date").dt.year())` |
|月 | `df["date"].dt.month` | `df.select(pl.col("date").dt.month())` |
|日 | `df["date"].dt.day` | `df.select(pl.col("date").dt.day())` |

### 缺失数据

|运营|熊猫 |极地 |
|------------|--------|--------|
|删除空值 | `df.dropna()` | `df.drop_nulls()` |
|填充空值 | `df.fillna(0)` | `df.fill_null(0)` |
|检查空 | `df["col"].isna()` | `df.select(pl.col("col").is_null())` |
|前向填充 | `df.fillna(method="ffill")` | `df.select(pl.col("col").fill_null(strategy="forward"))` |

### 其他操作

|运营|熊猫 |极地 |
|------------|--------|--------|
|独特的价值观| `df["col"].unique()` | `df["col"].unique()` |
|价值很重要| `df["col"].value_counts()` | `df["col"].value_counts()` |
|描述| `df.describe()` | `df.describe()` |
|样品| `df.sample(n=100)` | `df.sample(n=100)` |
|头| `df.head()` | `df.head()` |
|尾巴| `df.tail()` | `df.tail()` |

## 常见迁移模式

### 模式 1：链式操作

**熊猫：**
<<<代码块_4>>>

**极地：**
<<<代码块_5>>>

### 模式 2：应用函数

**熊猫：**
<<<代码块_6>>>

**极地：**
```python
# Use expressions instead
df = df.with_columns(result=pl.col("value") * 2)

# If custom function needed
df = df.with_columns(
    result=pl.col("value").map_elements(lambda x: x * 2, return_dtype=pl.Float64)
)
```

### 模式 3：条件列创建

**熊猫：**
```python
df["category"] = np.where(
    df["value"] > 100,
    "high",
    np.where(df["value"] > 50, "medium", "low")
)
```

**极地：**
```python
df = df.with_columns(
    category=pl.when(pl.col("value") > 100)
        .then("high")
        .when(pl.col("value") > 50)
        .then("medium")
        .otherwise("low")
)
```

### 模式 4：组变换

**熊猫：**
```python
df["group_mean"] = df.groupby("category")["value"].transform("mean")
```

**极地：**
```python
df = df.with_columns(
    group_mean=pl.col("value").mean().over("category")
)
```

### 模式 5：多重聚合

**熊猫：**
```python
result = df.groupby("category").agg({
    "value": ["mean", "sum", "count"],
    "price": ["min", "max"]
})
```

**极地：**
```python
result = df.group_by("category").agg(
    pl.col("value").mean().alias("value_mean"),
    pl.col("value").sum().alias("value_sum"),
    pl.col("value").count().alias("value_count"),
    pl.col("price").min().alias("price_min"),
    pl.col("price").max().alias("price_max")
)
```

## 要避免的性能反模式

### 反模式 1：顺序管道操作

**不好（禁用并行化）：**
```python
df = df.pipe(function1).pipe(function2).pipe(function3)
```

**好（启用并行化）：**
```python
df = df.with_columns(
    function1_result(),
    function2_result(),
    function3_result()
)
```

### 反模式 2：热路径中的 Python 函数

**不好：**
```python
df = df.with_columns(
    result=pl.col("value").map_elements(lambda x: x * 2)
)
```

**好：**
```python
df = df.with_columns(result=pl.col("value") * 2)
```

### 反模式 3：使用急切读取大文件

**不好：**
```python
df = pl.read_csv("large_file.csv")
result = df.filter(pl.col("age") > 25).select("name", "age")
```

**好：**
```python
lf = pl.scan_csv("large_file.csv")
result = lf.filter(pl.col("age") > 25).select("name", "age").collect()
```

### 反模式 4：行迭代

**不好：**
```python
for row in df.iter_rows():
    # Process row
    pass
```

**好：**
```python
# Use vectorized operations
df = df.with_columns(
    # Vectorized computation
)
```

## 迁移清单

从 pandas 迁移到 Polars 时：

1. **删除索引操作** - 使用整数位置或group_by
2. **用表达式替换apply/map** - 使用Polars原生操作
3. **更新列分配** - 使用 `with_columns()` 而不是直接分配
4. **将 groupby.transform 更改为 .over()** - 窗口函数的工作方式不同
5. **更新字符串操作** - 使用 `.str.to_uppercase()` 而不是 `.str.upper()`
6. **添加显式类型转换** - Polars 不会默默地转换类型
7. **考虑惰性求值** - 对于大数据使用 `scan_*` 而不是 `read_*`
8. **更新聚合语法** - 在 Polars 中更明确
9. **删除reset_index调用** - 在Polars中不需要
10. **更新条件逻辑** - 使用 `when().then().otherwise()` 模式

## 兼容层

对于逐步迁移，您可以使用这两个库：

```python
import pandas as pd
import polars as pl

# Convert pandas to Polars
pl_df = pl.from_pandas(pd_df)

# Convert Polars to pandas
pd_df = pl_df.to_pandas()

# Use Arrow for zero-copy (when possible)
pl_df = pl.from_arrow(pd_df)
pd_df = pl_df.to_arrow().to_pandas()
```

## 何时坚持使用 Pandas

在以下情况下考虑与熊猫待在一起：
- 处理需要复杂索引操作的时间序列
- 需要广泛的生态系统支持（一些库仅支持 pandas）
- 团队缺乏 Rust/Polars 专业知识
- 数据很小，性能并不重要
- 使用高级 pandas 功能，无需 Polars 等效功能

## 何时切换到 Polar

在以下情况下切换到 Polar：
- 性能至关重要
- 处理大型数据集（>1GB）
- 需要惰性评估和查询优化
- 想要更好的类型安全
- 默认需要并行执行
- 开始一个新项目