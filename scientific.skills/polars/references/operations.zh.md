<!-- 此文件由机器翻译自 operations.md -->

# Polars 操作参考

本参考资料涵盖了所有常见的 Polars 操作以及全面的示例。

## 选择操作

### 选择列

**基本选择：**
```python
# Select specific columns
df.select("name", "age", "city")

# Using expressions
df.select(pl.col("name"), pl.col("age"))
```

**基于模式的选择：**
<<<代码块_1>>>

**计算列：**
<<<代码块_2>>>

### 带列（添加/修改）

添加新列或修改现有列，同时保留所有其他列：

<<<代码块_3>>>

## 过滤操作

### 基本过滤

<<<代码块_4>>>

### 高级过滤

**字符串操作：**
<<<代码块_5>>>

**会员检查：**
<<<代码块_6>>>

**范围过滤器：**
```python
# Between values
df.filter(pl.col("age").is_between(25, 35))

# Date range
df.filter(
    pl.col("date") >= pl.date(2023, 1, 1),
    pl.col("date") <= pl.date(2023, 12, 31)
)
```

**空过滤：**
```python
# Filter out nulls
df.filter(pl.col("value").is_not_null())

# Keep only nulls
df.filter(pl.col("value").is_null())
```

## 分组和聚合

### 基本分组依据

```python
# Group by single column
df.group_by("department").agg(
    pl.col("salary").mean().alias("avg_salary"),
    pl.len().alias("employee_count")
)

# Group by multiple columns
df.group_by("department", "location").agg(
    pl.col("salary").sum()
)

# Maintain order
df.group_by("category", maintain_order=True).agg(
    pl.col("value").sum()
)
```

### 聚合函数

**数量和长度：**
```python
df.group_by("category").agg(
    pl.len().alias("count"),
    pl.col("id").count().alias("non_null_count"),
    pl.col("id").n_unique().alias("unique_count")
)
```

**统计汇总：**
```python
df.group_by("group").agg(
    pl.col("value").sum().alias("total"),
    pl.col("value").mean().alias("average"),
    pl.col("value").median().alias("median"),
    pl.col("value").std().alias("std_dev"),
    pl.col("value").var().alias("variance"),
    pl.col("value").min().alias("minimum"),
    pl.col("value").max().alias("maximum"),
    pl.col("value").quantile(0.95).alias("p95")
)
```

**第一个和最后一个：**
```python
df.group_by("user_id").agg(
    pl.col("timestamp").first().alias("first_seen"),
    pl.col("timestamp").last().alias("last_seen"),
    pl.col("event").first().alias("first_event")
)
```

**列表聚合：**
```python
# Collect values into lists
df.group_by("category").agg(
    pl.col("item").alias("all_items")  # Creates list column
)
```

### 条件聚合

聚合内的过滤器：

```python
df.group_by("department").agg(
    # Count high earners
    (pl.col("salary") > 100000).sum().alias("high_earners"),

    # Average of filtered values
    pl.col("salary").filter(pl.col("bonus") > 0).mean().alias("avg_with_bonus"),

    # Conditional sum
    pl.when(pl.col("active"))
      .then(pl.col("sales"))
      .otherwise(0)
      .sum()
      .alias("active_sales")
)
```

### 多重聚合

有效地组合多个聚合：

```python
df.group_by("store_id").agg(
    pl.col("transaction_id").count().alias("num_transactions"),
    pl.col("amount").sum().alias("total_sales"),
    pl.col("amount").mean().alias("avg_transaction"),
    pl.col("customer_id").n_unique().alias("unique_customers"),
    pl.col("amount").max().alias("largest_transaction"),
    pl.col("timestamp").min().alias("first_transaction_date"),
    pl.col("timestamp").max().alias("last_transaction_date")
)
```

## 窗口函数

窗口函数应用聚合，同时保留原始行计数。

### 基本窗口操作

**团体统计：**
```python
# Add group mean to each row
df.with_columns(
    avg_age_by_dept=pl.col("age").mean().over("department")
)

# Multiple group columns
df.with_columns(
    group_avg=pl.col("value").mean().over("category", "region")
)
```

**排名：**
```python
df.with_columns(
    # Rank within groups
    rank=pl.col("score").rank().over("team"),

    # Dense rank (no gaps)
    dense_rank=pl.col("score").rank(method="dense").over("team"),

    # Row number
    row_num=pl.col("timestamp").sort().rank(method="ordinal").over("user_id")
)
```

### 窗口映射策略

**group_to_rows（默认）：**
保留原始行顺序：
```python
df.with_columns(
    group_mean=pl.col("value").mean().over("category", mapping_strategy="group_to_rows")
)
```

**爆炸：**
更快地将行分组在一起：
```python
df.with_columns(
    group_mean=pl.col("value").mean().over("category", mapping_strategy="explode")
)
```

**加入：**
创建列表列：
```python
df.with_columns(
    group_values=pl.col("value").over("category", mapping_strategy="join")
)
```

### 滚动窗口

**基于时间的滚动：**
```python
df.with_columns(
    rolling_avg=pl.col("value").rolling_mean(
        window_size="7d",
        by="date"
    )
)
```

**基于行的滚动：**
```python
df.with_columns(
    rolling_sum=pl.col("value").rolling_sum(window_size=3),
    rolling_max=pl.col("value").rolling_max(window_size=5)
)
```

### 累积运算

```python
df.with_columns(
    cumsum=pl.col("value").cum_sum().over("group"),
    cummax=pl.col("value").cum_max().over("group"),
    cummin=pl.col("value").cum_min().over("group"),
    cumprod=pl.col("value").cum_prod().over("group")
)
```

### 平移和滞后/超前

```python
df.with_columns(
    # Previous value (lag)
    prev_value=pl.col("value").shift(1).over("user_id"),

    # Next value (lead)
    next_value=pl.col("value").shift(-1).over("user_id"),

    # Calculate difference from previous
    diff=pl.col("value") - pl.col("value").shift(1).over("user_id")
)
```

## 排序

### 基本排序

```python
# Sort by single column
df.sort("age")

# Sort descending
df.sort("age", descending=True)

# Sort by multiple columns
df.sort("department", "age")

# Mixed sorting order
df.sort(["department", "salary"], descending=[False, True])
```

### 高级排序

**空处理：**
```python
# Nulls first
df.sort("value", nulls_last=False)

# Nulls last
df.sort("value", nulls_last=True)
```

**按表达式排序：**
```python
# Sort by computed value
df.sort(pl.col("first_name").str.len())

# Sort by multiple expressions
df.sort(
    pl.col("last_name").str.to_lowercase(),
    pl.col("age").abs()
)
```

## 条件操作

### 当/然后/否则

```python
# Basic conditional
df.with_columns(
    status=pl.when(pl.col("age") >= 18)
        .then("adult")
        .otherwise("minor")
)

# Multiple conditions
df.with_columns(
    category=pl.when(pl.col("score") >= 90)
        .then("A")
        .when(pl.col("score") >= 80)
        .then("B")
        .when(pl.col("score") >= 70)
        .then("C")
        .otherwise("F")
)

# Conditional computation
df.with_columns(
    adjusted_price=pl.when(pl.col("is_member"))
        .then(pl.col("price") * 0.9)
        .otherwise(pl.col("price"))
)
```

## 字符串操作

### 常用字符串方法

```python
df.with_columns(
    # Case conversion
    upper=pl.col("name").str.to_uppercase(),
    lower=pl.col("name").str.to_lowercase(),
    title=pl.col("name").str.to_titlecase(),

    # Trimming
    trimmed=pl.col("text").str.strip_chars(),

    # Substring
    first_3=pl.col("name").str.slice(0, 3),

    # Replace
    cleaned=pl.col("text").str.replace("old", "new"),
    cleaned_all=pl.col("text").str.replace_all("old", "new"),

    # Split
    parts=pl.col("full_name").str.split(" "),

    # Length
    name_length=pl.col("name").str.len_chars()
)
```

### 字符串过滤

```python
# Contains
df.filter(pl.col("email").str.contains("@gmail.com"))

# Starts/ends with
df.filter(pl.col("name").str.starts_with("A"))
df.filter(pl.col("file").str.ends_with(".csv"))

# Regex matching
df.filter(pl.col("phone").str.contains(r"^\d{3}-\d{4}$"))
```

## 日期和时间操作

### 日期解析

```python
# Parse strings to dates
df.with_columns(
    date=pl.col("date_str").str.strptime(pl.Date, "%Y-%m-%d"),
    datetime=pl.col("dt_str").str.strptime(pl.Datetime, "%Y-%m-%d %H:%M:%S")
)
```

### 日期组件

```python
df.with_columns(
    year=pl.col("date").dt.year(),
    month=pl.col("date").dt.month(),
    day=pl.col("date").dt.day(),
    weekday=pl.col("date").dt.weekday(),
    hour=pl.col("datetime").dt.hour(),
    minute=pl.col("datetime").dt.minute()
)
```

### 日期算术

```python
# Add duration
df.with_columns(
    next_week=pl.col("date") + pl.duration(weeks=1),
    next_month=pl.col("date") + pl.duration(months=1)
)

# Difference between dates
df.with_columns(
    days_diff=(pl.col("end_date") - pl.col("start_date")).dt.total_days()
)
```

### 日期过滤

```python
# Filter by date range
df.filter(
    pl.col("date").is_between(pl.date(2023, 1, 1), pl.date(2023, 12, 31))
)

# Filter by year
df.filter(pl.col("date").dt.year() == 2023)

# Filter by month
df.filter(pl.col("date").dt.month().is_in([6, 7, 8]))  # Summer months
```

## 列表操作

### 使用列表列

```python
# Create list column
df.with_columns(
    items_list=pl.col("item1", "item2", "item3").to_list()
)

# List operations
df.with_columns(
    list_len=pl.col("items").list.len(),
    first_item=pl.col("items").list.first(),
    last_item=pl.col("items").list.last(),
    unique_items=pl.col("items").list.unique(),
    sorted_items=pl.col("items").list.sort()
)

# Explode lists to rows
df.explode("items")

# Filter list elements
df.with_columns(
    filtered=pl.col("items").list.eval(pl.element() > 10)
)
```

## 结构体操作

### 使用嵌套结构

```python
# Create struct column
df.with_columns(
    address=pl.struct(["street", "city", "zip"])
)

# Access struct fields
df.with_columns(
    city=pl.col("address").struct.field("city")
)

# Unnest struct to columns
df.unnest("address")
```

## 唯一和重复操作

```python
# Get unique rows
df.unique()

# Unique on specific columns
df.unique(subset=["name", "email"])

# Keep first/last duplicate
df.unique(subset=["id"], keep="first")
df.unique(subset=["id"], keep="last")

# Identify duplicates
df.with_columns(
    is_duplicate=pl.col("id").is_duplicated()
)

# Count duplicates
df.group_by("email").agg(
    pl.len().alias("count")
).filter(pl.col("count") > 1)
```

## 采样

```python
# Random sample
df.sample(n=100)

# Sample fraction
df.sample(fraction=0.1)

# Sample with seed for reproducibility
df.sample(n=100, seed=42)
```

## 列重命名

```python
# Rename specific columns
df.rename({"old_name": "new_name", "age": "years"})

# Rename with expression
df.select(pl.col("*").name.suffix("_renamed"))
df.select(pl.col("*").name.prefix("data_"))
df.select(pl.col("*").name.to_uppercase())
```