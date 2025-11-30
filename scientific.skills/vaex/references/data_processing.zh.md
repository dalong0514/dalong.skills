<!-- 此文件由机器翻译自 data_processing.md -->

# 数据处理和操作

本参考内容涵盖 Vaex 中的过滤、选择、虚拟列、表达式、聚合、groupby 操作和数据转换。

## 过滤和选择

Vaex 使用布尔表达式来高效过滤数据，无需复制：

### 基本过滤

```python
# Simple filter
df_filtered = df[df.age > 25]

# Multiple conditions
df_filtered = df[(df.age > 25) & (df.salary > 50000)]
df_filtered = df[(df.category == 'A') | (df.category == 'B')]

# Negation
df_filtered = df[~(df.age < 18)]
```

### 选择对象

Vaex 可以同时维护多个命名选择：

<<<代码块_1>>>

### 高级过滤

<<<代码块_2>>>

## 虚拟列和表达式

虚拟列是即时计算的，内存开销为零：

### 创建虚拟列

<<<代码块_3>>>

### 表达方法

<<<代码块_4>>>

### 条件表达式

<<<代码块_5>>>

## 字符串操作

通过 `.str` 访问器访问字符串方法：

### 基本字符串方法

<<<代码块_6>>>

### 高级字符串操作

```python
# Replacing
df['clean_text'] = df.text.str.replace('bad', 'good')

# Splitting (returns first part)
df['first_name'] = df.full_name.str.split(' ')[0]

# Concatenation
df['full_name'] = df.first_name + ' ' + df.last_name

# Padding
df['padded'] = df.code.str.pad(10, '0', 'left')  # Zero-padding
```

## 日期时间操作

通过 `.dt` 访问器访问日期时间方法：

### 日期时间属性

```python
# Parsing strings to datetime
df['date_parsed'] = df.date_string.astype('datetime64')

# Extracting components
df['year'] = df.timestamp.dt.year
df['month'] = df.timestamp.dt.month
df['day'] = df.timestamp.dt.day
df['hour'] = df.timestamp.dt.hour
df['minute'] = df.timestamp.dt.minute
df['second'] = df.timestamp.dt.second

# Day of week
df['weekday'] = df.timestamp.dt.dayofweek  # 0=Monday
df['day_name'] = df.timestamp.dt.day_name  # 'Monday', 'Tuesday', ...

# Date arithmetic
df['tomorrow'] = df.date + pd.Timedelta(days=1)
df['next_week'] = df.date + pd.Timedelta(weeks=1)
```

## 聚合

Vaex 在数十亿行中高效执行聚合：

### 基本聚合

```python
# Single column
mean_age = df.age.mean()
std_age = df.age.std()
min_age = df.age.min()
max_age = df.age.max()
sum_sales = df.sales.sum()
count_rows = df.count()

# With selections
mean_adult_age = df.age.mean(selection='adults')

# Multiple at once with delay
mean = df.age.mean(delay=True)
std = df.age.std(delay=True)
results = vaex.execute([mean, std])
```

### 可用的聚合函数

```python
# Central tendency
df.x.mean()
df.x.median_approx()  # Approximate median (fast)

# Dispersion
df.x.std()           # Standard deviation
df.x.var()           # Variance
df.x.min()
df.x.max()
df.x.minmax()        # Both min and max

# Count
df.count()           # Total rows
df.x.count()         # Non-missing values

# Sum and product
df.x.sum()
df.x.prod()

# Percentiles
df.x.quantile(0.5)           # Median
df.x.quantile([0.25, 0.75])  # Quartiles

# Correlation
df.correlation(df.x, df.y)
df.covar(df.x, df.y)

# Higher moments
df.x.kurtosis()
df.x.skew()

# Unique values
df.x.nunique()       # Count unique
df.x.unique()        # Get unique values (returns array)
```

## 分组操作

对数据进行分组并计算每组的聚合：

### 基本分组依据

```python
# Single column groupby
grouped = df.groupby('category')

# Aggregation
result = grouped.agg({'sales': 'sum'})
result = grouped.agg({'sales': 'sum', 'quantity': 'mean'})

# Multiple aggregations on same column
result = grouped.agg({
    'sales': ['sum', 'mean', 'std'],
    'quantity': 'sum'
})
```

### 高级分组依据

```python
# Multiple grouping columns
result = df.groupby(['category', 'region']).agg({
    'sales': 'sum',
    'quantity': 'mean'
})

# Custom aggregation functions
result = df.groupby('category').agg({
    'sales': lambda x: x.max() - x.min()
})

# Available aggregation functions
# 'sum', 'mean', 'std', 'min', 'max', 'count', 'first', 'last'
```

### GroupBy 与 Binning

```python
# Bin continuous variable and aggregate
result = df.groupby(vaex.vrange(0, 100, 10)).agg({
    'sales': 'sum'
})

# Datetime binning
result = df.groupby(df.timestamp.dt.year).agg({
    'sales': 'sum'
})
```

## 分箱和离散化

从连续变量创建箱：

### 简单分箱

```python
# Create bins
df['age_bin'] = df.age.digitize([18, 30, 50, 65, 100])

# Labeled bins
bins = [0, 18, 30, 50, 65, 100]
labels = ['child', 'young_adult', 'adult', 'middle_age', 'senior']
df['age_group'] = df.age.digitize(bins)
# Note: Apply labels using where() or mapping
```

### 统计分箱

```python
# Equal-width bins
df['value_bin'] = df.value.digitize(
    vaex.vrange(df.value.min(), df.value.max(), 10)
)

# Quantile-based bins
quantiles = df.value.quantile([0.25, 0.5, 0.75])
df['value_quartile'] = df.value.digitize(quantiles)
```

## 多维聚合

计算网格上的统计数据：

```python
# 2D histogram/heatmap data
counts = df.count(binby=[df.x, df.y], limits=[[0, 10], [0, 10]], shape=(100, 100))

# Mean on a grid
mean_z = df.mean(df.z, binby=[df.x, df.y], limits=[[0, 10], [0, 10]], shape=(50, 50))

# Multiple statistics on grid
stats = df.mean(df.z, binby=[df.x, df.y], shape=(50, 50), delay=True)
counts = df.count(binby=[df.x, df.y], shape=(50, 50), delay=True)
results = vaex.execute([stats, counts])
```

## 处理缺失数据

处理缺失值、空值和 NaN 值：

### 检测缺失数据

```python
# Check for missing
df['age_missing'] = df.age.isna()
df['age_present'] = df.age.notna()

# Count missing
missing_count = df.age.isna().sum()
missing_pct = df.age.isna().mean() * 100
```

### 处理缺失数据

```python
# Filter out missing
df_clean = df[df.age.notna()]

# Fill missing with value
df['age_filled'] = df.age.fillna(0)
df['age_filled'] = df.age.fillna(df.age.mean())

# Forward/backward fill (for time series)
df['age_ffill'] = df.age.fillna(method='ffill')
df['age_bfill'] = df.age.fillna(method='bfill')
```

### Vaex 中缺少数据类型

Vaex 区分：
- **NaN** - IEEE 浮点非数字
- **NA** - 箭头空类型
- **缺失** - 缺失数据的通用术语

```python
# Check which missing type
df.is_masked('column_name')  # True if uses Arrow null (NA)

# Convert between types
df['col_masked'] = df.col.as_masked()  # Convert to NA representation
```

## 排序

```python
# Sort by single column
df_sorted = df.sort('age')
df_sorted = df.sort('age', ascending=False)

# Sort by multiple columns
df_sorted = df.sort(['category', 'age'])

# Note: Sorting materializes a new column with indices
# For very large datasets, consider if sorting is necessary
```

## 连接数据框

根据键组合 DataFrame：

```python
# Inner join
df_joined = df1.join(df2, on='key_column')

# Left join
df_joined = df1.join(df2, on='key_column', how='left')

# Join on different column names
df_joined = df1.join(
    df2,
    left_on='id',
    right_on='user_id',
    how='left'
)

# Multiple key columns
df_joined = df1.join(df2, on=['key1', 'key2'])
```

## 添加和删除列

### 添加列

```python
# Virtual column (no memory)
df['new_col'] = df.x + df.y

# From external array (must match length)
import numpy as np
new_data = np.random.rand(len(df))
df['random'] = new_data

# Constant value
df['constant'] = 42
```

### 删除列

```python
# Drop single column
df = df.drop('column_name')

# Drop multiple columns
df = df.drop(['col1', 'col2', 'col3'])

# Select specific columns (drop others)
df = df[['col1', 'col2', 'col3']]
```

### 重命名列

```python
# Rename single column
df = df.rename('old_name', 'new_name')

# Rename multiple columns
df = df.rename({
    'old_name1': 'new_name1',
    'old_name2': 'new_name2'
})
```

## 常见模式

### 模式：复杂特征工程

```python
# Multiple derived features
df['log_price'] = df.price.log()
df['price_per_unit'] = df.price / df.quantity
df['is_discount'] = df.discount > 0
df['price_category'] = (df.price > 100).where('expensive', 'affordable')
df['revenue'] = df.price * df.quantity * (1 - df.discount)
```

### 模式：文本清理

```python
# Clean and standardize text
df['email_clean'] = df.email.str.lower().str.strip()
df['has_valid_email'] = df.email_clean.str.contains('@')
df['domain'] = df.email_clean.str.split('@')[1]
```

### 模式：基于时间的分析

```python
# Extract temporal features
df['year'] = df.timestamp.dt.year
df['month'] = df.timestamp.dt.month
df['day_of_week'] = df.timestamp.dt.dayofweek
df['is_weekend'] = df.day_of_week >= 5
df['quarter'] = ((df.month - 1) // 3) + 1
```

### 模式：分组统计

```python
# Compute statistics by group
monthly_sales = df.groupby(df.timestamp.dt.month).agg({
    'revenue': ['sum', 'mean', 'count'],
    'quantity': 'sum'
})

# Multiple grouping levels
category_region_sales = df.groupby(['category', 'region']).agg({
    'sales': 'sum',
    'profit': 'mean'
})
```

## 性能提示

1. **使用虚拟列** - 它们是即时计算的，没有内存成本
2. **延迟=True的批量操作** - 一次计算多个聚合
3. **避免 `.values` 或 `.to_pandas_df()`** - 尽可能保持惰性操作
4. **使用选择** - 多个命名选择比创建新的 DataFrame 更有效
5. **利用表达式** - 它们支持查询优化
6. **最小化排序** - 在大型数据集上排序的成本很高

## 相关资源

- 对于 DataFrame 创建：请参阅 `core_dataframes.md`
- 性能优化：参见`performance.md`
- 可视化：参见`visualization.md`
- 对于 ML 管道：请参阅 `machine_learning.md`