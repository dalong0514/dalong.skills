<!-- 此文件由机器翻译自 io_guide.md -->

# Polars 数据 I/O 指南

使用 Polars 读取和写入各种格式数据的综合指南。

## CSV 文件

### 读取 CSV

**Eager 模式（加载到内存中）：**
```python
import polars as pl

# Basic read
df = pl.read_csv("data.csv")

# With options
df = pl.read_csv(
    "data.csv",
    separator=",",
    has_header=True,
    columns=["col1", "col2"],  # Select specific columns
    n_rows=1000,  # Read only first 1000 rows
    skip_rows=10,  # Skip first 10 rows
    dtypes={"col1": pl.Int64, "col2": pl.Utf8},  # Specify types
    null_values=["NA", "null", ""],  # Define null values
    encoding="utf-8",
    ignore_errors=False
)
```

**惰性模式（扫描而不加载 - 建议用于大文件）：**
<<<代码块_1>>>

### 写入 CSV

<<<代码块_2>>>

### 多个 CSV 文件

**读取多个文件：**
<<<代码块_3>>>

## 镶木地板文件

Parquet 是出于性能和压缩的考虑而推荐的格式。

### 阅读镶木地板

**渴望：**
<<<代码块_4>>>

**懒惰（推荐）：**
<<<代码块_5>>>

### 写镶木地板

<<<代码块_6>>>

### 分区镶木地板（蜂巢式）

**写入分区：**
```python
# Write with partitioning
df.write_parquet(
    "output_dir",
    partition_by=["year", "month"]  # Creates directory structure
)
# Creates: output_dir/year=2023/month=01/data.parquet
```

**读取分区：**
```python
lf = pl.scan_parquet("output_dir/**/*.parquet")

# Hive partitioning columns are automatically added
result = lf.filter(pl.col("year") == 2023).collect()
```

## JSON 文件

### 读取 JSON

**NDJSON（换行符分隔的 JSON）- 推荐：**
```python
df = pl.read_ndjson("data.ndjson")

# Lazy
lf = pl.scan_ndjson("data.ndjson")
```

**标准 JSON：**
```python
df = pl.read_json("data.json")

# From JSON string
df = pl.read_json('{"col1": [1, 2], "col2": ["a", "b"]}')
```

### 编写 JSON

```python
# Write NDJSON
df.write_ndjson("output.ndjson")

# Write standard JSON
df.write_json("output.json")

# Pretty printed
df.write_json("output.json", pretty=True, row_oriented=False)
```

## Excel 文件

### 读取 Excel

```python
# Read first sheet
df = pl.read_excel("data.xlsx")

# Specific sheet
df = pl.read_excel("data.xlsx", sheet_name="Sheet1")
# Or by index
df = pl.read_excel("data.xlsx", sheet_id=0)

# With options
df = pl.read_excel(
    "data.xlsx",
    sheet_name="Sheet1",
    columns=["A", "B", "C"],  # Excel columns
    n_rows=100,
    skip_rows=5,
    has_header=True
)
```

### 编写 Excel

```python
# Write to Excel
df.write_excel("output.xlsx")

# Multiple sheets
with pl.ExcelWriter("output.xlsx") as writer:
    df1.write_excel(writer, worksheet="Sheet1")
    df2.write_excel(writer, worksheet="Sheet2")
```

## 数据库连接

### 从数据库读取

```python
import polars as pl

# Read entire table
df = pl.read_database("SELECT * FROM users", connection_uri="postgresql://...")

# Using connectorx for better performance
df = pl.read_database_uri(
    "SELECT * FROM users WHERE age > 25",
    uri="postgresql://user:pass@localhost/db"
)
```

### 写入数据库

```python
# Using SQLAlchemy
from sqlalchemy import create_engine

engine = create_engine("postgresql://user:pass@localhost/db")
df.write_database("table_name", connection=engine)

# With options
df.write_database(
    "table_name",
    connection=engine,
    if_exists="replace",  # or "append", "fail"
)
```

### 通用数据库连接器

**PostgreSQL：**
```python
uri = "postgresql://username:password@localhost:5432/database"
df = pl.read_database_uri("SELECT * FROM table", uri=uri)
```

**MySQL：**
```python
uri = "mysql://username:password@localhost:3306/database"
df = pl.read_database_uri("SELECT * FROM table", uri=uri)
```

**SQLite：**
```python
uri = "sqlite:///path/to/database.db"
df = pl.read_database_uri("SELECT * FROM table", uri=uri)
```

## 云存储

### AWS S3

```python
# Read from S3
df = pl.read_parquet("s3://bucket/path/to/file.parquet")
lf = pl.scan_parquet("s3://bucket/path/*.parquet")

# Write to S3
df.write_parquet("s3://bucket/path/output.parquet")

# With credentials
import os
os.environ["AWS_ACCESS_KEY_ID"] = "your_key"
os.environ["AWS_SECRET_ACCESS_KEY"] = "your_secret"
os.environ["AWS_REGION"] = "us-west-2"

df = pl.read_parquet("s3://bucket/file.parquet")
```

### Azure Blob 存储

```python
# Read from Azure
df = pl.read_parquet("az://container/path/file.parquet")

# Write to Azure
df.write_parquet("az://container/path/output.parquet")

# With credentials
os.environ["AZURE_STORAGE_ACCOUNT_NAME"] = "account"
os.environ["AZURE_STORAGE_ACCOUNT_KEY"] = "key"
```

### 谷歌云存储（GCS）

```python
# Read from GCS
df = pl.read_parquet("gs://bucket/path/file.parquet")

# Write to GCS
df.write_parquet("gs://bucket/path/output.parquet")

# With credentials
os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = "/path/to/credentials.json"
```

## 谷歌 BigQuery

```python
# Read from BigQuery
df = pl.read_database(
    "SELECT * FROM project.dataset.table",
    connection_uri="bigquery://project"
)

# Or using Google Cloud SDK
from google.cloud import bigquery
client = bigquery.Client()

query = "SELECT * FROM project.dataset.table WHERE date > '2023-01-01'"
df = pl.from_pandas(client.query(query).to_dataframe())
```

##阿帕奇箭

### IPC/羽毛格式

**阅读：**
```python
df = pl.read_ipc("data.arrow")
lf = pl.scan_ipc("data.arrow")
```

**写：**
```python
df.write_ipc("output.arrow")

# Compressed
df.write_ipc("output.arrow", compression="zstd")
```

### 箭头流

```python
# Write streaming format
df.write_ipc("output.arrows", compression="zstd")

# Read streaming
df = pl.read_ipc("output.arrows")
```

### 从/到箭头

```python
import pyarrow as pa

# From Arrow Table
arrow_table = pa.table({"col": [1, 2, 3]})
df = pl.from_arrow(arrow_table)

# To Arrow Table
arrow_table = df.to_arrow()
```

## 内存中格式

### Python 字典

```python
# From dict
df = pl.DataFrame({
    "col1": [1, 2, 3],
    "col2": ["a", "b", "c"]
})

# To dict
data_dict = df.to_dict()  # Column-oriented
data_dict = df.to_dict(as_series=False)  # Lists instead of Series
```

### NumPy 数组

```python
import numpy as np

# From NumPy
arr = np.array([[1, 2], [3, 4], [5, 6]])
df = pl.DataFrame(arr, schema=["col1", "col2"])

# To NumPy
arr = df.to_numpy()
```

### 熊猫数据框

```python
import pandas as pd

# From Pandas
pd_df = pd.DataFrame({"col": [1, 2, 3]})
pl_df = pl.from_pandas(pd_df)

# To Pandas
pd_df = pl_df.to_pandas()

# Zero-copy when possible
pl_df = pl.from_arrow(pd_df)
```

### 行列表

```python
# From list of dicts
data = [
    {"name": "Alice", "age": 25},
    {"name": "Bob", "age": 30}
]
df = pl.DataFrame(data)

# To list of dicts
rows = df.to_dicts()

# From list of tuples
data = [("Alice", 25), ("Bob", 30)]
df = pl.DataFrame(data, schema=["name", "age"])
```

## 流式传输大文件

对于大于内存的数据集，请使用带有流式传输的惰性模式：

```python
# Streaming mode
lf = pl.scan_csv("very_large.csv")
result = lf.filter(pl.col("value") > 100).collect(streaming=True)

# Streaming with multiple files
lf = pl.scan_parquet("data/*.parquet")
result = lf.group_by("category").agg(pl.col("value").sum()).collect(streaming=True)
```

## 最佳实践

### 格式选择

**在以下情况下使用镶木地板：**
- 需要压缩（最多比 CSV 小 10 倍）
- 想要快速读/写
- 需要保留数据类型
- 处理大型数据集
- 需要谓词下推

**在以下情况下使用 CSV：**
- 需要人类可读的格式
- 与遗留系统接口
- 数据量小
- 需要通用兼容性

**在以下情况下使用 JSON：**
- 使用嵌套/分层数据
- 需要Web API兼容性
- 数据具有灵活的模式

**在以下情况下使用 Arrow IPC：**
- 需要零拷贝数据共享
- 需要最快的序列化
- 在 Arrow 兼容系统之间工作

### 读取大文件

```python
# 1. Always use lazy mode
lf = pl.scan_csv("large.csv")  # NOT read_csv

# 2. Filter and select early (pushdown optimization)
result = (
    lf
    .select("col1", "col2", "col3")  # Only needed columns
    .filter(pl.col("date") > "2023-01-01")  # Filter early
    .collect()
)

# 3. Use streaming for very large data
result = lf.filter(...).select(...).collect(streaming=True)

# 4. Read only needed rows during development
df = pl.read_csv("large.csv", n_rows=10000)  # Sample for testing
```

### 写入大文件

```python
# 1. Use Parquet with compression
df.write_parquet("output.parquet", compression="zstd")

# 2. Use partitioning for very large datasets
df.write_parquet("output", partition_by=["year", "month"])

# 3. Write streaming
lf = pl.scan_csv("input.csv")
lf.sink_parquet("output.parquet")  # Streaming write
```

### 性能提示

```python
# 1. Specify dtypes when reading CSV
df = pl.read_csv(
    "data.csv",
    dtypes={"id": pl.Int64, "name": pl.Utf8}  # Avoids inference
)

# 2. Use appropriate compression
df.write_parquet("output.parquet", compression="snappy")  # Fast
df.write_parquet("output.parquet", compression="zstd")    # Better compression

# 3. Parallel reading
df = pl.read_csv("data.csv", parallel="auto")

# 4. Read multiple files in parallel
lf = pl.scan_parquet("data/*.parquet")  # Automatic parallel read
```

## 错误处理

```python
try:
    df = pl.read_csv("data.csv")
except pl.exceptions.ComputeError as e:
    print(f"Error reading CSV: {e}")

# Ignore errors during parsing
df = pl.read_csv("messy.csv", ignore_errors=True)

# Handle missing files
from pathlib import Path
if Path("data.csv").exists():
    df = pl.read_csv("data.csv")
else:
    print("File not found")
```

## 模式管理

```python
# Infer schema from sample
schema = pl.read_csv("data.csv", n_rows=1000).schema

# Use inferred schema for full read
df = pl.read_csv("data.csv", dtypes=schema)

# Define schema explicitly
schema = {
    "id": pl.Int64,
    "name": pl.Utf8,
    "date": pl.Date,
    "value": pl.Float64
}
df = pl.read_csv("data.csv", dtypes=schema)
```