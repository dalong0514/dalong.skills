<!-- 此文件由机器翻译自 data-access.md -->

# DrugBank 数据访问

## 身份验证和设置

### 帐户创建
DrugBank 需要用户身份验证才能访问数据：
1. 在 go.drugbank.com 创建帐户
2.接受许可协议（学术使用免费，商业使用付费）
3. 获取用户名和密码凭证

### 凭证管理

**环境变量（推荐）**
```bash
export DRUGBANK_USERNAME="your_username"
export DRUGBANK_PASSWORD="your_password"
```

**配置文件**
创建`~/.config/drugbank.ini`：
<<<代码块_1>>>

**直接规格**
<<<代码块_2>>>

## Python 包安装

### drugbank-下载器
编程访问的主要工具：
<<<代码块_3>>>

**要求：** Python >=3.9

### 可选依赖项
<<<代码块_4>>>

## 数据下载方法

### 下载完整数据库
<<<代码块_5>>>

### 自定义存储位置
<<<代码块_6>>>

### 验证下载
```python
import os
if os.path.exists(path):
    size_mb = os.path.getsize(path) / (1024 * 1024)
    print(f"Downloaded successfully: {size_mb:.1f} MB")
```

## 使用下载的数据

### 打开压缩的 XML 而不解压
```python
from drugbank_downloader import open_drugbank
import xml.etree.ElementTree as ET

# Open file directly from zip
with open_drugbank() as file:
    tree = ET.parse(file)
    root = tree.getroot()
```

### 解析 XML 树
```python
from drugbank_downloader import parse_drugbank, get_drugbank_root

# Get parsed tree
tree = parse_drugbank()

# Get root element directly
root = get_drugbank_root()
```

### CLI 用法
```bash
# Download using command line
drugbank_downloader --username USER --password PASS

# Download latest version
drugbank_downloader
```

## 数据格式和版本

### 可用格式
- **XML**：主要格式，最全面的数据
- **JSON**：通过 API 提供（需要单独的 API 密钥）
- **CSV/TSV**：从 Web 界面导出或解析 XML
- **SQL**：数据库转储可供下载

### 版本管理
```python
# Specify exact version for reproducibility
path = download_drugbank(version='5.1.10')

# List cached versions
from pathlib import Path
drugbank_dir = Path.home() / '.data' / 'drugbank'
if drugbank_dir.exists():
    versions = [d.name for d in drugbank_dir.iterdir() if d.is_dir()]
    print(f"Cached versions: {versions}")
```

### 版本历史
- **版本 6.0** (2024)：最新版本，扩展药物条目
- **版本 5.1.x** (2019-2023)：增量更新
- **版本 5.0** (2017)：约 9,591 个药物条目
- **版本 4.0** (2014)：添加了代谢物结构
- **版本 3.0** (2011)：添加了转运蛋白和途径数据
- **版本 2.0** (2009)：添加了交互和 ADMET

## API 访问

### REST API 端点
```python
import requests

# Query by DrugBank ID
drug_id = "DB00001"
url = f"https://go.drugbank.com/drugs/{drug_id}.json"
headers = {"Authorization": "Bearer YOUR_API_KEY"}

response = requests.get(url, headers=headers)
if response.status_code == 200:
    drug_data = response.json()
```

### 速率限制
- **开发关键**：3,000 个请求/月
- **生产密钥**：基于许可证的自定义限制
- **最佳实践**：在本地缓存结果以最大程度地减少 API 调用

### 区域范围界定
DrugBank API 按地区划分：
- **美国**：FDA 批准的药物
- **加拿大**：加拿大卫生部批准的药物
- **EU**：EMA 批准的药物

如果适用，请在 API 请求中指定区域。

## 数据缓存策略

### 中间结果
```python
import pickle
from pathlib import Path

# Cache parsed data
cache_file = Path("drugbank_parsed.pkl")

if cache_file.exists():
    with open(cache_file, 'rb') as f:
        data = pickle.load(f)
else:
    # Parse and process
    root = get_drugbank_root()
    data = process_drugbank_data(root)

    # Save cache
    with open(cache_file, 'wb') as f:
        pickle.dump(data, f)
```

### 版本特定的缓存
```python
version = "5.1.10"
cache_file = Path(f"drugbank_{version}_processed.pkl")
# Ensures cache invalidation when version changes
```

## 故障排除

### 常见问题

**身份验证失败**
- 验证凭据是否正确
- 检查许可协议是否被接受
- 确保帐户未过期

**下载失败**
- 检查互联网连接
- 验证足够的磁盘空间（需要约 1-2 GB）
- 如果最新版本失败，请尝试指定旧版本

**解析错误**
- 确保完整下载（检查文件大小）
- 验证 XML 未损坏
- 使用 lxml 解析器更好地处理错误

### 错误处理
```python
from drugbank_downloader import download_drugbank
import logging

logging.basicConfig(level=logging.INFO)

try:
    path = download_drugbank()
    print(f"Success: {path}")
except Exception as e:
    print(f"Download failed: {e}")
    # Fallback: specify older stable version
    path = download_drugbank(version='5.1.7')
```

## 最佳实践

1. **版本规范**：始终指定准确的版本以进行可重复的研究
2. **凭证安全**：使用环境变量，切勿硬编码凭证
3. **缓存**：缓存中间处理结果，避免重新解析
4. **文档**：记录分析中使用的 DrugBank 版本
5. **许可证合规性**：确保您的用例获得适当的许可
6. **本地存储**：保留本地副本以减少下载频率
7. **错误处理**：针对网络问题实施强大的错误处理