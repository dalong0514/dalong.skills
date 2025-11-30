<!-- 此文件由机器翻译自 advanced_integrations.md -->

# 高级集成参考

本文档提供有关 MarkItDown 高级功能的详细信息，包括 Azure 文档智能集成、LLM 支持的描述和插件系统。

## Azure 文档智能集成

Azure 文档智能（以前称为表单识别器）通过高级表格提取和布局分析提供卓越的 PDF 处理。

### 设置

**先决条件：**
1.Azure订阅
2.在Azure中创建的文档智能资源
3. 端点 URL 和 API 密钥

**创建 Azure 资源：**
```bash
# Using Azure CLI
az cognitiveservices account create \
  --name my-doc-intelligence \
  --resource-group my-resource-group \
  --kind FormRecognizer \
  --sku F0 \
  --location eastus
```

### 基本用法

<<<代码块_1>>>

### 来自环境变量的配置

<<<代码块_2>>>

### 何时使用 Azure 文档智能

**用于：**
- 带有复杂表格的复杂 PDF
- 多列布局
- 表格和结构化文件
- 需要OCR的扫描文档
- 具有混合内容类型的 PDF
- 格式复杂的文档

**相对于标准提取的优点：**
- **卓越的表格提取** - 更好地处理合并单元格、复杂布局
- **布局分析** - 了解文档结构（页眉、页脚、列）
- **表单字段** - 从表单中提取键值对
- **阅读顺序** - 在复杂的布局中保持正确的文本流
- **OCR 质量** - 从扫描文档中提取高质量文本

### 比较示例

**标准萃取：**
<<<代码块_3>>>

**Azure 文档智能：**
<<<代码块_4>>>

### 成本考虑

Azure 文档智能是一项付费服务：
- **免费套餐**：每月 500 页
- **付费等级**：按处理页面付费
- 监控使用情况以控制成本
- 对简单文档使用标准提取

### 错误处理

<<<代码块_5>>>

## LLM 支持的图像描述

使用大型语言模型为图像生成详细的上下文描述。

### 使用 OpenAI 设置

<<<代码块_6>>>

### 支持的用例

**文档中的图像：**
```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()
md = MarkItDown(llm_client=client, llm_model="gpt-4o")

# PowerPoint with images
result = md.convert("presentation.pptx")

# Word documents with images
result = md.convert("report.docx")

# Standalone images
result = md.convert("diagram.png")
```

### 自定义提示

根据特定需求自定义 LLM 提示：

```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()

# For diagrams
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    llm_prompt="Analyze this diagram and explain all components, connections, and relationships in detail"
)

# For charts
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    llm_prompt="Describe this chart, including the type, axes, data points, trends, and key insights"
)

# For UI screenshots
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    llm_prompt="Describe this user interface screenshot, listing all UI elements, their layout, and functionality"
)

# For scientific figures
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    llm_prompt="Describe this scientific figure in detail, including methodology, results shown, and significance"
)
```

### 型号选择

**GPT-4o（推荐）：**
- 最佳视觉能力
- 高质量的描述
- 善于理解上下文
- 每张图像的成本更高

**GPT-4o-迷你：**
- 成本更低的替代方案
- 适用于更简单的图像
- 处理速度更快
- 可能会错过微妙的细节

```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()

# High quality (more expensive)
md_quality = MarkItDown(llm_client=client, llm_model="gpt-4o")

# Budget option (less expensive)
md_budget = MarkItDown(llm_client=client, llm_model="gpt-4o-mini")
```

### 环境配置

```python
import os
from markitdown import MarkItDown
from openai import OpenAI

# Set API key in environment
os.environ['OPENAI_API_KEY'] = 'YOUR-API-KEY'

client = OpenAI()  # Uses env variable
md = MarkItDown(llm_client=client, llm_model="gpt-4o")
```

### 替代法学硕士提供商

**人类克劳德：**
```python
from markitdown import MarkItDown
from anthropic import Anthropic

# Note: Check current compatibility with MarkItDown
client = Anthropic(api_key="YOUR-API-KEY")
# May require adapter for MarkItDown compatibility
```

**Azure OpenAI：**
```python
from markitdown import MarkItDown
from openai import AzureOpenAI

client = AzureOpenAI(
    api_key="YOUR-AZURE-KEY",
    api_version="2024-02-01",
    azure_endpoint="https://YOUR-RESOURCE.openai.azure.com"
)

md = MarkItDown(llm_client=client, llm_model="gpt-4o")
```

### 成本管理

**降低法学硕士成本的策略：**

1. **选择性加工：**
```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()

# Only use LLM for important documents
if is_important_document(file):
    md = MarkItDown(llm_client=client, llm_model="gpt-4o")
else:
    md = MarkItDown()  # Standard processing

result = md.convert(file)
```

2. **图像过滤：**
```python
# Pre-process to identify images that need descriptions
# Only use LLM for complex/important images
```

3. **批量处理：**
```python
# Process multiple images in batches
# Monitor costs and set limits
```

4. **型号选择：**
```python
# Use gpt-4o-mini for simple images
# Reserve gpt-4o for complex visualizations
```

### 性能考虑因素

**LLM 处理会增加延迟：**
- 每个图像都需要一个API调用
- 处理时间：每张图像 1-5 秒
- 依赖网络
- 考虑对多个图像进行并行处理

**批量优化：**
```python
from markitdown import MarkItDown
from openai import OpenAI
import concurrent.futures

client = OpenAI()
md = MarkItDown(llm_client=client, llm_model="gpt-4o")

def process_image(image_path):
    return md.convert(image_path)

# Process multiple images in parallel
images = ["img1.jpg", "img2.jpg", "img3.jpg"]
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
    results = list(executor.map(process_image, images))
```

## 组合的高级功能

### Azure 文档智能 + LLM 描述

结合两者以获得最佳质量：

```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    docintel_endpoint="YOUR-AZURE-ENDPOINT",
    docintel_key="YOUR-AZURE-KEY"
)

# Best possible PDF conversion with image descriptions
result = md.convert("complex_report.pdf")
```

**使用案例：**
- 带图表的研究论文
- 带图表的业务报告
- 带图表的技术文档
- 带有视觉数据的演示

### 智能文档处理管道

```python
from markitdown import MarkItDown
from openai import OpenAI
import os

def smart_convert(file_path):
    """Intelligently choose processing method based on file type."""
    client = OpenAI()
    ext = os.path.splitext(file_path)[1].lower()

    # PDFs with complex tables: Use Azure
    if ext == '.pdf':
        md = MarkItDown(
            docintel_endpoint=os.getenv('AZURE_ENDPOINT'),
            docintel_key=os.getenv('AZURE_KEY')
        )

    # Documents/presentations with images: Use LLM
    elif ext in ['.pptx', '.docx']:
        md = MarkItDown(
            llm_client=client,
            llm_model="gpt-4o"
        )

    # Simple formats: Standard processing
    else:
        md = MarkItDown()

    return md.convert(file_path)

# Use it
result = smart_convert("document.pdf")
```

## 插件系统

MarkItDown 支持自定义插件来扩展功能。

### 插件架构

为了安全起见，插件默认被禁用：

```python
from markitdown import MarkItDown

# Enable plugins
md = MarkItDown(enable_plugins=True)
```

### 创建自定义插件

**插件结构：**
```python
class CustomConverter:
    """Custom converter plugin for MarkItDown."""

    def can_convert(self, file_path):
        """Check if this plugin can handle the file."""
        return file_path.endswith('.custom')

    def convert(self, file_path):
        """Convert file to Markdown."""
        # Your conversion logic here
        return {
            'text_content': '# Converted Content\n\n...'
        }
```

### 插件注册

```python
from markitdown import MarkItDown

md = MarkItDown(enable_plugins=True)

# Register custom plugin
md.register_plugin(CustomConverter())

# Use normally
result = md.convert("file.custom")
```

### 插件用例

**自定义格式：**
- 专有文档格式
- 专门的科学数据格式
- 旧文件格式

**增强处理：**
- 自定义 OCR 引擎
- 专门的表提取
- 特定领域的解析

**整合：**
- 企业文档系统
- 自定义数据库
- 专门的API

### 插件安全

**重要的安全考虑：**
- 插件以完全系统访问权限运行
- 仅对受信任的插件启用
- 使用前验证插件代码
- 除非需要，否则禁用生产中的插件

## 高级功能的错误处理

```python
from markitdown import MarkItDown
from openai import OpenAI

def robust_convert(file_path):
    """Convert with fallback strategies."""
    try:
        # Try with all advanced features
        client = OpenAI()
        md = MarkItDown(
            llm_client=client,
            llm_model="gpt-4o",
            docintel_endpoint=os.getenv('AZURE_ENDPOINT'),
            docintel_key=os.getenv('AZURE_KEY')
        )
        return md.convert(file_path)

    except Exception as azure_error:
        print(f"Azure failed: {azure_error}")

        try:
            # Fallback: LLM only
            client = OpenAI()
            md = MarkItDown(llm_client=client, llm_model="gpt-4o")
            return md.convert(file_path)

        except Exception as llm_error:
            print(f"LLM failed: {llm_error}")

            # Final fallback: Standard processing
            md = MarkItDown()
            return md.convert(file_path)

# Use it
result = robust_convert("document.pdf")
```

## 最佳实践

### Azure 文档智能
- 仅用于复杂的 PDF（成本优化）
- 监控使用情况和成本
- 安全地存储凭据
- 妥善处理配额限制
- 如果需要，回退到标准处理
### 法学硕士整合
- 针对任务复杂性使用适当的模型
- 针对特定用例自定义提示
- 监控API成本
- 实施速率限制
- 尽可能缓存结果
- 优雅地处理 API 错误

### 组合功能
- 测试成本/质量权衡
- 有选择地使用重要文件
- 实现智能路由
- 监控绩效和成本
- 有后备策略

### 安全
- 安全地存储 API 密钥（环境变量、秘密管理器）
- 切勿将凭据提交给代码
- 除非需要，否则禁用插件
- 验证所有输入
- 使用最小权限访问