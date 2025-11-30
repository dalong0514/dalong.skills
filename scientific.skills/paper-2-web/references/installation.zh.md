<!-- 此文件由机器翻译自 installation.md -->

# 安装和配置

## 系统要求

### 硬件要求
- **GPU**：具有头部说话功能的视频生成需要 NVIDIA A6000（最低 48GB）
- **CPU**：推荐用于 PDF 处理和文档转换的多核处理器
- **RAM**：最低 16GB，对于大型纸张建议 32GB

### 软件要求
- **Python**：3.11 或更高版本
- **LibreOffice**：文档格式转换所需（PDF 到 PPTX 等）
- **Poppler 实用程序**：PDF 处理和操作所需

## 安装步骤

### 1.克隆存储库
```bash
git clone https://github.com/YuhangChen1/Paper2All.git
cd Paper2All
```

### 2.安装依赖项
<<<代码块_1>>>

### 3.安装系统依赖项

**Ubuntu/Debian：**
<<<代码块_2>>>

**苹果系统：**
<<<代码块_3>>>

**Windows：**
- 从 https://www.libreoffice.org/ 下载并安装 LibreOffice
- 从 https://github.com/oschwartz10612/poppler-windows 下载并安装 Poppler

## API 配置

使用以下凭据在项目根目录中创建一个 `.env` 文件：

### 所需的 API 密钥

**选项 1：OpenAI API**
<<<代码块_4>>>

**选项 2：OpenRouter API**（OpenAI 的替代方案）
<<<代码块_5>>>

### 可选 API 密钥

**Google 搜索 API**（用于自动徽标发现）
<<<代码块_6>>>

## 型号配置

系统支持多个LLM后端：

### 支持的型号
- GPT-4（推荐最佳质量）
- GPT-4.1（最新版本）
- GPT-3.5-turbo（更快，成本更低）
- 通过 OpenRouter 的 Claude 模型
- 其他 OpenRouter 支持的型号

### 型号选择

使用 `--model-choice` 参数或 `--model_name_t` 和 `--model_name_v` 参数指定模型：
- 型号选择 1：适用于所有组件的 GPT-4
- 型号选择 2：适用于所有组件的 GPT-4.1
- 自定义：为文本和视觉处理指定单独的模型

## 验证

测试安装：

```bash
python pipeline_all.py --help
```

如果成功，您应该会看到包含所有可用选项的帮助菜单。

## 故障排除

### 常见问题

**1.未找到 LibreOffice**
- 确保 LibreOffice 已安装并位于您的系统路径中
- 尝试运行`libreoffice --version`来验证

**2.未找到 Poppler 实用程序**
- 使用 `pdftoppm -v` 验证安装
- 如果需要，将 Poppler bin 目录添加到 PATH

**3.视频生成的 GPU/CUDA 错误**
- 确保 NVIDIA 驱动程序是最新的
- 验证 CUDA 工具包是否已安装
- 使用 `nvidia-smi` 检查 GPU 内存

**4. API 密钥错误**
- 验证 `.env` 文件是否位于项目根目录中
- 检查API密钥是否有效并且有足够的积分
- 确保 `.env` 中的键周围没有多余的空格或引号

## 目录结构

安装后，整理您的工作空间：

```
Paper2All/
├── .env                  # API credentials
├── input/               # Place your paper files here
│   └── paper_name/      # Each paper in its own directory
│       └── main.tex     # LaTeX source or PDF
├── output/              # Generated outputs
│   └── paper_name/
│       ├── website/     # Generated website files
│       ├── video/       # Generated video files
│       └── poster/      # Generated poster files
└── ...
```