<!-- 此文件由机器翻译自 installation.md -->

# ToolUniverse 安装和设置

## 安装

```bash
uv pip install tooluniverse
```

## 基本设置

### Python SDK
<<<代码块_1>>>

## 模型上下文协议 (MCP) 设置

ToolUniverse 提供本机 MCP 支持，可与 Claude Desktop、Claude Code 和其他 MCP 兼容系统集成。

### 启动 MCP 服务器
<<<代码块_2>>>

这将启动一个 MCP 服务器，该服务器通过模型上下文协议公开 ToolUniverse 的 600 多个工具。

### 克劳德桌面集成

添加到 Claude 桌面配置 (~/.config/Claude/claude_desktop_config.json)：
<<<代码块_3>>>

### 克劳德代码集成

ToolUniverse MCP 服务器通过 MCP 协议与 Claude Code 本地工作。

## 与其他平台集成

### OpenRouter 集成
ToolUniverse 与 OpenRouter 集成，可通过单个 API 访问 100 多个 LLM：
- GPT-5，克劳德，双子座
- Qwen，Deepseek
- 开源模型

### 支持的 LLM 平台
- 克劳德桌面和克劳德代码
- 双子座 CLI
- Qwen 代码
- ChatGPT API
- GPT Codex CLI

## 要求

-Python 3.8+
- 对于 Tool_Finder（基于嵌入的搜索）：推荐使用 GPU
- 对于 Tool_Finder_LLM：不需要 GPU（使用基于 LLM 的搜索）

## 验证

测试安装：
<<<代码块_4>>>