<!-- 此文件由机器翻译自 media_processing.md -->

# 媒体处理参考

本文档提供有关使用 MarkItDown 处理图像和音频文件的详细信息。

## 图像处理

MarkItDown 可以使用 OCR 从图像中提取文本并检索 EXIF 元数据。

### 基本图像转换

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("photo.jpg")
print(result.text_content)
```

### 图像处理功能

**提取的内容：**
1. **EXIF 元数据** - 相机设置、日期、位置等。
2. **OCR 文本** - 在图像中检测到的文本（需要 tesseract）
3. **图像描述** - AI生成的描述（与LLM集成）

### EXIF 元数据提取

来自相机和智能手机的图像包含自动提取的 EXIF 元数据：

<<<代码块_1>>>

**示例输出包括：**
- 相机品牌和型号
- 捕获日期和时间
- GPS 坐标（如果有）
- 曝光设置（ISO、快门速度、光圈）
- 图像尺寸
- 方向

### OCR（光学字符识别）

从包含文本的图像中提取文本（屏幕截图、扫描文档、文本照片）：

**要求：**
- 安装tesseract OCR引擎：
  <<<代码块_2>>>

**用途：**
<<<代码块_3>>>

**OCR 最佳实践：**
- 使用高分辨率图像以获得更好的准确性
- 确保文本和背景之间有良好的对比度
- 如果可能的话，拉直倾斜的文本
- 使用光线充足、清晰的图像

### LLM 生成的图像描述

使用 GPT-4o 或其他视觉模型生成图像的详细上下文描述：

<<<代码块_4>>>

**针对特定需求的自定义提示：**

<<<代码块_5>>>

### 支持的图像格式

MarkItDown 支持所有常见图像格式：
- JPEG/JPG
- 巴布亚新几内亚
- 动图
-骨形态发生蛋白
- TIFF
- 网页版
- HEIC（在某些平台上需要额外的库）

## 音频处理

MarkItDown 可以使用语音识别将音频文件转录为文本。

### 基本音频转换

<<<代码块_6>>>

### 音频转录设置

**安装：**
```bash
pip install 'markitdown[audio]'
```

这将安装 `speech_recognition` 库和依赖项。

### 支持的音频格式

- 音频
- 亚洲国际电影节
- FLAC
- MP3（需要 ffmpeg 或 libav）
- OGG（需要 ffmpeg 或 libav）
- 语音识别支持的其他格式

### 音频转录引擎

MarkItDown 使用 `speech_recognition` 库，它支持多个后端：

**默认（谷歌语音识别）：**
```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("audio.wav")
```

**注意：** 默认 Google 语音识别需要互联网连接。

### 音频质量注意事项

为了获得最佳转录准确性：
- 使用清晰的音频和最小的背景噪音
- 更喜欢 WAV 或 FLAC 以获得更好的质量
- 确保讲话清晰且音量良好
- 避免多个扬声器重叠
- 尽可能使用单声道音频

### 音频预处理技巧

为了获得更好的结果，请考虑预处理音频：

```python
# Example: If you have pydub installed
from pydub import AudioSegment
from pydub.effects import normalize

# Load and normalize audio
audio = AudioSegment.from_file("recording.mp3")
audio = normalize(audio)
audio.export("normalized.wav", format="wav")

# Then convert with MarkItDown
from markitdown import MarkItDown
md = MarkItDown()
result = md.convert("normalized.wav")
```

## 组合媒体工作流程

### 批量处理多个图像

```python
from markitdown import MarkItDown
from openai import OpenAI
import os

client = OpenAI()
md = MarkItDown(llm_client=client, llm_model="gpt-4o")

# Process all images in directory
for filename in os.listdir("images"):
    if filename.lower().endswith(('.png', '.jpg', '.jpeg')):
        result = md.convert(f"images/{filename}")

        # Save markdown with same name
        output = filename.rsplit('.', 1)[0] + '.md'
        with open(f"output/{output}", "w") as f:
            f.write(result.text_content)
```

### 截图分析管道

```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    llm_prompt="Describe this screenshot comprehensively, including UI elements, text, and layout"
)

screenshots = ["screen1.png", "screen2.png", "screen3.png"]
analysis = []

for screenshot in screenshots:
    result = md.convert(screenshot)
    analysis.append({
        'file': screenshot,
        'content': result.text_content
    })

# Now ready for further processing
```

### 使用 OCR 记录图像

对于扫描文档或文档照片：

```python
from markitdown import MarkItDown

md = MarkItDown()

# Process scanned pages
pages = ["page1.jpg", "page2.jpg", "page3.jpg"]
full_text = []

for page in pages:
    result = md.convert(page)
    full_text.append(result.text_content)

# Combine into single document
document = "\n\n---\n\n".join(full_text)
print(document)
```

### 演示幻灯片图像

当您将演示幻灯片作为图像时：

```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    llm_prompt="Describe this presentation slide, including title, bullet points, and visual elements"
)

# Process slide images
for i in range(1, 21):  # 20 slides
    result = md.convert(f"slides/slide_{i}.png")
    print(f"## Slide {i}\n\n{result.text_content}\n\n")
```

## 错误处理

### 图像处理错误

```python
from markitdown import MarkItDown

md = MarkItDown()

try:
    result = md.convert("image.jpg")
    print(result.text_content)
except FileNotFoundError:
    print("Image file not found")
except Exception as e:
    print(f"Error processing image: {e}")
```

### 音频处理错误

```python
from markitdown import MarkItDown

md = MarkItDown()

try:
    result = md.convert("audio.mp3")
    print(result.text_content)
except Exception as e:
    print(f"Transcription failed: {e}")
    # Common issues: format not supported, no speech detected, network error
```

## 性能优化

### 图像处理

- **法学硕士描述**：速度较慢但信息更丰富
- **仅限 OCR**：文本提取速度更快
- **仅 EXIF**：最快，仅元数据
- **批处理**：并行处理多个图像

### 音频处理

- **文件大小**：较大的文件需要更长的时间
- **音频长度**：转录时间尺度与持续时间
- **格式转换**：WAV/FLAC 比 MP3/OGG 更快
- **网络依赖性**：默认转录需要互联网

## 用例

### 文档数字化
将扫描文档或文档照片转换为可搜索文本。

### 会议记录
将会议录音转录为文本以进行分析。

### 演示分析
从演示幻灯片图像中提取内容。

### 截图文档
生成 UI 屏幕截图的描述以供文档使用。

### 图像存档
从照片集中提取元数据和内容。

### 辅助功能
使用 LLM 集成生成图像的替代文本描述。

### 数据提取
来自包含表格、表单或结构化数据的图像的 OCR 文本。