<!-- 此文件由机器翻译自 usage_examples.md -->

# 使用示例和工作流程

## 完整的工作流程示例

### 示例 1：会议演示包

**场景**：准备一个包含网站、海报和视频的大型会议演示。

**用户请求**：“我需要为我的 NeurIPS 论文提交创建一个完整的演示文稿包。生成一个网站、海报和视频演示文稿。”

**工作流程**：

```bash
# Step 1: Organize paper files
mkdir -p input/neurips2025_paper
cp main.tex input/neurips2025_paper/
cp -r figures/ input/neurips2025_paper/
cp -r tables/ input/neurips2025_paper/
cp bibliography.bib input/neurips2025_paper/

# Step 2: Generate all components
python pipeline_all.py \
  --input-dir input/neurips2025_paper \
  --output-dir output/ \
  --model-choice 1 \
  --generate-website \
  --generate-poster \
  --generate-video \
  --poster-width-inches 48 \
  --poster-height-inches 36 \
  --enable-logo-search

# Step 3: Review outputs
ls -R output/neurips2025_paper/
# - website/index.html
# - poster/poster_final.pdf
# - video/final_video.mp4
```

**输出**：
- 展示研究的互动网站
- 4'×3'会议海报（可打印）
- 12 分钟的演示视频
- 处理时间：约 45 分钟（无头部说话）

---

### 示例 2：预印本快速网站

**场景**：为 bioRxiv 预印本创建一个可探索的主页。

**用户请求**：“将我的基因组学预印本转换为交互式网站，以伴随 bioRxiv 提交。”

**工作流程**：

<<<代码块_1>>>

**提示**：
- 包含 bioRxiv DOI 的链接
- 添加 GitHub 存储库链接
- 包括数据可用性部分
- 如果可能的话，嵌入交互式可视化

---

### 示例 3：用于期刊提交的视频摘要

**场景**：为鼓励多媒体提交的期刊创建视频摘要。

**用户请求**：“为我提交的《自然通讯》生成一个 5 分钟的视频摘要。”

**工作流程**：

<<<代码块_2>>>

**输出**：
- 5分钟视频摘要
- 注重视觉效果
- 清晰易懂的旁白
- 期刊就绪格式

---

### 示例 4：多纸网站生成

**场景**：为研究小组的多篇论文创建网站。

**用户请求**：“为我们实验室今年发表的所有 5 篇论文生成网站。”

**工作流程**：

<<<代码块_3>>>

**最佳实践**：
- 使用一致的命名约定
- 大批量加工过夜
- 检查每个网站的准确性
- 部署到统一实验室网站

---

### 示例 5：虚拟会议海报

**场景**：为具有交互元素的虚拟会议创建数字海报。

**用户请求**：“为虚拟 ISMB 会议创建海报，其中包含指向代码和数据的可点击链接。”

**工作流程**：

<<<代码块_4>>>

**数字增强**：
- 带有嵌入超链接的 PDF
- 适用于虚拟平台的高分辨率 PNG
- 带有视频链接的单独 PDF 可供下载

---

### 示例 6：宣传视频剪辑

**场景**：为社交媒体制作简短的宣传视频。

**用户请求**：“为 Twitter 生成我们的 Cell 论文的 2 分钟精彩视频。”

**工作流程**：

<<<代码块_5>>>

**社交媒体优化**：
- 适用于 Instagram 的方形格式 (1:1)
- Twitter/LinkedIn 的水平格式 (16:9)
- TikTok/故事的垂直格式 (9:16)
- 为关键发现添加文本叠加

---

## 常见用例模式

###图案1：乳胶纸→全包

**输入**：包含所有资源的 LaTeX 源
**输出**：网站+海报+视频
**时间**：45-90 分钟
**最适合**：主要出版物、会议演示

<<<代码块_6>>>

---

### 模式 2：PDF → 互动网站

**输入**：已发布的 PDF 论文
**输出**：可探索的网站
**时间**：15-30分钟
**最适合**：出版后推广、预印本增强

```bash
python pipeline_all.py \
  --input-dir [pdf_dir] \
  --output-dir [output_dir] \
  --model-choice 1 \
  --generate-website
```

---

### 模式 3：LaTeX → 会议海报

**输入**：乳胶纸
**输出**：打印海报（自定义尺寸）
**时间**：10-20分钟
**最适合**：会议海报会议

```bash
python pipeline_all.py \
  --input-dir [latex_dir] \
  --output-dir [output_dir] \
  --model-choice 1 \
  --generate-poster \
  --poster-width-inches [width] \
  --poster-height-inches [height]
```

---

### 模式 4：LaTeX → 演示视频

**输入**：乳胶纸
**输出**：解说演示视频
**时间**：20-60分钟（不含头部说话）
**最适合**：视频摘要、在线演示、课程材料

```bash
python pipeline_light.py \
  --model_name_t gpt-4.1 \
  --model_name_v gpt-4.1 \
  --result_dir [output_dir] \
  --paper_latex_root [latex_dir]
```

---

## 特定于平台的输出

### Twitter/X 促销内容

系统自动检测 Twitter 目标的数字文件夹名称：

```bash
# Create Twitter-optimized content
mkdir -p input/001_twitter_post/
# System generates English promotional content
```

**生成的输出**：
- 简短、引人入胜的摘要
- 关键人物亮点
- 标签推荐
- 线程就绪格式

---

###小红书内容

对于中文社交媒体，请使用字母数字文件夹名称：

```bash
# Create Xiaohongshu-optimized content
mkdir -p input/xhs_genomics/
# System generates Chinese promotional content
```

**生成的输出**：
- 中文内容
- 适合平台的格式
- 视觉优先的演示
- 参与度优化

---

## 常见场景故障排除

### 场景：大纸张（>50 页）

**挑战**：处理时间和内容选择
**解决方案**：
```bash
# Option 1: Focus on key sections
# Edit LaTeX to comment out less critical sections

# Option 2: Process in parts
# Generate website for overview
# Generate separate detailed videos for methods/results

# Option 3: Use faster model for initial pass
# Review and regenerate critical components with better model
```

---

### 场景：复杂的数学内容

**挑战**：方程可能无法完美呈现
**解决方案**：
- 使用 LaTeX 输入（而非 PDF）实现最佳方程处理
- 检查生成的内容以确保方程的准确性
- 如果需要，手动调整复杂的方程
- 考虑对关键方程使用图形屏幕截图

---

### 场景：非标准纸张结构

**挑战**：论文不遵循标准 IMRAD 格式
**解决方案**：
- 在论文元数据中提供自定义部分指导
- 检查生成的结构并进行调整
- 使用更强大的模型（GPT-4.1）以获得更好的适配
- 考虑在 LaTeX 注释中手动注释部分

---

### 场景：API 预算有限

**挑战**：在保持质量的同时降低成本
**解决方案**：
```bash
# Use GPT-3.5-turbo for simple papers
python pipeline_all.py \
  --input-dir [paper_dir] \
  --output-dir [output_dir] \
  --model-choice 3

# Generate only needed components
# Website-only (cheapest)
# Poster-only (moderate)
# Video without talking-head (moderate)
```

---

### 场景：截止日期紧迫

**挑战**：需要快速输出
**解决方案**：
```bash
# Parallel processing if multiple papers
# Use faster models (GPT-3.5-turbo)
# Generate only essential component first
# Skip optional features (logo search, talking-head)

python pipeline_light.py \
  --model_name_t gpt-3.5-turbo \
  --model_name_v gpt-3.5-turbo \
  --result_dir [output_dir] \
  --paper_latex_root [latex_dir]
```

**优先顺序**：
1. 网站（最快、最通用）
2.海报（速度适中，打印截止时间）
3.视频（最慢，可以稍后生成）

---

## 质量优化技巧

### 获得最佳网站结果
1. 对所有资源使用 LaTeX 输入
2. 包含高分辨率图形
3. 确保纸张具有清晰的截面结构
4.启用徽标搜索以获得专业外观
5.审查并测试所有交互元素

### 最佳海报效果
1.提供高分辨率图形（300+DPI）
2. 指定所需的确切海报尺寸
3. 包含机构品牌信息
4.使用专业的配色方案
5. 在完整海报之前测试打印小预览

### 获得最佳视频效果
1.使用LaTeX进行最清晰的内容提取
2. 适当指定目标持续时间
3. 视频生成前审查脚本
4.选择合适的展示风格
5. 测试音频质量和节奏

### 获得最佳总体结果
1. 从干净、组织良好的 LaTeX 源开始
2. 使用 GPT-4 或 GPT-4.1 以获得最高质量
3. 在最终确定之前审查所有输出
4. 迭代任何需要调整的组件
5. 组合组件以获得有凝聚力的演示包