<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：pptx
描述：“演示文稿工具包 (.pptx)。创建/编辑幻灯片、布局、内容、演讲者注释、评论，以进行程序化演示文稿创建和修改。”
许可证：专有。 LICENSE.txt 有完整的条款
---

# PPTX 创建、编辑和分析

## 概述

.pptx 文件是包含 XML 文件和资源的 ZIP 存档。使用文本提取、原始 XML 访问或 html2pptx 工作流程创建、编辑或分析 PowerPoint 演示文稿。将此技能应用于程序化演示文稿创建和修改。

## 阅读并分析内容

### 文本提取
要阅读演示文稿的文本内容，请将文档转换为 Markdown：

```bash
# Convert document to markdown
python -m markitdown path-to-file.pptx
```

### 原始 XML 访问
需要原始 XML 访问：评论、演讲者注释、幻灯片布局、动画、设计元素和复杂格式。对于任何这些功能，请解压缩演示文稿并读取其原始 XML 内容。

#### 解压文件
`python ooxml/scripts/unpack.py <office_file> <output_dir>`

**注意**：unpack.py 脚本位于相对于项目根目录的 `skills/pptx/ooxml/scripts/unpack.py` 处。如果此路径中不存在该脚本，请使用 `find . -name "unpack.py"` 来查找它。

#### 关键文件结构
* `ppt/presentation.xml` - 主要演示文稿元数据和幻灯片引用
* `ppt/slides/slide{N}.xml` - 各个幻灯片内容（slide1.xml、slide2.xml 等）
* `ppt/notesSlides/notesSlide{N}.xml` - 每张幻灯片的演讲者注释
* `ppt/comments/modernComment_*.xml` - 针对特定幻灯片的评论
* `ppt/slideLayouts/` - 幻灯片的布局模板
* `ppt/slideMasters/` - 幻灯片母版模板
* `ppt/theme/` - 主题和样式信息
* `ppt/media/` - 图像和其他媒体文件

#### 版式和颜色提取
**当给出要模拟的示例设计时**：始终首先使用以下方法分析演示文稿的版式和颜色：
1. **读取主题文件**：检查 `ppt/theme/theme1.xml` 的颜色 (`<a:clrScheme>`) 和字体 (`<a:fontScheme>`)
2. **幻灯片内容示例**：检查 `ppt/slides/slide1.xml` 的实际字体使用情况 (`<a:rPr>`) 和颜色
3. **搜索模式**：使用 grep 在所有 XML 文件中查找颜色（`<a:solidFill>`、`<a:srgbClr>`）和字体引用

## 创建新的 PowerPoint 演示文稿 **无需模板**

从头开始创建新的 PowerPoint 演示文稿时，使用 **html2pptx** 工作流程将 HTML 幻灯片转换为具有精确定位的 PowerPoint。

### 设计原则

**关键**：在创建任何演示文稿之前，请分析内容并选择适当的设计元素：
1. **考虑主题**：本次演示的内容是什么？它暗示了什么基调、行业或情绪？
2. **检查品牌**：如果用户提到公司/组织，请考虑其品牌颜色和标识
3. **将调色板与内容相匹配**：选择反映主题的颜色
4. **陈述你的方法**：在编写代码之前解释你的设计选择

**要求**：
- ✅ 在编写代码之前说明您的内容知情设计方法
- ✅ 仅使用网络安全字体：Arial、Helvetica、Times New Roman、Georgia、Courier New、Verdana、Tahoma、Trebuchet MS、Impact
- ✅ 通过尺寸、重量和颜色创建清晰的视觉层次
- ✅ 确保可读性：强烈的对比度、适当大小的文本、干净的对齐方式
- ✅ 保持一致：在幻灯片中重复图案、间距和视觉语言

#### 调色板选择

**创造性地选择颜色**：
- **超越默认值**：什么颜色真正匹配这个特定主题？避免自动驾驶选择。
- **考虑多个角度**：主题、行业、情绪、能量水平、目标受众、品牌标识（如果提及）
- **勇于冒险**：尝试意想不到的组合 - 医疗保健演示不一定是绿色的，金融不一定是海军的
- **构建你的调色板**：选择 3-5 种协同工作的颜色（主色 + 支持色调 + 强调色）
- **确保对比度**：文本在背景上必须清晰可读

**示例调色板**（使用它们来激发创造力 - 选择一个、调整它或创建您自己的调色板）：

1. **经典蓝**：深海军蓝 (#1C2833)、板岩灰 (#2E4053)、银色 (#AAB7B8)、灰白色 (#F4F6F6)
2. **青色和珊瑚色**：青色 (#5EA8A7)、深青色 (#277884)、珊瑚色 (#FE4447)、白色 (#FFFFFF)
3. **粗红色**：红色 (#C0392B)、亮红色 (#E74C3C)、橙色 (#F39C12)、黄色 (#F1C40F)、绿色 (#2ECC71)
4. **暖色腮红**：紫红色 (#A49393)、腮红 (#EED6D3)、玫瑰色 (#E8B4B8)、奶油色 (#FAF7F2)
5. **勃艮第奢华**：勃艮第 (#5D1D2E)、深红色 (#951233)、铁锈色 (#C15937)、金色 (#997929)
6. **深紫色和翡翠色**：紫色 (#B165FB)、深蓝色 (#181B24)、翡翠色 (#40695B)、白色 (#FFFFFF)
7. **奶油色和森林绿**：奶油色 (#FFE1C7)、森林绿 (#40695B)、白色 (#FCFCFC)
8. **粉色和紫色**：粉色 (#F8275B)、珊瑚色 (#FF574A)、玫瑰色 (#FF737D)、紫色 (#3D2F68)
9. **青柠色和李子色**：青柠色 (#C5DE82)、李子色 (#7C3A5F)、珊瑚色 (#FD8C6E)、蓝灰色 (#98ACB5)
10. **黑色和金色**：金色 (#BF9A4A)、黑色 (#000000)、奶油色 (#F4F6F6)
11. **鼠尾草和赤土色**：鼠尾草（#87A96B）、赤土色（#E07A5F）、奶油色（#F4F1DE）、木炭（#2C2C2C）
12. **木炭色和红色**：木炭色 (#292929)、红色 (#E33737)、浅灰色 (#CCCBCB)
13. **活力橙色**：橙色 (#F96D00)、浅灰色 (#F2F2F2)、木炭色 (#222831)
14. **森林绿**：黑色(#191A19)、绿色(#4E9F3D)、深绿色(#1E5128)、白色(#FFFFFF)
15. **复古彩虹**：紫色(#722880)、粉色(#D72D51)、橙色(#EB5C18)、琥珀色(#F08800)、金色(#DEB600)
16. **复古土色**：芥末色 (#E3B448)、鼠尾草色 (#CBD18F)、森林绿 (#3A6B35)、奶油色 (#F4F1DE)
17. **海岸玫瑰**：老玫瑰色（#AD7670）、海狸色（#B49886）、蛋壳色（#F3ECDC）、灰灰色（#BFD5BE）
18. **橙色和绿松石色**：浅橙色 (#FC993E)、灰绿松石色 (#667C6F)、白色 (#FCFCFC)

#### 视觉细节选项

**几何图案**：
- 对角部分分隔线而不是水平分隔线
- 不对称列宽（30/70、40/60、25/75）
- 文本标题旋转 90° 或 270°
- 圆形/六角形图像框架
- 角落里的三角形强调形状
- 重叠形状的深度

**边框和框架处理**：
- 仅一侧有粗单色边框 (10-20pt)
- 对比色双线边框
- 角括号代替全框
- L 形边框（上+左或下+右）
- 标题下方的下划线重音（3-5pt 厚）

**版式处理**：
- 极端尺寸对比（72pt 标题 vs 11pt 正文）
- 标题全部大写，字母间距宽
- 超大显示类型的编号部分
- 等宽字体（Courier New）用于数据/统计/技术内容
- 压缩字体（Arial Narrow）用于提供密集信息
- 强调的概述文本

**图表和数据样式**：
- 单色图表，关键数据采用单一强调色
- 水平条形图而不是垂直条形图
- 点图而不是条形图
- 网格线最少或根本没有
- 数据标签直接在元素上（无图例）
- 关键指标的数字过大

**布局创新**：
- 带文本覆盖的全出血图像
- 用于导航/上下文的侧边栏（20-30% 宽度）
- 模块化网格系统（3×3、4×4块）
- Z 型或 F 型内容流
- 彩色形状上的浮动文本框
- 杂志风格的多栏布局

**背景处理**：
- 纯色块占据幻灯片的 40-60%
- 渐变填充（仅垂直或对角线）
- 分割背景（两种颜色，对角线或垂直）
- 边到边色带
- 负空间作为设计元素

### 布局技巧
**对于带有图表或表格的幻灯片：**
- **两列布局（首选）**：使用跨越整个宽度的标题，然后使用下面两列 - 一列中的文本/项目符号，另一列中的特色内容。这提供了更好的平衡并使图表/表格更具可读性。使用列宽不等的 Flexbox（例如，40%/60% 分割）来优化每种内容类型的空间。
- **完整幻灯片布局**：让特色内容（图表/表格）占据整个幻灯片，以获得最大的影响力和可读性
- **切勿垂直堆叠**：不要将图表/表格放在单列中的文本下方 - 这会导致可读性差和布局问题

### 工作流程
1. **强制 - 阅读整个文件**：从头到尾完整阅读 [`html2pptx.md`](html2pptx.md)。 **阅读此文件时切勿设置任何范围限制。** 在继续创建演示文稿之前，请阅读完整的文件内容以了解详细语法、关键格式规则和最佳实践。
2. 为每张幻灯片创建一个具有适当尺寸的 HTML 文件（例如，16:9 为 720pt × 405pt）
   - 对所有文本内容使用 `<p>`、`<h1>`-`<h6>`、`<ul>`、`<ol>`
   - 对将添加图表/表格的区域使用 `class="placeholder"`（以灰色背景渲染以提高可见性）
   - **关键**：首先使用 Sharp 将渐变和图标栅格化为 PNG 图像，然后在 HTML 中引用
- **布局**：对于带有图表/表格/图像的幻灯片，使用全幻灯片布局或两列布局以获得更好的可读性
3. 使用 [`html2pptx.js`](scripts/html2pptx.js) 库创建并运行 JavaScript 文件，将 HTML 幻灯片转换为 PowerPoint 并保存演示文稿
   - 使用`html2pptx()`函数处理每个HTML文件
   - 使用 PptxGenJS API 将图表和表格添加到占位符区域
   - 使用 `pptx.writeFile()` 保存演示文稿
4. **视觉验证**：生成缩略图并检查布局问题
   - 创建缩略图网格：`python scripts/thumbnail.py output.pptx workspace/thumbnails --cols 4`
   - 阅读并仔细检查缩略图：
     - **文本截断**：文本被标题栏、形状或幻灯片边缘截断
     - **文本重叠**：文本与其他文本或形状重叠
     - **定位问题**：内容距离幻灯片边界或其他元素太近
     - **对比度问题**：文本和背景之间的对比度不足
   - 如果发现问题，调整 HTML 边距/间距/颜色并重新生成演示文稿
   - 重复直到所有幻灯片视觉正确

## 编辑现有的 PowerPoint 演示文稿

要编辑现有 PowerPoint 演示文稿中的幻灯片，请使用原始 Office Open XML (OOXML) 格式。这涉及解压 .pptx 文件、编辑 XML 内容以及重新打包。

### 工作流程
1. **强制 - 读取整个文件**：从头到尾完整读取 [`ooxml.md`](ooxml.md)（约 500 行）。  **阅读此文件时切勿设置任何范围限制。** 在进行任何演示文稿编辑之前，请阅读完整的文件内容，以获取有关 OOXML 结构和编辑工作流程的详细指导。
2. 解压演示文稿：`python ooxml/scripts/unpack.py <office_file> <output_dir>`
3. 编辑XML文件（主要是`ppt/slides/slide{N}.xml`及相关文件）
4. **关键**：每次编辑后立即验证并修复任何验证错误，然后再继续：`python ooxml/scripts/validate.py <dir> --original <file>`
5. 打包最终演示文稿：`python ooxml/scripts/pack.py <input_directory> <office_file>`

## 使用模板创建新的 PowerPoint 演示文稿**

要创建遵循现有模板设计的演示文稿，请在替换占位符上下文之前复制并重新排列模板幻灯片。

### 工作流程
1. **提取模板文本并创建可视化缩略图网格**：
   * 提取文本：`python -m markitdown template.pptx > template-content.md`
   * 阅读`template-content.md`：阅读整个文件以了解模板演示的内容。 **读取此文件时切勿设置任何范围限制。**
   * 创建缩略图网格：`python scripts/thumbnail.py template.pptx`
   * 有关更多详细信息，请参阅[创建缩略图网格](#creating-thumbnail-grids) 部分

2. **分析模板并将库存保存到文件**：
   * **视觉分析**：查看缩略图网格以了解幻灯片布局、设计模式和视觉结构
   * 在`template-inventory.md`处创建并保存模板清单文件，其中包含：
     <<<代码块_1>>>
   * **使用缩略图网格**：参考视觉缩略图来识别：
     - 布局模式（标题幻灯片、内容布局、部分分隔符）
     - 图像占位符位置和计数
     - 幻灯片组之间的设计一致性
     - 视觉层次和结构
   * 需要此清单文件才能在下一步中选择适当的模板

3. **根据模板库存创建演示大纲**：
   * 查看第 2 步中的可用模板。
   * 为第一张幻灯片选择简介或标题模板。这应该是第一个模板。
   * 为其他幻灯片选择安全的、基于文本的布局。
   * **关键：将布局结构与实际内容相匹配**：
     - 单栏布局：用于统一叙述或单一主题
     - 两列布局：仅当恰好有 2 个不同的项目/概念时使用
     - 三列布局：仅当恰好有 3 个不同的项目/概念时使用
     - 图像 + 文本布局：仅在可插入实际图像时使用
     - 引用布局：仅用于人们的实际引用（带有归属），切勿用于强调
     - 切勿使用占位符多于可用内容的布局
     - 如果有 2 个项目，不要强迫它们采用 3 列布局
     - 如果有 4 个以上项目，请考虑分成多张幻灯片或使用列表格式
   * 在选择布局之前计算实际内容片段
   * 验证所选布局中的每个占位符都将填充有意义的内容
   * 选择一个代表每个内容部分的**最佳**布局的选项。
* 使用内容和利用可用设计的模板映射来保存 `outline.md`
   * 模板映射示例：
      <<<代码块_2>>>

4. **使用 `rearrange.py`** 复制、重新排序和删除幻灯片：
   * 使用 `scripts/rearrange.py` 脚本创建一个新的演示文稿，其中幻灯片按所需顺序排列：
     <<<代码块_3>>>
   * 该脚本可以处理重复幻灯片的复制、删除未使用的幻灯片以及自动重新排序
   * 幻灯片索引从 0 开始（第一张幻灯片是 0，第二张幻灯片是 1，等等）
   * 相同的幻灯片索引可以出现多次以复制该幻灯片

5. **使用 `inventory.py` 脚本提取所有文本**：
   * **运行库存提取**：
     <<<代码块_4>>>
   * **阅读 text-inventory.json**：阅读整个 text-inventory.json 文件以了解所有形状及其属性。 **读取此文件时切勿设置任何范围限制。**

   * 库存JSON结构：
      <<<代码块_5>>>

   * 主要特点：
     - **幻灯片**：命名为“slide-0”、“slide-1”等。
     - **形状**：按视觉位置（从上到下、从左到右）排序，如“shape-0”、“shape-1”等。
     - **占位符类型**：TITLE、CENTER_TITLE、SUBTITLE、BODY、OBJECT 或 null
     - **默认字体大小**：`default_font_size`，以从布局占位符提取的点为单位（如果可用）
     - **幻灯片编号已过滤**：具有 SLIDE_NUMBER 占位符类型的形状会自动从库存中排除
     - **项目符号**：当 `bullet: true` 时，始终包含 `level`（即使为 0）
     - **间距**：`space_before`、`space_after` 和 `line_spacing`（仅在设置时包含）
     - **颜色**：`color` 用于 RGB（例如，“FF0000”），`theme_color` 用于主题颜色（例如，“DARK_1”）
     - **属性**：输出中仅包含非默认值

6. **生成替换文本并将数据保存到 JSON 文件**
   基于上一步的文本清单：
   - **关键**：首先验证库存中存在哪些形状 - 仅参考实际存在的形状
   - **验证**：replace.py 脚本将验证替换 JSON 中的所有形状是否存在于库存中
     - 如果引用了不存在的形状，则会出现错误，显示可用的形状
     - 如果引用不存在的幻灯片，则会出现错误，指示该幻灯片不存在
     - 在脚本退出之前立即显示所有验证错误
   - **重要**：replace.py 脚本在内部使用 inventory.py 来识别所有文本形状
   - **自动清除**：库存中的所有文本形状都将被清除，除非您为它们提供“段落”
   - 将“段落”字段添加到需要内容的形状（而不是“replacement_paragraphs”）
   - 替换 JSON 中没有“段落”的形状将自动清除其文本
   - 带项目符号的段落将自动左对齐。当 `"bullet": true` 时，不要设置 `alignment` 属性
   - 为占位符文本生成适当的替换内容
   - 使用形状大小来确定适当的内容长度
   - **关键**：包括原始清单中的段落属性 - 不要只提供文本
   - **重要**：当项目符号：true 时，请勿在文本中包含项目符号（•、-、*） - 它们会自动添加
   - **基本格式规则**：
     - 标题/标题通常应具有 `"bold": true`
     - 列表项应具有 `"bullet": true, "level": 0` （当项目符号为 true 时，需要级别）
     - 保留任何对齐属性（例如，`"alignment": "CENTER"` 用于居中文本）
     - 包含与默认值不同的字体属性（例如，`"font_size": 14.0`、`"font_name": "Lora"`）
     - 颜色：使用 `"color": "FF0000"` 表示 RGB 或 `"theme_color": "DARK_1"` 表示主题颜色
     - 替换脚本需要**格式正确的段落**，而不仅仅是文本字符串
     - **重叠形状**：首选具有较大 default_font_size 或更合适的 placeholder_type 的形状
   - 将更新后的库存和替换品保存到 `replacement-text.json`
   - **警告**：不同的模板布局具有不同的形状数量 - 在创建替换件之前始终检查实际库存

   显示正确格式的示例段落字段：
   <<<代码块_6>>>

   **替换 JSON 中未列出的形状将自动清除**：
   ```json
   {
     "slide-0": {
       "shape-0": {
         "paragraphs": [...] // This shape gets new text
       }
       // shape-1 and shape-2 from inventory will be cleared automatically
     }
   }
   ```
**演示文稿的常见格式模式**：
   - 标题幻灯片：粗体文本，有时居中
   - 幻灯片中的节标题：粗体文本
   - 项目符号列表：每个项目需要`"bullet": true, "level": 0`
   - 正文：通常不需要特殊属性
   - 引号：可能有特殊的对齐方式或字体属性

7. **使用 `replace.py` 脚本应用替换**
   ```bash
   python scripts/replace.py working.pptx replacement-text.json output.pptx
   ```

   该脚本将：
   - 首先使用 inventory.py 中的函数提取所有文本形状的清单
   - 验证替换 JSON 中的所有形状是否存在于库存中
   - 库存中识别的所有形状的清晰文本
   - 仅将新文本应用于替换 JSON 中定义的“段落”形状
   - 通过应用 JSON 中的段落属性来保留格式
   - 自动处理项目符号、对齐方式、字体属性和颜色
   - 保存更新的演示文稿

   验证错误示例：
   ```
   ERROR: Invalid shapes in replacement JSON:
     - Shape 'shape-99' not found on 'slide-0'. Available shapes: shape-0, shape-1, shape-4
     - Slide 'slide-999' not found in inventory
   ```

   ```
   ERROR: Replacement text made overflow worse in these shapes:
     - slide-0/shape-2: overflow worsened by 1.25" (was 0.00", now 1.25")
   ```

## 创建缩略图网格

要创建 PowerPoint 幻灯片的可视化缩略图网格以供快速分析和参考：

```bash
python scripts/thumbnail.py template.pptx [output_prefix]
```

**特点**：
- 创建：`thumbnails.jpg`（或 `thumbnails-1.jpg`、`thumbnails-2.jpg` 等，用于大型牌组）
- 默认：5 列，每个网格最多 30 张幻灯片 (5×6)
- 自定义前缀：`python scripts/thumbnail.py template.pptx my-grid`
  - 注意：如果您希望在特定目录中输出，则输出前缀应包含路径（例如，`workspace/my-grid`）
- 调整列：`--cols 4`（范围：3-6，影响每个网格的幻灯片）
- 网格限制：3 列 = 12 张幻灯片/网格、4 列 = 20、5 列 = 30、6 列 = 42
- 幻灯片为零索引（幻灯片 0、幻灯片 1 等）

**用例**：
- 模板分析：快速了解幻灯片布局和设计模式
- 内容审查：整个演示文稿的视觉概述
- 导航参考：通过视觉外观查找特定幻灯片
- 质量检查：验证所有幻灯片的格式是否正确

**示例**：
```bash
# Basic usage
python scripts/thumbnail.py presentation.pptx

# Combine options: custom name, columns
python scripts/thumbnail.py template.pptx analysis --cols 4
```

## 将幻灯片转换为图像

要直观地分析 PowerPoint 幻灯片，请使用两步过程将其转换为图像：

1. **将 PPTX 转换为 PDF**:
   ```bash
   soffice --headless --convert-to pdf template.pptx
   ```

2. **将 PDF 页面转换为 JPEG 图像**：
   ```bash
   pdftoppm -jpeg -r 150 template.pdf slide
   ```
   这将创建诸如 `slide-1.jpg`、`slide-2.jpg` 等文件。

选项：
- `-r 150`：将分辨率设置为 150 DPI（调整质量/尺寸平衡）
- `-jpeg`：输出 JPEG 格式（如果愿意，请使用 `-png` 作为 PNG）
- `-f N`：要转换的第一页（例如，`-f 2` 从第 2 页开始）
- `-l N`：要转换的最后一页（例如，`-l 5` 在第 5 页停止）
- `slide`：输出文件的前缀

具体范围示例：
```bash
pdftoppm -jpeg -r 150 -f 2 -l 5 template.pdf slide  # Converts only pages 2-5
```

## 代码风格指南
**重要**：生成 PPTX 操作的代码时：
- 编写简洁的代码
- 避免冗长的变量名称和冗余操作
- 避免不必要的打印语句

## 依赖关系

所需的依赖项（应该已安装）：

- **markitdown**：`uv pip install "markitdown[pptx]"`（用于从演示文稿中提取文本）
- **pptxgenjs**：`npm install -g pptxgenjs`（用于通过 html2pptx 创建演示文稿）
- **playwright**：`npm install -g playwright`（用于 html2pptx 中的 HTML 渲染）
- **react-icons**：`npm install -g react-icons react react-dom`（对于图标）
- **锐利**：`npm install -g sharp`（用于 SVG 光栅化和图像处理）
- **LibreOffice**：`sudo apt-get install libreoffice`（用于 PDF 转换）
- **Poppler**：`sudo apt-get install poppler-utils`（用于 pdftoppm 将 PDF 转换为图像）
- **defusedxml**：`uv pip install defusedxml`（用于安全 XML 解析）