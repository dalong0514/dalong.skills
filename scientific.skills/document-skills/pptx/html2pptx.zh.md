<!-- 此文件由机器翻译自 html2pptx.md -->

# HTML 到 PowerPoint 指南

使用 `html2pptx.js` 库将 HTML 幻灯片转换为具有精确定位的 PowerPoint 演示文稿。

## 目录

1. [创建 HTML 幻灯片](#creating-html-slides)
2. [使用 html2pptx 库](#using-the-html2pptx-library)
3. [使用 PptxGenJS](#using-pptxgenjs)

---

## 创建 HTML 幻灯片

每张 HTML 幻灯片都必须包含正确的正文尺寸：

### 布局尺寸

- **16:9**（默认）：`width: 720pt; height: 405pt`
- **4:3**: `width: 720pt; height: 540pt`
- **16:10**：`width: 720pt; height: 450pt`

### 支持的元素

- `<p>`、`<h1>`-`<h6>` - 带样式的文本
- `<ul>`、`<ol>` - 列表（切勿使用手动项目符号•、-、*）
- `<b>`、`<strong>` - 粗体文本（内联格式）
- `<i>`、`<em>` - 斜体文本（内联格式）
- `<u>` - 带下划线的文本（内联格式）
- `<span>` - 使用 CSS 样式进行内联格式化（粗体、斜体、下划线、颜色）
- `<br>` - 换行符
- `<div>` 与 bg/border - 变成形状
- `<img>` - 图像
- `class="placeholder"` - 为图表保留空间（返回`{ id, x, y, w, h }`）

### 关键文本规则

**所有文本必须位于 `<p>`、`<h1>`-`<h6>`、`<ul>` 或 `<ol>` 标签内：**
- ✅ 正确：`<div><p>Text here</p></div>`
- ❌ 错误：`<div>Text here</div>` - **文本不会出现在 PowerPoint 中**
- ❌ 错误：`<span>Text</span>` - **文本不会出现在 PowerPoint 中**
- `<div>` 或 `<span>` 中没有文本标记的文本将被静默忽略

**切勿使用手动项目符号（•、-、* 等）** - 使用 `<ul>` 或 `<ol>` 列表代替

**仅使用普遍可用的网络安全字体：**
- ✅ 网络安全字体：`Arial`、`Helvetica`、`Times New Roman`、`Georgia`、`Courier New`、`Verdana`、`Tahoma`、`Trebuchet MS`、 `Impact`、`Comic Sans MS`
- ❌ 错误：`'Segoe UI'`、`'SF Pro'`、`'Roboto'`、自定义字体 - **可能会导致渲染问题**

### 造型

- 在正文上使用 `display: flex` 来防止边距崩溃破坏溢出验证
- 使用 `margin` 作为间距（尺寸中包含填充）
- 内联格式：使用 `<b>`、`<i>`、`<u>` 标签或 `<span>` 与 CSS 样式
  - `<span>` 支持：`font-weight: bold`、`font-style: italic`、`text-decoration: underline`、`color: #rrggbb`
  - `<span>` 不支持：`margin`、`padding`（在 PowerPoint 文本运行中不受支持）
  - 示例：`<span style="font-weight: bold; color: #667eea;">Bold blue text</span>`
- Flexbox 工作 - 根据渲染布局计算位置
- 在 CSS 中使用带有 `#` 前缀的十六进制颜色
- **文本对齐**：如果文本长度稍有偏差，请在需要时使用 CSS `text-align`（`center`、`right` 等）作为 PptxGenJS 文本格式的提示

### 形状样式（仅限 DIV 元素）

**重要提示：背景、边框和阴影仅适用于 `<div>` 元素，不适用于文本元素（`<p>`、`<h1>`-`<h6>`、`<ul>`、`<ol>`）**

- **背景**：仅 `<div>` 元素上的 CSS `background` 或 `background-color`
  - 示例：`<div style="background: #f0f0f0;">` - 创建带背景的形状
- **边框**：`<div>` 元素上的 CSS `border` 转换为 PowerPoint 形状边框
  - 支持统一边框：`border: 2px solid #333333`
  - 支持部分边框：`border-left`、`border-right`、`border-top`、`border-bottom`（呈现为线条形状）
  - 示例：`<div style="border-left: 8pt solid #E76F51;">`
- **边框半径**：圆角 `<div>` 元素上的 CSS `border-radius`
  - `border-radius: 50%` 或更高版本创建圆形
  - 相对于形状较小尺寸计算的百分比 <50%
  - 支持 px 和 pt 单位（例如，`border-radius: 8pt;`、`border-radius: 12px;`）
  - 示例：`<div style="border-radius: 25%;">` 在 100x200px 框上 = 100px 的 25% = 25px 半径
- **盒子阴影**：`<div>` 元素上的 CSS `box-shadow` 转换为 PowerPoint 阴影
  - 仅支持外部阴影（忽略内部阴影以防止损坏）
- 示例：`<div style="box-shadow: 2px 2px 8px rgba(0, 0, 0, 0.3);">`
  - 注意：PowerPoint 不支持插入/内部阴影，将被跳过

### 图标和渐变

- **关键：切勿使用 CSS 渐变（`linear-gradient`、`radial-gradient`）** - 它们不会转换为 PowerPoint
- **始终首先使用 Sharp 创建渐变/图标 PNG，然后在 HTML 中引用**
- 对于渐变：将 SVG 栅格化为 PNG 背景图像
- 对于图标：将react-icons SVG 栅格化为PNG 图像
- 所有视觉效果必须在 HTML 渲染之前预先渲染为光栅图像

**使用 Sharp 光栅化图标：**

```javascript
const React = require('react');
const ReactDOMServer = require('react-dom/server');
const sharp = require('sharp');
const { FaHome } = require('react-icons/fa');

async function rasterizeIconPng(IconComponent, color, size = "256", filename) {
  const svgString = ReactDOMServer.renderToStaticMarkup(
    React.createElement(IconComponent, { color: `#${color}`, size: size })
  );

  // Convert SVG to PNG using Sharp
  await sharp(Buffer.from(svgString))
    .png()
    .toFile(filename);

  return filename;
}

// Usage: Rasterize icon before using in HTML
const iconPath = await rasterizeIconPng(FaHome, "4472c4", "256", "home-icon.png");
// Then reference in HTML: <img src="home-icon.png" style="width: 40pt; height: 40pt;">
```

**使用锐利光栅化渐变：**

<<<代码块_1>>>

### 示例

<<<代码块_2>>>

## 使用 html2pptx 库

### 依赖关系

这些库已全局安装并可供使用：
- `pptxgenjs`
- `playwright`
- `sharp`

### 基本用法

<<<代码块_3>>>

### API 参考

#### 函数签名
<<<代码块_4>>>

#### 参数
- `htmlFile`（字符串）：HTML 文件的路径（绝对或相对）
- `pres` (pptxgen)：已设置布局的 PptxGenJS 演示文稿实例
- `options`（对象，可选）：
  - `tmpDir`（字符串）：生成文件的临时目录（默认值：`process.env.TMPDIR || '/tmp'`）
  - `slide`（对象）：要重用的现有幻灯片（默认值：创建新幻灯片）

#### 返回
<<<代码块_5>>>

### 验证

该库在抛出之前自动验证并收集所有错误：

1. **HTML 尺寸必须匹配演示文稿布局** - 报告尺寸不匹配
2. **内容不得溢出主体** - 以精确的测量值报告溢出
3. **CSS 渐变** - 报告不支持的渐变使用
4. **文本元素样式** - 报告文本元素上的背景/边框/阴影（仅允许在 div 上）

**所有验证错误都会收集并在一条错误消息中一起报告**，使您可以一次修复所有问题，而不是一次解决一个问题。

### 使用占位符

<<<代码块_6>>>

### 完整示例

```javascript
const pptxgen = require('pptxgenjs');
const html2pptx = require('./html2pptx');

async function createPresentation() {
    const pptx = new pptxgen();
    pptx.layout = 'LAYOUT_16x9';
    pptx.author = 'Your Name';
    pptx.title = 'My Presentation';

    // Slide 1: Title
    const { slide: slide1 } = await html2pptx('slides/title.html', pptx);

    // Slide 2: Content with chart
    const { slide: slide2, placeholders } = await html2pptx('slides/data.html', pptx);

    const chartData = [{
        name: 'Sales',
        labels: ['Q1', 'Q2', 'Q3', 'Q4'],
        values: [4500, 5500, 6200, 7100]
    }];

    slide2.addChart(pptx.charts.BAR, chartData, {
        ...placeholders[0],
        showTitle: true,
        title: 'Quarterly Sales',
        showCatAxisTitle: true,
        catAxisTitle: 'Quarter',
        showValAxisTitle: true,
        valAxisTitle: 'Sales ($000s)'
    });

    // Save
    await pptx.writeFile({ fileName: 'presentation.pptx' });
    console.log('Presentation created successfully!');
}

createPresentation().catch(console.error);
```

## 使用 PptxGenJS

使用 `html2pptx` 将 HTML 转换为幻灯片后，您将使用 PptxGenJS 添加动态内容，例如图表、图像和其他元素。

### ⚠️ 关键规则

#### 颜色
- **切勿在 PptxGenJS 中使用带有十六进制颜色的 `#` 前缀** - 导致文件损坏
- ✅ 正确：`color: "FF0000"`、`fill: { color: "0066CC" }`
- ❌错误：`color: "#FF0000"`（破坏文档）

### 添加图像

始终根据实际图像尺寸计算纵横比：

```javascript
// Get image dimensions: identify image.png | grep -o '[0-9]* x [0-9]*'
const imgWidth = 1860, imgHeight = 1519;  // From actual file
const aspectRatio = imgWidth / imgHeight;

const h = 3;  // Max height
const w = h * aspectRatio;
const x = (10 - w) / 2;  // Center on 16:9 slide

slide.addImage({ path: "chart.png", x, y: 1.5, w, h });
```

### 添加文本

```javascript
// Rich text with formatting
slide.addText([
    { text: "Bold ", options: { bold: true } },
    { text: "Italic ", options: { italic: true } },
    { text: "Normal" }
], {
    x: 1, y: 2, w: 8, h: 1
});
```

### 添加形状

```javascript
// Rectangle
slide.addShape(pptx.shapes.RECTANGLE, {
    x: 1, y: 1, w: 3, h: 2,
    fill: { color: "4472C4" },
    line: { color: "000000", width: 2 }
});

// Circle
slide.addShape(pptx.shapes.OVAL, {
    x: 5, y: 1, w: 2, h: 2,
    fill: { color: "ED7D31" }
});

// Rounded rectangle
slide.addShape(pptx.shapes.ROUNDED_RECTANGLE, {
    x: 1, y: 4, w: 3, h: 1.5,
    fill: { color: "70AD47" },
    rectRadius: 0.2
});
```

### 添加图表

**大多数图表必需：** 使用 `catAxisTitle`（类别）和 `valAxisTitle`（值）的轴标签。

**图表数据格式：**
- 使用**带有所有标签的单一系列**来制作简单的条形图/折线图
- 每个系列都会创建一个单独的图例条目
- 标签数组定义 X 轴值

**时间序列数据 - 选择正确的粒度：**
- **< 30 天**：使用每日分组（例如，“10-01”、“10-02”）- 避免创建单点图表的每月聚合
- **30-365 天**：使用每月分组（例如“2024-01”、“2024-02”）
- **> 365 天**：使用年度分组（例如“2023”、“2024”）
- **验证**：只有 1 个数据点的图表可能表明该时间段的聚合不正确

```javascript
const { slide, placeholders } = await html2pptx('slide.html', pptx);

// CORRECT: Single series with all labels
slide.addChart(pptx.charts.BAR, [{
    name: "Sales 2024",
    labels: ["Q1", "Q2", "Q3", "Q4"],
    values: [4500, 5500, 6200, 7100]
}], {
    ...placeholders[0],  // Use placeholder position
    barDir: 'col',       // 'col' = vertical bars, 'bar' = horizontal
    showTitle: true,
    title: 'Quarterly Sales',
    showLegend: false,   // No legend needed for single series
    // Required axis labels
    showCatAxisTitle: true,
    catAxisTitle: 'Quarter',
    showValAxisTitle: true,
    valAxisTitle: 'Sales ($000s)',
    // Optional: Control scaling (adjust min based on data range for better visualization)
    valAxisMaxVal: 8000,
    valAxisMinVal: 0,  // Use 0 for counts/amounts; for clustered data (e.g., 4500-7100), consider starting closer to min value
    valAxisMajorUnit: 2000,  // Control y-axis label spacing to prevent crowding
    catAxisLabelRotate: 45,  // Rotate labels if crowded
    dataLabelPosition: 'outEnd',
    dataLabelColor: '000000',
    // Use single color for single-series charts
    chartColors: ["4472C4"]  // All bars same color
});
```

#### 散点图

**重要**：散点图数据格式不寻常 - 第一个系列包含 X 轴值，后续系列包含 Y 值：

```javascript
// Prepare data
const data1 = [{ x: 10, y: 20 }, { x: 15, y: 25 }, { x: 20, y: 30 }];
const data2 = [{ x: 12, y: 18 }, { x: 18, y: 22 }];

const allXValues = [...data1.map(d => d.x), ...data2.map(d => d.x)];

slide.addChart(pptx.charts.SCATTER, [
    { name: 'X-Axis', values: allXValues },  // First series = X values
    { name: 'Series 1', values: data1.map(d => d.y) },  // Y values only
    { name: 'Series 2', values: data2.map(d => d.y) }   // Y values only
], {
    x: 1, y: 1, w: 8, h: 4,
    lineSize: 0,  // 0 = no connecting lines
    lineDataSymbol: 'circle',
    lineDataSymbolSize: 6,
    showCatAxisTitle: true,
    catAxisTitle: 'X Axis',
    showValAxisTitle: true,
    valAxisTitle: 'Y Axis',
    chartColors: ["4472C4", "ED7D31"]
});
```

#### 折线图

```javascript
slide.addChart(pptx.charts.LINE, [{
    name: "Temperature",
    labels: ["Jan", "Feb", "Mar", "Apr"],
    values: [32, 35, 42, 55]
}], {
    x: 1, y: 1, w: 8, h: 4,
    lineSize: 4,
    lineSmooth: true,
    // Required axis labels
    showCatAxisTitle: true,
    catAxisTitle: 'Month',
    showValAxisTitle: true,
    valAxisTitle: 'Temperature (°F)',
    // Optional: Y-axis range (set min based on data range for better visualization)
    valAxisMinVal: 0,     // For ranges starting at 0 (counts, percentages, etc.)
    valAxisMaxVal: 60,
    valAxisMajorUnit: 20,  // Control y-axis label spacing to prevent crowding (e.g., 10, 20, 25)
    // valAxisMinVal: 30,  // PREFERRED: For data clustered in a range (e.g., 32-55 or ratings 3-5), start axis closer to min value to show variation
    // Optional: Chart colors
    chartColors: ["4472C4", "ED7D31", "A5A5A5"]
});
```

#### 饼图（无需轴标签）

**关键**：饼图需要一个**单个数据系列**，其中包含 `labels` 数组中的所有类别以及 `values` 数组中的相应值。

```javascript
slide.addChart(pptx.charts.PIE, [{
    name: "Market Share",
    labels: ["Product A", "Product B", "Other"],  // All categories in one array
    values: [35, 45, 20]  // All values in one array
}], {
    x: 2, y: 1, w: 6, h: 4,
    showPercent: true,
    showLegend: true,
    legendPos: 'r',  // right
    chartColors: ["4472C4", "ED7D31", "A5A5A5"]
});
```

#### 多个数据系列

```javascript
slide.addChart(pptx.charts.LINE, [
    {
        name: "Product A",
        labels: ["Q1", "Q2", "Q3", "Q4"],
        values: [10, 20, 30, 40]
    },
    {
        name: "Product B",
        labels: ["Q1", "Q2", "Q3", "Q4"],
        values: [15, 25, 20, 35]
    }
], {
    x: 1, y: 1, w: 8, h: 4,
    showCatAxisTitle: true,
    catAxisTitle: 'Quarter',
    showValAxisTitle: true,
    valAxisTitle: 'Revenue ($M)'
});
```

### 图表颜色

**关键**：使用十六进制颜色**不带** `#` 前缀 - 包括 `#` 会导致文件损坏。

**将图表颜色与您选择的设计调色板对齐**，确保数据可视化有足够的对比度和独特性。调整颜色：
- 相邻系列之间的强烈对比
- 幻灯片背景的可读性
- 辅助功能（避免仅红绿组合）

```javascript
// Example: Ocean palette-inspired chart colors (adjusted for contrast)
const chartColors = ["16A085", "FF6B9D", "2C3E50", "F39C12", "9B59B6"];

// Single-series chart: Use one color for all bars/points
slide.addChart(pptx.charts.BAR, [{
    name: "Sales",
    labels: ["Q1", "Q2", "Q3", "Q4"],
    values: [4500, 5500, 6200, 7100]
}], {
    ...placeholders[0],
    chartColors: ["16A085"],  // All bars same color
    showLegend: false
});

// Multi-series chart: Each series gets a different color
slide.addChart(pptx.charts.LINE, [
    { name: "Product A", labels: ["Q1", "Q2", "Q3"], values: [10, 20, 30] },
    { name: "Product B", labels: ["Q1", "Q2", "Q3"], values: [15, 25, 20] }
], {
    ...placeholders[0],
    chartColors: ["16A085", "FF6B9D"]  // One color per series
});
```

### 添加表

可以使用基本或高级格式添加表格：

#### 基本表

```javascript
slide.addTable([
    ["Header 1", "Header 2", "Header 3"],
    ["Row 1, Col 1", "Row 1, Col 2", "Row 1, Col 3"],
    ["Row 2, Col 1", "Row 2, Col 2", "Row 2, Col 3"]
], {
    x: 0.5,
    y: 1,
    w: 9,
    h: 3,
    border: { pt: 1, color: "999999" },
    fill: { color: "F1F1F1" }
});
```
#### 具有自定义格式的表格

```javascript
const tableData = [
    // Header row with custom styling
    [
        { text: "Product", options: { fill: { color: "4472C4" }, color: "FFFFFF", bold: true } },
        { text: "Revenue", options: { fill: { color: "4472C4" }, color: "FFFFFF", bold: true } },
        { text: "Growth", options: { fill: { color: "4472C4" }, color: "FFFFFF", bold: true } }
    ],
    // Data rows
    ["Product A", "$50M", "+15%"],
    ["Product B", "$35M", "+22%"],
    ["Product C", "$28M", "+8%"]
];

slide.addTable(tableData, {
    x: 1,
    y: 1.5,
    w: 8,
    h: 3,
    colW: [3, 2.5, 2.5],  // Column widths
    rowH: [0.5, 0.6, 0.6, 0.6],  // Row heights
    border: { pt: 1, color: "CCCCCC" },
    align: "center",
    valign: "middle",
    fontSize: 14
});
```

#### 具有合并单元格的表格

```javascript
const mergedTableData = [
    [
        { text: "Q1 Results", options: { colspan: 3, fill: { color: "4472C4" }, color: "FFFFFF", bold: true } }
    ],
    ["Product", "Sales", "Market Share"],
    ["Product A", "$25M", "35%"],
    ["Product B", "$18M", "25%"]
];

slide.addTable(mergedTableData, {
    x: 1,
    y: 1,
    w: 8,
    h: 2.5,
    colW: [3, 2.5, 2.5],
    border: { pt: 1, color: "DDDDDD" }
});
```

### 表选项

常用表选项：
- `x, y, w, h` - 位置和大小
- `colW` - 列宽数组（以英寸为单位）
- `rowH` - 行高数组（以英寸为单位）
- `border` - 边框样式：`{ pt: 1, color: "999999" }`
- `fill` - 背景颜色（无 # 前缀）
- `align` - 文本对齐方式：“左”、“中”、“右”
- `valign` - 垂直对齐：“顶部”、“中间”、“底部”
- `fontSize` - 文本大小
- `autoPage` - 如果内容溢出，自动创建新幻灯片