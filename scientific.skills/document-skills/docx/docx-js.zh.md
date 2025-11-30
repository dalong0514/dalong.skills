<!-- 此文件由机器翻译自 docx-js.md -->

# DOCX 库教程

使用 JavaScript/TypeScript 生成 .docx 文件。

**重要提示：开始之前请阅读整个文档。** 全文涵盖了关键格式规则和常见陷阱 - 跳过部分可能会导致文件损坏或呈现问题。

## 设置
假设 docx 已全局安装
如果未安装：`npm install -g docx`

```javascript
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, ImageRun, Media, 
        Header, Footer, AlignmentType, PageOrientation, LevelFormat, ExternalHyperlink, 
        InternalHyperlink, TableOfContents, HeadingLevel, BorderStyle, WidthType, TabStopType, 
        TabStopPosition, UnderlineType, ShadingType, VerticalAlign, SymbolRun, PageNumber,
        FootnoteReferenceRun, Footnote, PageBreak } = require('docx');

// Create & Save
const doc = new Document({ sections: [{ children: [/* content */] }] });
Packer.toBuffer(doc).then(buffer => fs.writeFileSync("doc.docx", buffer)); // Node.js
Packer.toBlob(doc).then(blob => { /* download logic */ }); // Browser
```

## 文本和格式
<<<代码块_1>>>

## 风格和专业格式

<<<代码块_2>>>

**专业字体组合：**
- **Arial（标题）+ Arial（正文）** - 最普遍支持、干净且专业
- **Times New Roman（标题）+ Arial（正文）** - 经典衬线标题与现代无衬线正文
- **Georgia（标题）+ Verdana（正文）** - 针对屏幕阅读进行了优化，对比度优雅

**关键造型原则：**
- **覆盖内置样式**：使用“Heading1”、“Heading2”、“Heading3”等精确ID来覆盖Word的内置标题样式
- **HeadingLevel 常量**：`HeadingLevel.HEADING_1` 使用“Heading1”样式，`HeadingLevel.HEADING_2` 使用“Heading2”样式等。
- **包括outlineLevel**：为H1设置`outlineLevel: 0`，为H2设置`outlineLevel: 1`等，以确保TOC正常工作
- **使用自定义样式**而不是内联格式以保持一致性
- **使用 `styles.default.document.run.font` 设置默认字体** - Arial 得到普遍支持
- **建立具有不同字体大小的视觉层次结构**（标题>标题>正文）
- **添加适当的间距**，使用 `before` 和 `after` 段落间距
- **谨慎使用颜色**：标题和标题默认为黑色 (000000) 和灰色阴影（标题 1、标题 2 等）
- **设置一致的边距**（1440 = 1 英寸为标准）


## 列表（始终使用正确的列表 - 切勿使用 UNICODE 项目符号）
<<<代码块_3>>>

## 表格
<<<代码块_4>>>

**重要：表格宽度和边框**
- 在每个单元格上同时使用 `columnWidths: [width1, width2, ...]` 数组和 `width: { size: X, type: WidthType.DXA }`
- DXA 值（二十分之一点）：1440 = 1 英寸，字母可用宽度 = 9360 DXA（边距为 1 英寸）
- 将边框应用于各个 `TableCell` 元素，而不是 `Table` 本身

**预先计算的列宽度（具有 1 英寸边距的字母大小 = 9360 DXA 总计）：**
- **2 列：** `columnWidths: [4680, 4680]`（等宽）
- **3 列：** `columnWidths: [3120, 3120, 3120]`（等宽）

## 链接和导航
<<<代码块_5>>>

## 图片与媒体
<<<代码块_6>>>

## 分页符
```javascript
// Manual page break
new Paragraph({ children: [new PageBreak()] }),

// Page break before paragraph
new Paragraph({
  pageBreakBefore: true,
  children: [new TextRun("This starts on a new page")]
})

// ⚠️ CRITICAL: NEVER use PageBreak standalone - it will create invalid XML that Word cannot open
// ❌ WRONG: new PageBreak() 
// ✅ CORRECT: new Paragraph({ children: [new PageBreak()] })
```

## 页眉/页脚和页面设置
```javascript
const doc = new Document({
  sections: [{
    properties: {
      page: {
        margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 }, // 1440 = 1 inch
        size: { orientation: PageOrientation.LANDSCAPE },
        pageNumbers: { start: 1, formatType: "decimal" } // "upperRoman", "lowerRoman", "upperLetter", "lowerLetter"
      }
    },
    headers: {
      default: new Header({ children: [new Paragraph({ 
        alignment: AlignmentType.RIGHT,
        children: [new TextRun("Header Text")]
      })] })
    },
    footers: {
      default: new Footer({ children: [new Paragraph({ 
        alignment: AlignmentType.CENTER,
        children: [new TextRun("Page "), new TextRun({ children: [PageNumber.CURRENT] }), new TextRun(" of "), new TextRun({ children: [PageNumber.TOTAL_PAGES] })]
      })] })
    },
    children: [/* content */]
  }]
});
```

## 标签
```javascript
new Paragraph({
  tabStops: [
    { type: TabStopType.LEFT, position: TabStopPosition.MAX / 4 },
    { type: TabStopType.CENTER, position: TabStopPosition.MAX / 2 },
    { type: TabStopType.RIGHT, position: TabStopPosition.MAX * 3 / 4 }
  ],
  children: [new TextRun("Left\tCenter\tRight")]
})
```

## 常量和快速参考
- **下划线：** `SINGLE`、`DOUBLE`、`WAVY`、`DASH`
- **边框：** `SINGLE`、`DOUBLE`、`DASHED`、`DOTTED`  
- **编号：** `DECIMAL` (1,2,3), `UPPER_ROMAN` (I,II,III), `LOWER_LETTER` (a,b,c)
- **制表符：** `LEFT`、`CENTER`、`RIGHT`、`DECIMAL`
- **符号：** `"2022"` (•)、`"00A9"` (©)、`"00AE"` (®)、`"2122"` (™)、`"00B0"` (°)、`"F070"` (✓)、`"F0FC"` (✗)

## 关键问题和常见错误
- **关键：PageBreak 必须始终位于段落内** - 独立的 PageBreak 会创建 Word 无法打开的无效 XML
- **始终使用 ShadingType.CLEAR 进行表格单元格着色** - 切勿使用 ShadingType.SOLID（导致黑色背景）。
- DXA 测量（1440 = 1 英寸）|每个表格单元格需要≥1 段 | TOC 仅需要 HeadingLevel 样式
- **始终使用自定义样式**和 Arial 字体，以获得专业的外观和适当的视觉层次结构
- **始终使用 `styles.default.document.run.font` 设置默认字体** - 推荐 Arial
- **始终使用表的 columnWidths 数组** + 单个单元格宽度以实现兼容性
- **切勿对项目符号使用 unicode 符号** - 始终使用正确的编号配置和 `LevelFormat.BULLET` 常量（不是字符串“bullet”）
- **切勿在任何地方使用 \n 进行换行** - 始终为每一行使用单独的段落元素
- **始终在段落子项中使用 TextRun 对象** - 切勿直接在段落上使用文本属性
- **对于图像至关重要**：ImageRun 需要 `type` 参数 - 始终指定“png”、“jpg”、“jpeg”、“gif”、“bmp”或“svg”
- **对于项目符号至关重要**：必须使用 `LevelFormat.BULLET` 常量，而不是字符串“bullet”，并包含 `text: "•"` 作为项目符号字符
- **对于编号至关重要**：每个编号参考都会创建一个独立列表。相同参考 = 继续编号（1、2、3，然后 4、5、6）。不同的参考 = 从 1 重新开始（1,2,3，然后 1,2,3）。为每个单独的编号部分使用唯一的参考名称！
- **目录的关键**：使用目录时，标题必须仅使用 HeadingLevel - 不要向标题段落添加自定义样式，否则目录将中断
- **表格**：设置`columnWidths`数组+单个单元格宽度，将边框应用于单元格而不是表格
- **在表格级别设置表格边距**以保持单元格填充一致（避免每个单元格重复）