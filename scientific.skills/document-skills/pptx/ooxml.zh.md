<!-- 此文件由机器翻译自 ooxml.md -->

# PowerPoint 的 Office Open XML 技术参考

**重要提示：开始之前请阅读整个文档。** 全文涵盖了关键的 XML 架构规则和格式要求。不正确的实施可能会创建 PowerPoint 无法打开的无效 PPTX 文件。

## 技术指南

### 架构合规性
- **`<p:txBody>`** 中的元素排序：`<a:bodyPr>`、`<a:lstStyle>`、`<a:p>`
- **空白**：将 `xml:space='preserve'` 添加到带有前导/尾随空格的 `<a:t>` 元素
- **Unicode**：ASCII 内容中的转义字符：`"` 变为 `&#8220;`
- **图像**：添加到 `ppt/media/`，在幻灯片 XML 中引用，设置尺寸以适合幻灯片边界
- **关系**：更新每张幻灯片资源的 `ppt/slides/_rels/slideN.xml.rels`
- **脏属性**：将 `dirty="0"` 添加到 `<a:rPr>` 和 `<a:endParaRPr>` 元素以指示干净状态

## 演示结构

### 基本幻灯片结构
```xml
<!-- ppt/slides/slide1.xml -->
<p:sld>
  <p:cSld>
    <p:spTree>
      <p:nvGrpSpPr>...</p:nvGrpSpPr>
      <p:grpSpPr>...</p:grpSpPr>
      <!-- Shapes go here -->
    </p:spTree>
  </p:cSld>
</p:sld>
```

### 文本框/带文本的形状
<<<代码块_1>>>

### 文本格式
<<<代码块_2>>>

### 列表
<<<代码块_3>>>

### 形状
<<<代码块_4>>>

### 图片
<<<代码块_5>>>

### 表格
<<<代码块_6>>>

### 幻灯片布局

```xml
<!-- Title Slide Layout -->
<p:sp>
  <p:nvSpPr>
    <p:nvPr>
      <p:ph type="ctrTitle"/>
    </p:nvPr>
  </p:nvSpPr>
  <!-- Title content -->
</p:sp>

<p:sp>
  <p:nvSpPr>
    <p:nvPr>
      <p:ph type="subTitle" idx="1"/>
    </p:nvPr>
  </p:nvSpPr>
  <!-- Subtitle content -->
</p:sp>

<!-- Content Slide Layout -->
<p:sp>
  <p:nvSpPr>
    <p:nvPr>
      <p:ph type="title"/>
    </p:nvPr>
  </p:nvSpPr>
  <!-- Slide title -->
</p:sp>

<p:sp>
  <p:nvSpPr>
    <p:nvPr>
      <p:ph type="body" idx="1"/>
    </p:nvPr>
  </p:nvSpPr>
  <!-- Content body -->
</p:sp>
```

## 文件更新

添加内容时，更新这些文件：

**`ppt/_rels/presentation.xml.rels`：**
```xml
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slide" Target="slides/slide1.xml"/>
<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideMaster" Target="slideMasters/slideMaster1.xml"/>
```

**`ppt/slides/_rels/slide1.xml.rels`：**
```xml
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout" Target="../slideLayouts/slideLayout1.xml"/>
<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="../media/image1.png"/>
```

**`[Content_Types].xml`：**
```xml
<Default Extension="png" ContentType="image/png"/>
<Default Extension="jpg" ContentType="image/jpeg"/>
<Override PartName="/ppt/slides/slide1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slide+xml"/>
```

**`ppt/presentation.xml`：**
```xml
<p:sldIdLst>
  <p:sldId id="256" r:id="rId1"/>
  <p:sldId id="257" r:id="rId2"/>
</p:sldIdLst>
```

**`docProps/app.xml`:** 更新幻灯片计数和统计信息
```xml
<Slides>2</Slides>
<Paragraphs>10</Paragraphs>
<Words>50</Words>
```

## 幻灯片操作

### 添加新幻灯片
将幻灯片添加到演示文稿末尾时：

1. **创建幻灯片文件** (`ppt/slides/slideN.xml`)
2. **更新`[Content_Types].xml`**：为新幻灯片添加覆盖
3. **更新`ppt/_rels/presentation.xml.rels`**：为新幻灯片添加关系
4. **更新`ppt/presentation.xml`**：将幻灯片ID添加到`<p:sldIdLst>`
5. **创建幻灯片关系** (`ppt/slides/_rels/slideN.xml.rels`)（如果需要）
6. **更新 `docProps/app.xml`**：增加幻灯片计数并更新统计信息（如果存在）

### 复制幻灯片
1. 使用新名称复制源幻灯片 XML 文件
2. 将新幻灯片中的所有 ID 更新为唯一
3. 按照上面的“添加新幻灯片”步骤操作
4. **关键**：删除或更新 `_rels` 文件中的任何注释幻灯片引用
5.删除对未使用的媒体文件的引用

### 重新排序幻灯片
1. **更新 `ppt/presentation.xml`**：重新排序 `<p:sldIdLst>` 中的 `<p:sldId>` 元素
2. `<p:sldId>`元素的顺序决定幻灯片顺序
3. 保持幻灯片 ID 和关系 ID 不变

示例：
```xml
<!-- Original order -->
<p:sldIdLst>
  <p:sldId id="256" r:id="rId2"/>
  <p:sldId id="257" r:id="rId3"/>
  <p:sldId id="258" r:id="rId4"/>
</p:sldIdLst>

<!-- After moving slide 3 to position 2 -->
<p:sldIdLst>
  <p:sldId id="256" r:id="rId2"/>
  <p:sldId id="258" r:id="rId4"/>
  <p:sldId id="257" r:id="rId3"/>
</p:sldIdLst>
```

### 删除幻灯片
1. **从 `ppt/presentation.xml`** 中删除：删除 `<p:sldId>` 条目
2. **Remove from `ppt/_rels/presentation.xml.rels`**：删除关系
3. **从`[Content_Types].xml`**中删除：删除覆盖条目
4. **删除文件**：删除`ppt/slides/slideN.xml`和`ppt/slides/_rels/slideN.xml.rels`
5. **更新`docProps/app.xml`**：减少幻灯片计数并更新统计信息
6. **清理未使用的媒体**：从 `ppt/media/` 中删除孤立的图像

注意：不要对剩余的幻灯片重新编号 - 保留其原始 ID 和文件名。


## 要避免的常见错误

- **编码**：转义 ASCII 内容中的 unicode 字符：`"` 变为 `&#8220;`
- **图像**：添加到 `ppt/media/` 并更新关系文件
- **列表**：省略列表标题中的项目符号
- **ID**：对 UUID 使用有效的十六进制值
- **主题**：检查 `theme` 目录中的所有主题的颜色

## 基于模板的演示文稿的验证清单

### 打包前，请务必：
- **清理未使用的资源**：删除未引用的媒体、字体和注释目录
- **修复 Content_Types.xml**：声明包中存在的所有幻灯片、布局和主题
- **修复关系 ID**： 
   - 如果不使用嵌入字体，请删除字体嵌入引用
- **删除损坏的引用**：检查所有 `_rels` 文件中对已删除资源的引用

### 常见模板重复陷阱：
- 复制后多张幻灯片引用相同的笔记幻灯片
- 模板幻灯片中的图像/媒体引用不再存在
- 不包含字体时的字体嵌入参考
- 缺少布局12-25的slideLayout声明
- docProps 目录可能无法解压 - 这是可选的