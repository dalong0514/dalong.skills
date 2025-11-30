<!-- 此文件由机器翻译自 ooxml.md -->

# Office Open XML 技术参考

**重要提示：开始之前请阅读整个文档。**本文档涵盖：
- [技术指南](#technical-guidelines) - 架构合规性规则和验证要求
- [文档内容模式](#document-content-patterns) - 用于标题、列表、表格、格式等的 XML 模式。
- [文档库 (Python)](#document-library-python) - 通过自动基础设施设置进行 OOXML 操作的推荐方法
- [跟踪更改 (Redlined)](#tracked-changes-redlined) - 用于实现跟踪更改的 XML 模式

## 技术指南

### 架构合规性
- **`<w:pPr>`** 中的元素排序：`<w:pStyle>`、`<w:numPr>`、`<w:spacing>`、`<w:ind>`、`<w:jc>`
- **空白**：将 `xml:space='preserve'` 添加到带有前导/尾随空格的 `<w:t>` 元素
- **Unicode**：ASCII 内容中的转义字符：`"` 变为 `&#8220;`
  - **字符编码参考**：大引号 `""` 变为 `&#8220;&#8221;`，撇号 `'` 变为 `&#8217;`，破折号 `—` 变为 `&#8212;`
- **跟踪的更改**：在 `<w:r>` 元素外部使用 `<w:del>` 和 `<w:ins>` 标签以及 `w:author="Claude"` 元素
  - **关键**：`<w:ins>` 以 `</w:ins>` 关闭，`<w:del>` 以 `</w:del>` 关闭 - 切勿混合
  - **RSID 必须是 8 位十六进制**：使用诸如 `00AB1234` 之类的值（仅限 0-9、A-F 字符）
  - **trackRevisions 放置**：在 settings.xml 中的 `<w:proofState>` 之后添加 `<w:trackRevisions/>`
- **图像**：添加到`word/media/`，在`document.xml`中引用，设置尺寸以防止溢出

## 文档内容模式

### 基本结构
```xml
<w:p>
  <w:r><w:t>Text content</w:t></w:r>
</w:p>
```

### 标题和样式
<<<代码块_1>>>

### 文本格式
<<<代码块_2>>>

### 列表
<<<代码块_3>>>

### 表格
<<<代码块_4>>>

### 布局
<<<代码块_5>>>

## 文件更新

添加内容时，更新这些文件：

**`word/_rels/document.xml.rels`：**
<<<代码块_6>>>

**`[Content_Types].xml`：**
```xml
<Default Extension="png" ContentType="image/png"/>
<Override PartName="/word/numbering.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.numbering+xml"/>
```

### 图片
**关键**：计算尺寸以防止页面溢出并保持宽高比。

```xml
<!-- Minimal required structure -->
<w:p>
  <w:r>
    <w:drawing>
      <wp:inline>
        <wp:extent cx="2743200" cy="1828800"/>
        <wp:docPr id="1" name="Picture 1"/>
        <a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">
          <a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">
            <pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">
              <pic:nvPicPr>
                <pic:cNvPr id="0" name="image1.png"/>
                <pic:cNvPicPr/>
              </pic:nvPicPr>
              <pic:blipFill>
                <a:blip r:embed="rId5"/>
                <!-- Add for stretch fill with aspect ratio preservation -->
                <a:stretch>
                  <a:fillRect/>
                </a:stretch>
              </pic:blipFill>
              <pic:spPr>
                <a:xfrm>
                  <a:ext cx="2743200" cy="1828800"/>
                </a:xfrm>
                <a:prstGeom prst="rect"/>
              </pic:spPr>
            </pic:pic>
          </a:graphicData>
        </a:graphic>
      </wp:inline>
    </w:drawing>
  </w:r>
</w:p>
```

### 链接（超链接）

**重要**：所有超链接（内部和外部）都需要在 styles.xml 中定义超链接样式。如果没有这种样式，链接将看起来像常规文本，而不是带有蓝色下划线的可点击链接。

**外部链接：**
```xml
<!-- In document.xml -->
<w:hyperlink r:id="rId5">
  <w:r>
    <w:rPr><w:rStyle w:val="Hyperlink"/></w:rPr>
    <w:t>Link Text</w:t>
  </w:r>
</w:hyperlink>

<!-- In word/_rels/document.xml.rels -->
<Relationship Id="rId5" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink" 
              Target="https://www.example.com/" TargetMode="External"/>
```

**内部链接：**

```xml
<!-- Link to bookmark -->
<w:hyperlink w:anchor="myBookmark">
  <w:r>
    <w:rPr><w:rStyle w:val="Hyperlink"/></w:rPr>
    <w:t>Link Text</w:t>
  </w:r>
</w:hyperlink>

<!-- Bookmark target -->
<w:bookmarkStart w:id="0" w:name="myBookmark"/>
<w:r><w:t>Target content</w:t></w:r>
<w:bookmarkEnd w:id="0"/>
```

**超链接样式（styles.xml 中必需）：**
```xml
<w:style w:type="character" w:styleId="Hyperlink">
  <w:name w:val="Hyperlink"/>
  <w:basedOn w:val="DefaultParagraphFont"/>
  <w:uiPriority w:val="99"/>
  <w:unhideWhenUsed/>
  <w:rPr>
    <w:color w:val="467886" w:themeColor="hyperlink"/>
    <w:u w:val="single"/>
  </w:rPr>
</w:style>
```

## 文档库 (Python)

使用 `scripts/document.py` 中的 Document 类来获取所有跟踪的更改和注释。它自动处理基础设施设置（people.xml、RSID、settings.xml、注释文件、关系、内容类型）。仅对库不支持的复杂场景使用直接 XML 操作。

**使用 Unicode 和实体：**
- **搜索**：实体表示法和 Unicode 字符均有效 - `contains="&#8220;Company"` 和 `contains="\u201cCompany"` 查找相同的文本
- **替换**：使用实体 (`&#8220;`) 或 Unicode (`\u201c`) - 两者都可以工作，并将根据文件的编码进行适当转换（ascii → 实体、utf-8 → Unicode）

### 初始化

**找到 docx 技能根**（包含 `scripts/` 和 `ooxml/` 的目录）：
```bash
# Search for document.py to locate the skill root
# Note: /mnt/skills is used here as an example; check your context for the actual location
find /mnt/skills -name "document.py" -path "*/docx/scripts/*" 2>/dev/null | head -1
# Example output: /mnt/skills/docx/scripts/document.py
# Skill root is: /mnt/skills/docx
```

**使用设置为 docx 技能根的 PYTHONPATH** 运行脚本：
```bash
PYTHONPATH=/mnt/skills/docx python your_script.py
```

**在您的脚本中**，从技能根导入：
```python
from scripts.document import Document, DocxXMLEditor

# Basic initialization (automatically creates temp copy and sets up infrastructure)
doc = Document('unpacked')

# Customize author and initials
doc = Document('unpacked', author="John Doe", initials="JD")

# Enable track revisions mode
doc = Document('unpacked', track_revisions=True)

# Specify custom RSID (auto-generated if not provided)
doc = Document('unpacked', rsid="07DC5ECB")
```

### 创建跟踪更改

**关键**：仅标记实际更改的文本。保留 `<w:del>`/`<w:ins>` 标记之外的所有未更改文本。标记未更改的文本会使编辑变得不专业并且更难以审核。

**属性处理**：Document 类自动将属性（w:id、w:date、w:rsidR、w:rsidDel、w16du:dateUtc、xml:space）注入到新元素中。保留原始文档中未更改的文本时，复制原始 `<w:r>` 元素及其现有属性以保持文档完整性。

**方法选择指南**：
- **添加您自己对常规文本的更改**：使用 `replace_node()` 和 `<w:del>`/`<w:ins>` 标签，或 `suggest_deletion()` 来删除整个 `<w:r>` 或 `<w:p>` 元素
- **部分修改其他作者跟踪的更改**：使用 `replace_node()` 将您的更改嵌套在他们的 `<w:ins>`/`<w:del>` 中
- **完全拒绝其他作者的插入**：在 `<w:ins>` 元素上使用 `revert_insertion()`（不是 `suggest_deletion()`）
- **完全拒绝其他作者的删除**：在 `<w:del>` 元素上使用 `revert_deletion()` 来使用跟踪更改恢复已删除的内容

```python
# Minimal edit - change one word: "The report is monthly" → "The report is quarterly"
# Original: <w:r w:rsidR="00AB12CD"><w:rPr><w:rFonts w:ascii="Calibri"/></w:rPr><w:t>The report is monthly</w:t></w:r>
node = doc["word/document.xml"].get_node(tag="w:r", contains="The report is monthly")
rpr = tags[0].toxml() if (tags := node.getElementsByTagName("w:rPr")) else ""
replacement = f'<w:r w:rsidR="00AB12CD">{rpr}<w:t>The report is </w:t></w:r><w:del><w:r>{rpr}<w:delText>monthly</w:delText></w:r></w:del><w:ins><w:r>{rpr}<w:t>quarterly</w:t></w:r></w:ins>'
doc["word/document.xml"].replace_node(node, replacement)

# Minimal edit - change number: "within 30 days" → "within 45 days"
# Original: <w:r w:rsidR="00XYZ789"><w:rPr><w:rFonts w:ascii="Calibri"/></w:rPr><w:t>within 30 days</w:t></w:r>
node = doc["word/document.xml"].get_node(tag="w:r", contains="within 30 days")
rpr = tags[0].toxml() if (tags := node.getElementsByTagName("w:rPr")) else ""
replacement = f'<w:r w:rsidR="00XYZ789">{rpr}<w:t>within </w:t></w:r><w:del><w:r>{rpr}<w:delText>30</w:delText></w:r></w:del><w:ins><w:r>{rpr}<w:t>45</w:t></w:r></w:ins><w:r w:rsidR="00XYZ789">{rpr}<w:t> days</w:t></w:r>'
doc["word/document.xml"].replace_node(node, replacement)

# Complete replacement - preserve formatting even when replacing all text
node = doc["word/document.xml"].get_node(tag="w:r", contains="apple")
rpr = tags[0].toxml() if (tags := node.getElementsByTagName("w:rPr")) else ""
replacement = f'<w:del><w:r>{rpr}<w:delText>apple</w:delText></w:r></w:del><w:ins><w:r>{rpr}<w:t>banana orange</w:t></w:r></w:ins>'
doc["word/document.xml"].replace_node(node, replacement)

# Insert new content (no attributes needed - auto-injected)
node = doc["word/document.xml"].get_node(tag="w:r", contains="existing text")
doc["word/document.xml"].insert_after(node, '<w:ins><w:r><w:t>new text</w:t></w:r></w:ins>')

# Partially delete another author's insertion
# Original: <w:ins w:author="Jane Smith" w:date="..."><w:r><w:t>quarterly financial report</w:t></w:r></w:ins>
# Goal: Delete only "financial" to make it "quarterly report"
node = doc["word/document.xml"].get_node(tag="w:ins", attrs={"w:id": "5"})
# IMPORTANT: Preserve w:author="Jane Smith" on the outer <w:ins> to maintain authorship
replacement = '''<w:ins w:author="Jane Smith" w:date="2025-01-15T10:00:00Z">
  <w:r><w:t>quarterly </w:t></w:r>
  <w:del><w:r><w:delText>financial </w:delText></w:r></w:del>
  <w:r><w:t>report</w:t></w:r>
</w:ins>'''
doc["word/document.xml"].replace_node(node, replacement)

# Change part of another author's insertion
# Original: <w:ins w:author="Jane Smith"><w:r><w:t>in silence, safe and sound</w:t></w:r></w:ins>
# Goal: Change "safe and sound" to "soft and unbound"
node = doc["word/document.xml"].get_node(tag="w:ins", attrs={"w:id": "8"})
replacement = f'''<w:ins w:author="Jane Smith" w:date="2025-01-15T10:00:00Z">
  <w:r><w:t>in silence, </w:t></w:r>
</w:ins>
<w:ins>
  <w:r><w:t>soft and unbound</w:t></w:r>
</w:ins>
<w:ins w:author="Jane Smith" w:date="2025-01-15T10:00:00Z">
  <w:del><w:r><w:delText>safe and sound</w:delText></w:r></w:del>
</w:ins>'''
doc["word/document.xml"].replace_node(node, replacement)

# Delete entire run (use only when deleting all content; use replace_node for partial deletions)
node = doc["word/document.xml"].get_node(tag="w:r", contains="text to delete")
doc["word/document.xml"].suggest_deletion(node)

# Delete entire paragraph (in-place, handles both regular and numbered list paragraphs)
para = doc["word/document.xml"].get_node(tag="w:p", contains="paragraph to delete")
doc["word/document.xml"].suggest_deletion(para)

# Add new numbered list item
target_para = doc["word/document.xml"].get_node(tag="w:p", contains="existing list item")
pPr = tags[0].toxml() if (tags := target_para.getElementsByTagName("w:pPr")) else ""
new_item = f'<w:p>{pPr}<w:r><w:t>New item</w:t></w:r></w:p>'
tracked_para = DocxXMLEditor.suggest_paragraph(new_item)
doc["word/document.xml"].insert_after(target_para, tracked_para)
# Optional: add spacing paragraph before content for better visual separation
# spacing = DocxXMLEditor.suggest_paragraph('<w:p><w:pPr><w:pStyle w:val="ListParagraph"/></w:pPr></w:p>')
# doc["word/document.xml"].insert_after(target_para, spacing + tracked_para)
```

### 添加评论

```python
# Add comment spanning two existing tracked changes
# Note: w:id is auto-generated. Only search by w:id if you know it from XML inspection
start_node = doc["word/document.xml"].get_node(tag="w:del", attrs={"w:id": "1"})
end_node = doc["word/document.xml"].get_node(tag="w:ins", attrs={"w:id": "2"})
doc.add_comment(start=start_node, end=end_node, text="Explanation of this change")

# Add comment on a paragraph
para = doc["word/document.xml"].get_node(tag="w:p", contains="paragraph text")
doc.add_comment(start=para, end=para, text="Comment on this paragraph")

# Add comment on newly created tracked change
# First create the tracked change
node = doc["word/document.xml"].get_node(tag="w:r", contains="old")
new_nodes = doc["word/document.xml"].replace_node(
    node,
    '<w:del><w:r><w:delText>old</w:delText></w:r></w:del><w:ins><w:r><w:t>new</w:t></w:r></w:ins>'
)
# Then add comment on the newly created elements
# new_nodes[0] is the <w:del>, new_nodes[1] is the <w:ins>
doc.add_comment(start=new_nodes[0], end=new_nodes[1], text="Changed old to new per requirements")

# Reply to existing comment
doc.reply_to_comment(parent_comment_id=0, text="I agree with this change")
```

### 拒绝跟踪更改

**重要**：使用 `revert_insertion()` 拒绝插入，使用 `revert_deletion()` 使用跟踪更改恢复删除。仅对常规未标记内容使用 `suggest_deletion()`。

```python
# Reject insertion (wraps it in deletion)
# Use this when another author inserted text that you want to delete
ins = doc["word/document.xml"].get_node(tag="w:ins", attrs={"w:id": "5"})
nodes = doc["word/document.xml"].revert_insertion(ins)  # Returns [ins]

# Reject deletion (creates insertion to restore deleted content)
# Use this when another author deleted text that you want to restore
del_elem = doc["word/document.xml"].get_node(tag="w:del", attrs={"w:id": "3"})
nodes = doc["word/document.xml"].revert_deletion(del_elem)  # Returns [del_elem, new_ins]

# Reject all insertions in a paragraph
para = doc["word/document.xml"].get_node(tag="w:p", contains="paragraph text")
nodes = doc["word/document.xml"].revert_insertion(para)  # Returns [para]

# Reject all deletions in a paragraph
para = doc["word/document.xml"].get_node(tag="w:p", contains="paragraph text")
nodes = doc["word/document.xml"].revert_deletion(para)  # Returns [para]
```

### 插入图像

**关键**：Document 类使用 `doc.unpacked_path` 处的临时副本。始终将图像复制到此临时目录，而不是原始解压文件夹。

```python
from PIL import Image
import shutil, os

# Initialize document first
doc = Document('unpacked')

# Copy image and calculate full-width dimensions with aspect ratio
media_dir = os.path.join(doc.unpacked_path, 'word/media')
os.makedirs(media_dir, exist_ok=True)
shutil.copy('image.png', os.path.join(media_dir, 'image1.png'))
img = Image.open(os.path.join(media_dir, 'image1.png'))
width_emus = int(6.5 * 914400)  # 6.5" usable width, 914400 EMUs/inch
height_emus = int(width_emus * img.size[1] / img.size[0])

# Add relationship and content type
rels_editor = doc['word/_rels/document.xml.rels']
next_rid = rels_editor.get_next_rid()
rels_editor.append_to(rels_editor.dom.documentElement,
    f'<Relationship Id="{next_rid}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/image1.png"/>')
doc['[Content_Types].xml'].append_to(doc['[Content_Types].xml'].dom.documentElement,
    '<Default Extension="png" ContentType="image/png"/>')

# Insert image
node = doc["word/document.xml"].get_node(tag="w:p", line_number=100)
doc["word/document.xml"].insert_after(node, f'''<w:p>
  <w:r>
    <w:drawing>
      <wp:inline distT="0" distB="0" distL="0" distR="0">
        <wp:extent cx="{width_emus}" cy="{height_emus}"/>
        <wp:docPr id="1" name="Picture 1"/>
        <a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">
          <a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">
            <pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">
              <pic:nvPicPr><pic:cNvPr id="1" name="image1.png"/><pic:cNvPicPr/></pic:nvPicPr>
              <pic:blipFill><a:blip r:embed="{next_rid}"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill>
              <pic:spPr><a:xfrm><a:ext cx="{width_emus}" cy="{height_emus}"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr>
            </pic:pic>
          </a:graphicData>
        </a:graphic>
      </wp:inline>
    </w:drawing>
  </w:r>
</w:p>''')
```

### 获取节点

```python
# By text content
node = doc["word/document.xml"].get_node(tag="w:p", contains="specific text")

# By line range
para = doc["word/document.xml"].get_node(tag="w:p", line_number=range(100, 150))

# By attributes
node = doc["word/document.xml"].get_node(tag="w:del", attrs={"w:id": "1"})

# By exact line number (must be line number where tag opens)
para = doc["word/document.xml"].get_node(tag="w:p", line_number=42)

# Combine filters
node = doc["word/document.xml"].get_node(tag="w:r", line_number=range(40, 60), contains="text")

# Disambiguate when text appears multiple times - add line_number range
node = doc["word/document.xml"].get_node(tag="w:r", contains="Section", line_number=range(2400, 2500))
```

### 保存

```python
# Save with automatic validation (copies back to original directory)
doc.save()  # Validates by default, raises error if validation fails

# Save to different location
doc.save('modified-unpacked')

# Skip validation (debugging only - needing this in production indicates XML issues)
doc.save(validate=False)
```

### 直接 DOM 操作

对于库未涵盖的复杂场景：

```python
# Access any XML file
editor = doc["word/document.xml"]
editor = doc["word/comments.xml"]

# Direct DOM access (defusedxml.minidom.Document)
node = doc["word/document.xml"].get_node(tag="w:p", line_number=5)
parent = node.parentNode
parent.removeChild(node)
parent.appendChild(node)  # Move to end

# General document manipulation (without tracked changes)
old_node = doc["word/document.xml"].get_node(tag="w:p", contains="original text")
doc["word/document.xml"].replace_node(old_node, "<w:p><w:r><w:t>replacement text</w:t></w:r></w:p>")

# Multiple insertions - use return value to maintain order
node = doc["word/document.xml"].get_node(tag="w:r", line_number=100)
nodes = doc["word/document.xml"].insert_after(node, "<w:r><w:t>A</w:t></w:r>")
nodes = doc["word/document.xml"].insert_after(nodes[-1], "<w:r><w:t>B</w:t></w:r>")
nodes = doc["word/document.xml"].insert_after(nodes[-1], "<w:r><w:t>C</w:t></w:r>")
# Results in: original_node, A, B, C
```

## 跟踪更改（红线）

**对所有跟踪的更改使用上面的 Document 类。** 下面的模式供构建替换 XML 字符串时参考。

### 验证规则
验证器在恢复 Claude 的更改后检查文档文本是否与原始文本匹配。这意味着：
- **切勿修改其他作者的 `<w:ins>` 或 `<w:del>` 标签内的文本**
- **始终使用嵌套删除**来删除其他作者的插入内容
- **必须使用 `<w:ins>` 或 `<w:del>` 标签正确跟踪每次编辑**

### 跟踪变更模式

**关键规则**：
1. 切勿修改其他作者跟踪更改中的内容。始终使用嵌套删除。
2. **XML 结构**：始终将 `<w:del>` 和 `<w:ins>` 放置在包含完整 `<w:r>` 元素的段落级别。切勿嵌套在 `<w:r>` 元素内 - 这会创建无效的 XML，从而中断文档处理。

**文本插入：**
```xml
<w:ins w:id="1" w:author="Claude" w:date="2025-07-30T23:05:00Z" w16du:dateUtc="2025-07-31T06:05:00Z">
  <w:r w:rsidR="00792858">
    <w:t>inserted text</w:t>
  </w:r>
</w:ins>
```

**文本删除：**
```xml
<w:del w:id="2" w:author="Claude" w:date="2025-07-30T23:05:00Z" w16du:dateUtc="2025-07-31T06:05:00Z">
  <w:r w:rsidDel="00792858">
    <w:delText>deleted text</w:delText>
  </w:r>
</w:del>
```

**删除其他作者的插入（必须使用嵌套结构）：**
```xml
<!-- Nest deletion inside the original insertion -->
<w:ins w:author="Jane Smith" w:id="16">
  <w:del w:author="Claude" w:id="40">
    <w:r><w:delText>monthly</w:delText></w:r>
  </w:del>
</w:ins>
<w:ins w:author="Claude" w:id="41">
  <w:r><w:t>weekly</w:t></w:r>
</w:ins>
```

**恢复其他作者的删除：**
```xml
<!-- Leave their deletion unchanged, add new insertion after it -->
<w:del w:author="Jane Smith" w:id="50">
  <w:r><w:delText>within 30 days</w:delText></w:r>
</w:del>
<w:ins w:author="Claude" w:id="51">
  <w:r><w:t>within 30 days</w:t></w:r>
</w:ins>
```