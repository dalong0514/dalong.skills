<!-- 此文件由机器翻译自 pdf_features.md -->

# PDF 功能参考

高级 PDF 功能：链接、书签、表单、加密和元数据。

## 文档元数据

设置可在 PDF 阅读器中查看的 PDF 文档属性。

```python
from reportlab.pdfgen import canvas

c = canvas.Canvas("output.pdf")

# Set metadata
c.setAuthor("John Doe")
c.setTitle("Annual Report 2024")
c.setSubject("Financial Analysis")
c.setKeywords("finance, annual, report, 2024")
c.setCreator("MyApp v1.0")

# ... draw content ...

c.save()
```

与鸭嘴兽：

<<<代码块_1>>>

## 书签和目的地

创建内部导航结构。

### 简单书签

<<<代码块_2>>>

### 书签级别

<<<代码块_3>>>

### 目的地适配模式

控制导航时页面的显示方式：

<<<代码块_4>>>

## 超链接

### 外部链接

<<<代码块_5>>>

### 内部链接

链接到文档中添加书签的位置：

<<<代码块_6>>>

### 段落中的链接

对于 Platypus 文档：

```python
from reportlab.platypus import Paragraph

# External link
text = '<link href="https://example.com" color="blue">Visit our website</link>'
para = Paragraph(text, style)

# Internal link (to anchor)
text = '<link href="#section1" color="blue">Go to Section 1</link>'
para1 = Paragraph(text, style)

# Create anchor
text = '<a name="section1"/>Section 1 Heading'
para2 = Paragraph(text, heading_style)

story.append(para1)
story.append(para2)
```

## 互动表格

创建可填写的 PDF 表单。

### 文本字段

```python
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfform
from reportlab.lib.colors import black, white

c = canvas.Canvas("form.pdf")

# Create text field
c.acroForm.textfield(
    name="name",
    tooltip="Enter your name",
    x=100,
    y=700,
    width=200,
    height=20,
    borderColor=black,
    fillColor=white,
    textColor=black,
    forceBorder=True,
    fontSize=12,
    maxlen=100,  # Maximum character length
)

# Label
c.drawString(100, 725, "Name:")

c.save()
```

### 复选框

```python
# Create checkbox
c.acroForm.checkbox(
    name="agree",
    tooltip="I agree to terms",
    x=100,
    y=650,
    size=20,
    buttonStyle='check',  # 'check', 'circle', 'cross', 'diamond', 'square', 'star'
    borderColor=black,
    fillColor=white,
    textColor=black,
    forceBorder=True,
    checked=False,  # Initial state
)

c.drawString(130, 655, "I agree to the terms and conditions")
```

### 单选按钮

```python
# Radio button group - only one can be selected
c.acroForm.radio(
    name="payment",  # Same name for group
    tooltip="Credit Card",
    value="credit",  # Value when selected
    x=100,
    y=600,
    size=15,
    selected=False,
)
c.drawString(125, 603, "Credit Card")

c.acroForm.radio(
    name="payment",  # Same name
    tooltip="PayPal",
    value="paypal",
    x=100,
    y=580,
    size=15,
    selected=False,
)
c.drawString(125, 583, "PayPal")
```

### 列表框

```python
# Listbox with multiple options
c.acroForm.listbox(
    name="country",
    tooltip="Select your country",
    value="US",  # Default selected
    x=100,
    y=500,
    width=150,
    height=80,
    borderColor=black,
    fillColor=white,
    textColor=black,
    forceBorder=True,
    options=[
        ("United States", "US"),
        ("Canada", "CA"),
        ("Mexico", "MX"),
        ("Other", "OTHER"),
    ],  # List of (label, value) tuples
    multiple=False,  # Allow multiple selections
)
```

### 选择（下拉菜单）

```python
# Dropdown menu
c.acroForm.choice(
    name="state",
    tooltip="Select state",
    value="CA",
    x=100,
    y=450,
    width=150,
    height=20,
    borderColor=black,
    fillColor=white,
    textColor=black,
    forceBorder=True,
    options=[
        ("California", "CA"),
        ("New York", "NY"),
        ("Texas", "TX"),
    ],
)
```

### 完整表格示例

```python
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from reportlab.lib.colors import black, white, lightgrey
from reportlab.lib.units import inch

def create_registration_form(filename):
    c = canvas.Canvas(filename, pagesize=letter)
    c.setFont("Helvetica-Bold", 16)
    c.drawString(inch, 10*inch, "Registration Form")

    y = 9*inch
    c.setFont("Helvetica", 12)

    # Name field
    c.drawString(inch, y, "Full Name:")
    c.acroForm.textfield(
        name="fullname",
        x=2*inch, y=y-5, width=4*inch, height=20,
        borderColor=black, fillColor=lightgrey, forceBorder=True
    )

    # Email field
    y -= 0.5*inch
    c.drawString(inch, y, "Email:")
    c.acroForm.textfield(
        name="email",
        x=2*inch, y=y-5, width=4*inch, height=20,
        borderColor=black, fillColor=lightgrey, forceBorder=True
    )

    # Age dropdown
    y -= 0.5*inch
    c.drawString(inch, y, "Age Group:")
    c.acroForm.choice(
        name="age_group",
        x=2*inch, y=y-5, width=2*inch, height=20,
        borderColor=black, fillColor=lightgrey, forceBorder=True,
        options=[("18-25", "18-25"), ("26-35", "26-35"),
                ("36-50", "36-50"), ("51+", "51+")]
    )

    # Newsletter checkbox
    y -= 0.5*inch
    c.acroForm.checkbox(
        name="newsletter",
        x=inch, y=y-5, size=15,
        buttonStyle='check', borderColor=black, forceBorder=True
    )
    c.drawString(inch + 25, y, "Subscribe to newsletter")

    c.save()

create_registration_form("registration.pdf")
```

## 加密和安全

使用密码和权限保护 PDF。

### 基本加密

```python
from reportlab.pdfgen import canvas

c = canvas.Canvas("secure.pdf")

# Encrypt with user password
c.encrypt(
    userPassword="user123",    # Password to open
    ownerPassword="owner456",  # Password to change permissions
    canPrint=1,                # Allow printing
    canModify=0,               # Disallow modifications
    canCopy=1,                 # Allow text copying
    canAnnotate=0,             # Disallow annotations
    strength=128,              # 40 or 128 bit encryption
)

# ... draw content ...

c.save()
```

### 权限设置

```python
c.encrypt(
    userPassword="user123",
    ownerPassword="owner456",
    canPrint=1,        # 1 = allow, 0 = deny
    canModify=0,       # Prevent content modification
    canCopy=1,         # Allow text/graphics copying
    canAnnotate=0,     # Prevent comments/annotations
    strength=128,      # Use 128-bit encryption
)
```

### 高级加密

```python
from reportlab.lib.pdfencrypt import StandardEncryption

# Create encryption object
encrypt = StandardEncryption(
    userPassword="user123",
    ownerPassword="owner456",
    canPrint=1,
    canModify=0,
    canCopy=1,
    canAnnotate=1,
    strength=128,
)

# Use with canvas
c = canvas.Canvas("secure.pdf")
c._doc.encrypt = encrypt

# ... draw content ...

c.save()
```

### 带加密的鸭嘴兽

```python
from reportlab.platypus import SimpleDocTemplate

doc = SimpleDocTemplate("secure.pdf")

# Set encryption
doc.encrypt = True
doc.canPrint = 1
doc.canModify = 0

# Or use encrypt() method
doc.encrypt = encrypt_object

doc.build(story)
```

## 页面转换

为演示文稿添加视觉效果。

```python
from reportlab.pdfgen import canvas

c = canvas.Canvas("presentation.pdf")

# Set transition for current page
c.setPageTransition(
    effectname="Wipe",  # Transition effect
    duration=1,         # Duration in seconds
    direction=0         # Direction (effect-specific)
)

# Available effects:
# "Split", "Blinds", "Box", "Wipe", "Dissolve",
# "Glitter", "R" (Replace), "Fly", "Push", "Cover",
# "Uncover", "Fade"

# Direction values (effect-dependent):
# 0, 90, 180, 270 for most directional effects

# Example: Slide with fade transition
c.setFont("Helvetica-Bold", 24)
c.drawString(100, 400, "Slide 1")
c.setPageTransition("Fade", 0.5)
c.showPage()

c.drawString(100, 400, "Slide 2")
c.setPageTransition("Wipe", 1, 90)
c.showPage()

c.save()
```

## PDF/A 合规性

创建档案质量的 PDF。

```python
from reportlab.pdfgen import canvas

c = canvas.Canvas("pdfa.pdf")

# Enable PDF/A-1b compliance
c.setPageCompression(0)  # PDF/A requires uncompressed
# Note: Full PDF/A requires additional XMP metadata
# This is simplified - full compliance needs more setup

# ... draw content ...

c.save()
```

## 压缩

控制文件大小与生成速度。

```python
# Enable page compression
c = canvas.Canvas("output.pdf", pageCompression=1)

# Compression reduces file size but slows generation
# 0 = no compression (faster, larger files)
# 1 = compression (slower, smaller files)
```

## 表单和 XObject

可重复使用的图形元素。

```python
from reportlab.pdfgen import canvas

c = canvas.Canvas("output.pdf")

# Begin form (reusable object)
c.beginForm("logo")
c.setFillColorRGB(0, 0, 1)
c.rect(0, 0, 100, 50, fill=1)
c.setFillColorRGB(1, 1, 1)
c.drawString(10, 20, "LOGO")
c.endForm()

# Use form multiple times
c.doForm("logo")  # At current position
c.translate(200, 0)
c.doForm("logo")  # At translated position
c.translate(200, 0)
c.doForm("logo")

c.save()

# Benefits: Smaller file size, faster rendering
```

## 最佳实践

1. **始终为专业文档设置元数据**
2. **对超过 10 页的文档使用书签**
3. **使链接在视觉上清晰**（蓝色，下划线）
4. **在多个 PDF 阅读器中测试表单**（行为各不相同）
5. **对敏感数据使用强加密（128 位）**
6. **设置用户和所有者密码**以确保完全安全
7. **启用打印**，除非有特别限制
8. **测试页面转换** - 一些读者不支持所有效果
9. **使用有意义的书签标题**进行导航
10. **考虑使用 PDF/A** 来满足长期存档需求
11. **验证表单字段名称** - 必须是唯一且有效的标识符
12. **添加工具提示**以形成字段以获得更好的用户体验