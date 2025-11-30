<!-- 此文件由机器翻译自 barcodes_reference.md -->

# 条形码参考

在 ReportLab 中创建条形码和 QR 码的综合指南。

## 可用的条形码类型

ReportLab 支持多种一维和二维条形码格式。

### 一维条形码（线性）

- **Code128** - 紧凑，编码完整 ASCII
- **Code39** (Standard39) - 字母数字，广泛支持
- **Code93** (Standard93) - 压缩 Code39
- **EAN-13** - 欧洲商品编号（零售）
- **EAN-8** - EAN 的缩写形式
- **EAN-5** - 5 位数附加（定价）
- **UPC-A** - 通用产品代码（北美）
- **ISBN** - 国际标准书号
- **代码11** - 电信
- **Codabar** - 血库、联邦快递、图书馆
- **I2of5** (Interleaved 2 of 5) - 仓库/配送
- **MSI** - 库存控制
- **POSTNET** - 美国邮政服务
- **USPS_4State** - 美国邮政服务
- **FIM** (A、B、C、D) - 正面识别标记（邮件分拣）

### 二维条形码

- **QR** - QR 码（广泛用于 URL、联系信息）
- **ECC200DataMatrix** - 数据矩阵格式

## 将条形码与画布结合使用

### Code128（推荐一般用途）

Code128 用途广泛且结构紧凑 - 使用强制校验和对完整 ASCII 字符集进行编码。

```python
from reportlab.pdfgen import canvas
from reportlab.graphics.barcode import code128
from reportlab.lib.units import inch

c = canvas.Canvas("barcode.pdf")

# Create barcode
barcode = code128.Code128("HELLO123")

# Draw on canvas
barcode.drawOn(c, 1*inch, 5*inch)

c.save()
```

### Code128 选项

<<<代码块_1>>>

### Code39（标准39）

支持：0-9、A-Z（大写）、空格和特殊字符 (-.$/+%*)。

<<<代码块_2>>>

### 扩展代码39

对完整 ASCII（Code39 字符对）进行编码。

<<<代码块_3>>>

### 代码93

<<<代码块_4>>>

### EAN-13（欧洲商品编号）

零售产品的 13 位条形码。

<<<代码块_5>>>

### EAN-8

缩写形式，8 位数字。

<<<代码块_6>>>

### UPC-A

北美使用的 12 位条形码。

```python
from reportlab.graphics.barcode import usps

# 11 digits (12th is checksum)
barcode = usps.UPCA(
    value="01234567890"
)

barcode.drawOn(canvas, x, y)
```

### ISBN（书籍）

```python
from reportlab.graphics.barcode.widgets import ISBNBarcodeWidget

# 10 or 13 digit ISBN
barcode = ISBNBarcodeWidget(
    value="978-0-123456-78-9"
)

# With pricing (EAN-5 add-on)
barcode = ISBNBarcodeWidget(
    value="978-0-123456-78-9",
    price=True,
)
```

### 二维码

最通用的二维条形码 - 可以编码 URL、文本、联系信息等。

```python
from reportlab.graphics.barcode.qr import QrCodeWidget
from reportlab.graphics.shapes import Drawing
from reportlab.graphics import renderPDF

# Create QR code
qr = QrCodeWidget("https://example.com")

# Size in pixels (QR codes are square)
qr.barWidth = 100  # Width in points
qr.barHeight = 100  # Height in points

# Error correction level
# L = 7% recovery, M = 15%, Q = 25%, H = 30%
qr.qrVersion = 1  # Auto-size (1-40, or None for auto)
qr.errorLevel = 'M'  # L, M, Q, H

# Draw
d = Drawing()
d.add(qr)
renderPDF.draw(d, canvas, x, y)
```

### 二维码 - 更多选项

```python
# URL QR Code
qr = QrCodeWidget("https://example.com")

# Contact information (vCard)
vcard_data = """BEGIN:VCARD
VERSION:3.0
FN:John Doe
TEL:+1-555-1234
EMAIL:john@example.com
END:VCARD"""
qr = QrCodeWidget(vcard_data)

# WiFi credentials
wifi_data = "WIFI:T:WPA;S:NetworkName;P:Password;;"
qr = QrCodeWidget(wifi_data)

# Plain text
qr = QrCodeWidget("Any text here")
```

### 数据矩阵 (ECC200)

适用于小物品的紧凑型二维条码。

```python
from reportlab.graphics.barcode.datamatrix import DataMatrixWidget

barcode = DataMatrixWidget(
    value="DATA123"
)

d = Drawing()
d.add(barcode)
renderPDF.draw(d, canvas, x, y)
```

### 邮政条形码

```python
from reportlab.graphics.barcode import usps

# POSTNET (older format)
barcode = usps.POSTNET(
    value="55555-1234",  # ZIP or ZIP+4
)

# USPS 4-State (newer)
barcode = usps.USPS_4State(
    value="12345678901234567890",  # 20-digit routing code
    routing="12345678901"
)

barcode.drawOn(canvas, x, y)
```

### FIM（正面识别标记）

用于邮件分拣。

```python
from reportlab.graphics.barcode import usps

# FIM-A, FIM-B, FIM-C, or FIM-D
barcode = usps.FIM(
    value="A"  # A, B, C, or D
)

barcode.drawOn(canvas, x, y)
```

## 在 Platypus 中使用条形码

对于流动文档，请将条形码包装在 Flowables 中。

### 简单方法 - 绘制 Flowable

```python
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.barcode.qr import QrCodeWidget
from reportlab.lib.units import inch

# Create drawing
d = Drawing(2*inch, 2*inch)

# Create barcode
qr = QrCodeWidget("https://example.com")
qr.barWidth = 2*inch
qr.barHeight = 2*inch
qr.x = 0
qr.y = 0

d.add(qr)

# Add to story
story.append(d)
```

### 定制流动包装

```python
from reportlab.platypus import Flowable
from reportlab.graphics.barcode import code128
from reportlab.lib.units import inch

class BarcodeFlowable(Flowable):
    def __init__(self, code, barcode_type='code128', width=2*inch, height=0.5*inch):
        Flowable.__init__(self)
        self.code = code
        self.barcode_type = barcode_type
        self.width_val = width
        self.height_val = height

        # Create barcode
        if barcode_type == 'code128':
            self.barcode = code128.Code128(code, barWidth=width/100, barHeight=height)
        # Add other types as needed

    def draw(self):
        self.barcode.drawOn(self.canv, 0, 0)

    def wrap(self, availWidth, availHeight):
        return (self.barcode.width, self.barcode.height)

# Use in story
story.append(BarcodeFlowable("PRODUCT123"))
```

## 完整示例

### 带条形码的产品标签

```python
from reportlab.pdfgen import canvas
from reportlab.graphics.barcode import code128
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch

def create_product_label(filename, product_code, product_name):
    c = canvas.Canvas(filename, pagesize=(4*inch, 2*inch))

    # Product name
    c.setFont("Helvetica-Bold", 14)
    c.drawCentredString(2*inch, 1.5*inch, product_name)

    # Barcode
    barcode = code128.Code128(product_code)
    barcode_width = barcode.width
    barcode_height = barcode.height

    # Center barcode
    x = (4*inch - barcode_width) / 2
    y = 0.5*inch

    barcode.drawOn(c, x, y)

    # Code text
    c.setFont("Courier", 10)
    c.drawCentredString(2*inch, 0.3*inch, product_code)

    c.save()

create_product_label("label.pdf", "ABC123456789", "Premium Widget")
```

### 二维码联系卡

```python
from reportlab.pdfgen import canvas
from reportlab.graphics.barcode.qr import QrCodeWidget
from reportlab.graphics.shapes import Drawing
from reportlab.graphics import renderPDF
from reportlab.lib.units import inch

def create_contact_card(filename, name, phone, email):
    c = canvas.Canvas(filename, pagesize=(3.5*inch, 2*inch))

    # Contact info
    c.setFont("Helvetica-Bold", 12)
    c.drawString(0.5*inch, 1.5*inch, name)
    c.setFont("Helvetica", 10)
    c.drawString(0.5*inch, 1.3*inch, phone)
    c.drawString(0.5*inch, 1.1*inch, email)

    # Create vCard data
    vcard = f"""BEGIN:VCARD
VERSION:3.0
FN:{name}
TEL:{phone}
EMAIL:{email}
END:VCARD"""

    # QR code
    qr = QrCodeWidget(vcard)
    qr.barWidth = 1.5*inch
    qr.barHeight = 1.5*inch

    d = Drawing()
    d.add(qr)

    renderPDF.draw(d, c, 1.8*inch, 0.2*inch)

    c.save()

create_contact_card("contact.pdf", "John Doe", "+1-555-1234", "john@example.com")
```

### 带有多个条形码的运输标签

```python
from reportlab.pdfgen import canvas
from reportlab.graphics.barcode import code128
from reportlab.lib.units import inch

def create_shipping_label(filename, tracking_code, zip_code):
    c = canvas.Canvas(filename, pagesize=(6*inch, 4*inch))

    # Title
    c.setFont("Helvetica-Bold", 16)
    c.drawString(0.5*inch, 3.5*inch, "SHIPPING LABEL")

    # Tracking barcode
    c.setFont("Helvetica", 10)
    c.drawString(0.5*inch, 2.8*inch, "Tracking Number:")

    tracking_barcode = code128.Code128(tracking_code, barHeight=0.5*inch)
    tracking_barcode.drawOn(c, 0.5*inch, 2*inch)

    c.setFont("Courier", 9)
    c.drawString(0.5*inch, 1.8*inch, tracking_code)

    # Additional info can be added

    c.save()

create_shipping_label("shipping.pdf", "1Z999AA10123456784", "12345")
```

## 条码选择指南

**在以下情况下选择 Code128：**
- 通用编码
- 需要对数字和字母进行编码
- 想要紧凑的尺寸
- 广泛支持

**在以下情况下选择 Code39：**
- 较旧的系统需要它
- 不需要小写字母
- 想要最大的兼容性

**在以下情况下选择二维码：**
- 需要对 URL 进行编码
- 想要移动设备扫描
- 需要高数据容量
- 想要纠错

**在以下情况下选择 EAN/UPC：**
- 零售产品标识
- 需要行业标准格式
- 全球分布

**在以下情况下选择 Data Matrix：**
- 空间非常有限
- 小件物品（PCB、电子产品）
- 需要 2D 紧凑格式

## 最佳实践

1. **尽早使用实际的条形码扫描仪/读取器测试扫描**
2. **在条形码周围添加安静区域**（空白区域） - 设置 `quiet=1`
3. **选择合适的高度** - 较高的条码更容易扫描
4. **在条形码下方包含人类可读的文本**以供手动输入
5. **使用 Code128** 作为通用用途的默认值 - 它紧凑且多功能
6. **对于 URL，使用 QR 码** - 对于移动用户来说更容易
7. **检查您所在行业的条形码标准**（零售业使用 EAN/UPC）
8. **测试打印质量** - 低 DPI 会使条形码无法扫描
9. **编码前验证数据** - 错误的校验位会导致问题
10. **考虑对 QR 码进行纠错** - 对重要数据使用“M”或“H”