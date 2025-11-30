<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：matplotlib
描述：“基础绘图库。创建线图、散点图、条形图、直方图、热图、3D、子图，导出 PNG/PDF/SVG，用于科学可视化和出版数据。”
---

# Matplotlib

## 概述

Matplotlib 是 Python 的基础可视化库，用于创建静态、动画和交互式绘图。该技能提供了有效使用 matplotlib 的指导，涵盖 pyplot 界面（MATLAB 风格）和面向对象的 API（图/轴），以及创建出版质量可视化的最佳实践。

## 何时使用此技能

该技能应该在以下情况下使用：
- 创建任何类型的绘图或图表（折线图、散点图、条形图、直方图、热图、轮廓图等）
- 生成科学或统计可视化
- 自定义绘图外观（颜色、样式、标签、图例）
- 创建带有子图的多面板图形
- 将可视化导出为各种格式（PNG、PDF、SVG 等）
- 构建交互式情节或动画
- 使用 3D 可视化
- 将绘图集成到 Jupyter 笔记本或 GUI 应用程序中

## 核心概念

### Matplotlib 层次结构

Matplotlib 使用对象的层次结构：

1. **图** - 所有绘图元素的顶级容器
2. **Axes** - 显示数据的实际绘图区域（一个图可以包含多个Axes）
3. **艺术家** - 图形上可见的所有内容（线条、文本、刻度线等）
4. **Axis** - 处理刻度和标签的数轴对象（x轴、y轴）

### 两个接口

**1. pyplot 接口（隐式，MATLAB 风格）**
```python
import matplotlib.pyplot as plt

plt.plot([1, 2, 3, 4])
plt.ylabel('some numbers')
plt.show()
```
- 方便快速、简单地绘图
- 自动维护状态
- 适合交互式工作和简单脚本

**2.面向对象的接口（显式）**
<<<代码块_1>>>
- **推荐用于大多数用例**
- 对图形和轴进行更明确的控制
- 更适合具有多个子图的复杂图形
- 更容易维护和调试

## 常见工作流程

### 1. 基本情节创建

**单图工作流程：**
<<<代码块_2>>>

### 2. 多个子图

**创建子图布局：**
<<<代码块_3>>>

### 3. 绘图类型和用例

**线图** - 时间序列、连续数据、趋势
<<<代码块_4>>>

**散点图** - 变量之间的关系、相关性
<<<代码块_5>>>

**条形图** - 分类比较
<<<代码块_6>>>

**直方图** - 分布
```python
ax.hist(data, bins=30, edgecolor='black', alpha=0.7)
```

**热图** - 矩阵数据、相关性
```python
im = ax.imshow(matrix, cmap='coolwarm', aspect='auto')
plt.colorbar(im, ax=ax)
```

**等高线图** - 2D 平面上的 3D 数据
```python
contour = ax.contour(X, Y, Z, levels=10)
ax.clabel(contour, inline=True, fontsize=8)
```

**箱线图** - 统计分布
```python
ax.boxplot([data1, data2, data3], labels=['A', 'B', 'C'])
```

**小提琴图** - 分布密度
```python
ax.violinplot([data1, data2, data3], positions=[1, 2, 3])
```

有关全面的绘图类型示例和变体，请参阅`references/plot_types.md`。

### 4. 样式和定制

**颜色指定方法：**
- 命名颜色：`'red'`、`'blue'`、`'steelblue'`
- 十六进制代码：`'#FF5733'`
- RGB 元组：`(0.1, 0.2, 0.3)`
- 颜色图：`cmap='viridis'`、`cmap='plasma'`、`cmap='coolwarm'`

**使用样式表：**
```python
plt.style.use('seaborn-v0_8-darkgrid')  # Apply predefined style
# Available styles: 'ggplot', 'bmh', 'fivethirtyeight', etc.
print(plt.style.available)  # List all available styles
```

**使用 rcParams 进行定制：**
```python
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['figure.titlesize'] = 18
```

**文字和注释：**
```python
ax.text(x, y, 'annotation', fontsize=12, ha='center')
ax.annotate('important point', xy=(x, y), xytext=(x+1, y+1),
            arrowprops=dict(arrowstyle='->', color='red'))
```

有关详细的样式选项和颜色图指南，请参阅`references/styling_guide.md`。

### 5. 保存数字

**导出为各种格式：**
```python
# High-resolution PNG for presentations/papers
plt.savefig('figure.png', dpi=300, bbox_inches='tight', facecolor='white')

# Vector format for publications (scalable)
plt.savefig('figure.pdf', bbox_inches='tight')
plt.savefig('figure.svg', bbox_inches='tight')

# Transparent background
plt.savefig('figure.png', dpi=300, bbox_inches='tight', transparent=True)
```

**重要参数：**
- `dpi`：分辨率（出版物为 300，网络为 150，屏幕为 72）
- `bbox_inches='tight'`：删除多余的空格
- `facecolor='white'`：确保白色背景（对于透明主题有用）
- `transparent=True`：透明背景

### 6. 使用 3D 绘图

```python
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Surface plot
ax.plot_surface(X, Y, Z, cmap='viridis')

# 3D scatter
ax.scatter(x, y, z, c=colors, marker='o')

# 3D line plot
ax.plot(x, y, z, linewidth=2)

# Labels
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
```

## 最佳实践

### 1. 接口选择
- **使用面向对象的接口**（fig, ax = plt.subplots()）用于生产代码
- 保留 pyplot 接口仅用于快速交互式探索
- 始终显式创建图形而不是依赖隐式状态

### 2. 图形尺寸和 DPI
- 在创建时设置Figsize：`fig, ax = plt.subplots(figsize=(10, 6))`
- 对输出介质使用适当的 DPI：
  - 屏幕/笔记本电脑：72-100 dpi
  - 网页：150 dpi
  - 印刷/出版物：300 dpi

### 3.布局管理
- 使用`constrained_layout=True`或`tight_layout()`来防止元素重叠
- 建议使用 `fig, ax = plt.subplots(constrained_layout=True)` 自动间距

### 4. 颜色图选择
- **顺序**（viridis、plasma、inferno）：具有一致进展的有序数据
- **发散**（coolwarm、RdBu）：具有有意义的中心点（例如零）的数据
- **定性**（tab10，Set3）：分类/名义数据
- 避免彩虹色图（jet）——它们在感知上并不统一

### 5. 辅助功能
- 使用色盲友好的色彩图（viridis、cividis）
- 除了颜色之外还为条形图添加图案/阴影线
- 确保元素之间有足够的对比度
- 包括描述性标签和图例

### 6. 性能
- 对于大型数据集，在绘图调用中使用 `rasterized=True` 来减小文件大小
- 在绘图之前使用适当的数据缩减（例如，下采样密集时间序列）
- 对于动画，使用位图传输以获得更好的性能

### 7. 代码组织
```python
# Good practice: Clear structure
def create_analysis_plot(data, title):
    """Create standardized analysis plot."""
    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)

    # Plot data
    ax.plot(data['x'], data['y'], linewidth=2)

    # Customize
    ax.set_xlabel('X Axis Label', fontsize=12)
    ax.set_ylabel('Y Axis Label', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)

    return fig, ax

# Use the function
fig, ax = create_analysis_plot(my_data, 'My Analysis')
plt.savefig('analysis.png', dpi=300, bbox_inches='tight')
```

## 快速参考脚本

此技能包括 `scripts/` 目录中的帮助程序脚本：

### `plot_template.py`
模板脚本展示了各种绘图类型和最佳实践。以此作为创建新可视化的起点。

**用途：**
```bash
python scripts/plot_template.py
```

### `style_configurator.py`
用于配置 matplotlib 样式首选项并生成自定义样式表的交互式实用程序。

**用途：**
```bash
python scripts/style_configurator.py
```

## 详细参考资料

如需全面信息，请参阅参考文档：

- **`references/plot_types.md`** - 包含代码示例和用例的完整绘图类型目录
- **`references/styling_guide.md`** - 详细的样式选项、颜色图和自定义
- **`references/api_reference.md`** - 核心类和方法参考
- **`references/common_issues.md`** - 常见问题的故障排除指南

## 与其他工具集成

Matplotlib 与以下各项集成良好：
- **NumPy/Pandas** - 从数组和 DataFrame 直接绘图
- **Seaborn** - 基于 matplotlib 构建的高级统计可视化
- **Jupyter** - 使用 `%matplotlib inline` 或 `%matplotlib widget` 进行交互式绘图
- **GUI 框架** - 嵌入 Tkinter、Qt、wxPython 应用程序

## 常见问题

1. **重叠元素**：使用 `constrained_layout=True` 或 `tight_layout()`
2. **状态混乱**：使用OO接口避免pyplot状态机问题
3. **许多图形的内存问题**：使用 `plt.close(fig)` 显式关闭图形
4. **字体警告**：安装字体或使用 `plt.rcParams['font.sans-serif']` 抑制警告
5. **DPI 混淆**：请记住，figsize 的单位是英寸，而不是像素：`pixels = dpi * inches`

## 其他资源

- 官方文档：https://matplotlib.org/
- 图库：https://matplotlib.org/stable/gallery/index.html
- 备忘单：https://matplotlib.org/cheatsheets/
- 教程：https://matplotlib.org/stable/tutorials/index.html