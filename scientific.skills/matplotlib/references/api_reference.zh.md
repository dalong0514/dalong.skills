<!-- 此文件由机器翻译自 api_reference.md -->

# Matplotlib API 参考

本文档提供了最常用的 matplotlib 类和方法的快速参考。

## 核心课程

### 图

所有绘图元素的顶级容器。

**创作：**
```python
fig = plt.figure(figsize=(10, 6), dpi=100, facecolor='white')
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
```

**关键方法：**
- `fig.add_subplot(nrows, ncols, index)` - 添加子图
- `fig.add_axes([left, bottom, width, height])` - 在特定位置添加轴
- `fig.savefig(filename, dpi=300, bbox_inches='tight')` - 保存图形
- `fig.tight_layout()` - 调整间距以防止重叠
- `fig.suptitle(title)` - 设置图形标题
- `fig.legend()` - 创建图形级图例
- `fig.colorbar(mappable)` - 将颜色条添加到图中
- `plt.close(fig)` - 关闭图形以释放内存

**关键属性：**
- `fig.axes` - 图中所有轴的列表
- `fig.dpi` - 分辨率（以每英寸点数为单位）
- `fig.figsize` - 图形尺寸以英寸为单位（宽度、高度）

### 轴

数据可视化的实际绘图区域。

**创作：**
<<<代码块_1>>>

**绘图方法：**

**线图：**
- `ax.plot(x, y, **kwargs)` - 线图
- `ax.step(x, y, where='pre'/'mid'/'post')` - 步进图
- `ax.errorbar(x, y, yerr, xerr)` - 误差线

**散点图：**
- `ax.scatter(x, y, s=size, c=color, marker='o', alpha=0.5)` - 散点图

**条形图：**
- `ax.bar(x, height, width=0.8, align='center')` - 垂直条形图
- `ax.barh(y, width)` - 水平条形图

**统计图：**
- `ax.hist(data, bins=10, density=False)` - 直方图
- `ax.boxplot(data, labels=None)` - 箱线图
- `ax.violinplot(data)` - 小提琴情节

**二维图：**
- `ax.imshow(array, cmap='viridis', aspect='auto')` - 显示图像/矩阵
- `ax.contour(X, Y, Z, levels=10)` - 轮廓线
- `ax.contourf(X, Y, Z, levels=10)` - 填充轮廓
- `ax.pcolormesh(X, Y, Z)` - 伪彩色图

**填充物：**
- `ax.fill_between(x, y1, y2, alpha=0.3)` - 在曲线之间填充
- `ax.fill_betweenx(y, x1, x2)` - 在垂直曲线之间填充

**文字和注释：**
- `ax.text(x, y, text, fontsize=12)` - 添加文本
- `ax.annotate(text, xy=(x, y), xytext=(x2, y2), arrowprops={})` - 用箭头注释

**定制方式：**

**标签和标题：**
- `ax.set_xlabel(label, fontsize=12)` - 设置 x 轴标签
- `ax.set_ylabel(label, fontsize=12)` - 设置 y 轴标签
- `ax.set_title(title, fontsize=14)` - 设置轴标题

**限制和范围：**
- `ax.set_xlim(left, right)` - 设置 x 轴限制
- `ax.set_ylim(bottom, top)` - 设置 y 轴限制
- `ax.set_xscale('linear'/'log'/'symlog')` - 设置 x 轴比例
- `ax.set_yscale('linear'/'log'/'symlog')` - 设置 y 轴刻度

**刻度线：**
- `ax.set_xticks(positions)` - 设置 x 刻度位置
- `ax.set_xticklabels(labels)` - 设置 x 刻度标签
- `ax.tick_params(axis='both', labelsize=10)` - 自定义刻度外观

**网格和脊柱：**
- `ax.grid(True, alpha=0.3, linestyle='--')` - 添加网格
- `ax.spines['top'].set_visible(False)` - 隐藏顶部书脊
- `ax.spines['right'].set_visible(False)` - 隐藏右书脊

**图例：**
- `ax.legend(loc='best', fontsize=10, frameon=True)` - 添加图例
- `ax.legend(handles, labels)` - 自定义图例

**外观和布局：**
- `ax.set_aspect('equal'/'auto'/ratio)` - 设置宽高比
- `ax.invert_xaxis()` - 反转 x 轴
- `ax.invert_yaxis()` - 反转 y 轴

### pyplot 模块

用于快速绘图的高级界面。

**人物创作：**
- `plt.figure()` - 创建新图形
- `plt.subplots()` - 创建图形和轴
- `plt.subplot()` - 将子图添加到当前图形

**绘图（使用当前轴）：**
- `plt.plot()` - 线图
- `plt.scatter()` - 散点图
- `plt.bar()` - 条形图
- `plt.hist()` - 直方图
-（所有轴方法均可用）

**显示并保存：**
- `plt.show()` - 显示图形
- `plt.savefig()` - 保存图形
- `plt.close()` - 关闭图

**款式：**
- `plt.style.use(style_name)` - 应用样式表
- `plt.style.available` - 列出可用样式

**状态管理：**
- `plt.gca()` - 获取当前轴
- `plt.gcf()` - 获取当前数字
- `plt.sca(ax)` - 设置当前轴
- `plt.clf()` - 清除当前图形
- `plt.cla()` - 清除当前轴

## 线条和标记样式

### 线条样式
- `'-'` 或 `'solid'` - 实线
- `'--'` 或 `'dashed'` - 虚线
- `'-.'` 或 `'dashdot'` - 点划线
- `':'` 或 `'dotted'` - 虚线
- `''` 或 `' '` 或 `'None'` - 无行

### 标记样式
- `'.'` - 点标记
- `'o'` - 圆形标记
- `'v'`、`'^'`、`'<'`、`'>'` - 三角形标记
- `'s'` - 方形标记
- `'p'` - 五边形标记
- `'*'` - 星形标记
- `'h'`、`'H'` - 六边形标记
- `'+'` - 加号标记
- `'x'` - X 标记
- `'D'`、`'d'` - 菱形标记

### 颜色规格
**单字符快捷键：**
- `'b'` - 蓝色
- `'g'` - 绿色
- `'r'` - 红色
- `'c'` - 青色
- `'m'` - 洋红色
- `'y'` - 黄色
- `'k'` - 黑色
- `'w'` - 白色

**命名颜色：**
- `'steelblue'`、`'coral'`、`'teal'` 等
- 查看完整列表：https://matplotlib.org/stable/gallery/color/named_colors.html

**其他格式：**
- 十六进制：`'#FF5733'`
- RGB 元组：`(0.1, 0.2, 0.3)`
- RGBA 元组：`(0.1, 0.2, 0.3, 0.5)`

## 常用参数

### 绘制函数参数

<<<代码块_2>>>

### 分散函数参数

<<<代码块_3>>>

### 文本参数

<<<代码块_4>>>

## rcParams 配置

全局自定义的常用 rcParams 设置：

<<<代码块_5>>>

## 复杂布局的 GridSpec

<<<代码块_6>>>

## 3D 绘图

```python
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot types
ax.plot(x, y, z)                    # 3D line
ax.scatter(x, y, z)                 # 3D scatter
ax.plot_surface(X, Y, Z)            # 3D surface
ax.plot_wireframe(X, Y, Z)          # 3D wireframe
ax.contour(X, Y, Z)                 # 3D contour
ax.bar3d(x, y, z, dx, dy, dz)       # 3D bar

# Customization
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.view_init(elev=30, azim=45)      # Set viewing angle
```

## 动画

```python
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()
line, = ax.plot([], [])

def init():
    ax.set_xlim(0, 2*np.pi)
    ax.set_ylim(-1, 1)
    return line,

def update(frame):
    x = np.linspace(0, 2*np.pi, 100)
    y = np.sin(x + frame/10)
    line.set_data(x, y)
    return line,

anim = FuncAnimation(fig, update, init_func=init,
                     frames=100, interval=50, blit=True)

# Save animation
anim.save('animation.gif', writer='pillow', fps=20)
anim.save('animation.mp4', writer='ffmpeg', fps=20)
```

## 图像操作

```python
# Read and display image
img = plt.imread('image.png')
ax.imshow(img)

# Display matrix as image
ax.imshow(matrix, cmap='viridis', aspect='auto',
          interpolation='nearest', origin='lower')

# Colorbar
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('Values')

# Image extent (set coordinates)
ax.imshow(img, extent=[x_min, x_max, y_min, y_max])
```

## 事件处理

```python
# Mouse click event
def on_click(event):
    if event.inaxes:
        print(f'Clicked at x={event.xdata:.2f}, y={event.ydata:.2f}')

fig.canvas.mpl_connect('button_press_event', on_click)

# Key press event
def on_key(event):
    print(f'Key pressed: {event.key}')

fig.canvas.mpl_connect('key_press_event', on_key)
```

## 有用的实用程序

```python
# Get current axis limits
xlims = ax.get_xlim()
ylims = ax.get_ylim()

# Set equal aspect ratio
ax.set_aspect('equal', adjustable='box')

# Share axes between subplots
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

# Twin axes (two y-axes)
ax2 = ax1.twinx()

# Remove tick labels
ax.set_xticklabels([])
ax.set_yticklabels([])

# Scientific notation
ax.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))

# Date formatting
import matplotlib.dates as mdates
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=7))
```