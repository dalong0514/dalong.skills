<!-- 此文件由机器翻译自 code-generation-printing.md -->

# SymPy 代码生成和打印

本文档介绍了 SymPy 生成各种语言的可执行代码、将表达式转换为不同输出格式以及自定义打印行为的功能。

## 代码生成

### 转换为 NumPy 函数

```python
from sympy import symbols, sin, cos, lambdify
import numpy as np

x, y = symbols('x y')
expr = sin(x) + cos(y)

# Create NumPy function
f = lambdify((x, y), expr, 'numpy')

# Use with NumPy arrays
x_vals = np.linspace(0, 2*np.pi, 100)
y_vals = np.linspace(0, 2*np.pi, 100)
result = f(x_vals, y_vals)
```

### Lambdify 选项

<<<代码块_1>>>

### 生成 C/C++ 代码

<<<代码块_2>>>

### 生成 Fortran 代码

<<<代码块_3>>>

### 高级代码生成

<<<代码块_4>>>

### 条码打印机

<<<代码块_5>>>

## 打印和输出格式

### 漂亮的印刷

<<<代码块_6>>>

### LaTeX 输出

```python
from sympy import latex, symbols, Integral, sin, sqrt

x, y = symbols('x y')
expr = Integral(sin(x)**2, (x, 0, pi))

# Convert to LaTeX
latex_str = latex(expr)
print(latex_str)
# \int\limits_{0}^{\pi} \sin^{2}{\left(x \right)}\, dx

# Custom LaTeX formatting
latex_str = latex(expr, mode='equation')  # Wrapped in equation environment
latex_str = latex(expr, mode='inline')    # Inline math

# For matrices
from sympy import Matrix
M = Matrix([[1, 2], [3, 4]])
latex(M)  # \left[\begin{matrix}1 & 2\\3 & 4\end{matrix}\right]
```

### MathML 输出

```python
from sympy.printing.mathml import mathml, print_mathml
from sympy import sin, pi

expr = sin(pi/4)

# Content MathML
mathml_str = mathml(expr)

# Presentation MathML
mathml_str = mathml(expr, printer='presentation')

# Print to console
print_mathml(expr)
```

### 字符串表示

```python
from sympy import symbols, sin, pi, srepr, sstr

x = symbols('x')
expr = sin(x)**2

# Standard string (what you see in Python)
str(expr)  # 'sin(x)**2'

# String representation (prettier)
sstr(expr)  # 'sin(x)**2'

# Reproducible representation
srepr(expr)  # "Pow(sin(Symbol('x')), Integer(2))"
# This can be eval()'ed to recreate the expression
```

### 定制印刷

```python
from sympy.printing.str import StrPrinter

class MyPrinter(StrPrinter):
    def _print_Symbol(self, expr):
        return f"<{expr.name}>"

    def _print_Add(self, expr):
        return " PLUS ".join(self._print(arg) for arg in expr.args)

printer = MyPrinter()
x, y = symbols('x y')
print(printer.doprint(x + y))  # "<x> PLUS <y>"
```

## Python 代码生成

### autowrap - 编译和导入

```python
from sympy.utilities.autowrap import autowrap
from sympy import symbols

x, y = symbols('x y')
expr = x**2 + y**2

# Automatically compile C code and create Python wrapper
f = autowrap(expr, backend='cython')
# or backend='f2py' for Fortran

# Use like a regular function
result = f(3, 4)  # 25
```

### ufuncify - 创建 NumPy ufunc

```python
from sympy.utilities.autowrap import ufuncify
import numpy as np

x, y = symbols('x y')
expr = x**2 + y**2

# Create universal function
f = ufuncify((x, y), expr)

# Works with NumPy broadcasting
x_arr = np.array([1, 2, 3])
y_arr = np.array([4, 5, 6])
result = f(x_arr, y_arr)  # [17, 29, 45]
```

## 表达式树操作

### 行走表达式树

```python
from sympy import symbols, sin, cos, preorder_traversal, postorder_traversal

x, y = symbols('x y')
expr = sin(x) + cos(y)

# Preorder traversal (parent before children)
for arg in preorder_traversal(expr):
    print(arg)

# Postorder traversal (children before parent)
for arg in postorder_traversal(expr):
    print(arg)

# Get all subexpressions
subexprs = list(preorder_traversal(expr))
```

### 树中的表达式替换

```python
from sympy import Wild, symbols, sin, cos

x, y = symbols('x y')
a = Wild('a')

expr = sin(x) + cos(y)

# Pattern matching and replacement
new_expr = expr.replace(sin(a), a**2)  # sin(x) -> x**2
```

## Jupyter 笔记本集成

### 显示数学

```python
from sympy import init_printing, display
from IPython.display import display as ipy_display

# Initialize printing for Jupyter
init_printing(use_latex='mathjax')  # or 'png', 'svg'

# Display expressions beautifully
expr = Integral(sin(x)**2, x)
display(expr)  # Renders as LaTeX in notebook

# Multiple outputs
ipy_display(expr1, expr2, expr3)
```

### 互动小部件

```python
from sympy import symbols, sin
from IPython.display import display
from ipywidgets import interact, FloatSlider
import matplotlib.pyplot as plt
import numpy as np

x = symbols('x')
expr = sin(x)

@interact(a=FloatSlider(min=0, max=10, step=0.1, value=1))
def plot_expr(a):
    f = lambdify(x, a * expr, 'numpy')
    x_vals = np.linspace(-np.pi, np.pi, 100)
    plt.plot(x_vals, f(x_vals))
    plt.show()
```

## 表示形式之间的转换

### 字符串到 SymPy

```python
from sympy.parsing.sympy_parser import parse_expr
from sympy import symbols

x, y = symbols('x y')

# Parse string to expression
expr = parse_expr('x**2 + 2*x + 1')
expr = parse_expr('sin(x) + cos(y)')

# With transformations
from sympy.parsing.sympy_parser import (
    standard_transformations,
    implicit_multiplication_application
)

transformations = standard_transformations + (implicit_multiplication_application,)
expr = parse_expr('2x', transformations=transformations)  # Treats '2x' as 2*x
```

### LaTeX 到 SymPy

```python
from sympy.parsing.latex import parse_latex

# Parse LaTeX
expr = parse_latex(r'\frac{x^2}{y}')
# Returns: x**2/y

expr = parse_latex(r'\int_0^\pi \sin(x) dx')
```

### Mathematica 到 SymPy

```python
from sympy.parsing.mathematica import parse_mathematica

# Parse Mathematica code
expr = parse_mathematica('Sin[x]^2 + Cos[y]^2')
# Returns SymPy expression
```

## 导出结果

### 导出到文件

```python
from sympy import symbols, sin
import json

x = symbols('x')
expr = sin(x)**2

# Export as LaTeX to file
with open('output.tex', 'w') as f:
    f.write(latex(expr))

# Export as string
with open('output.txt', 'w') as f:
    f.write(str(expr))

# Export as Python code
with open('output.py', 'w') as f:
    f.write(f"from numpy import sin\n")
    f.write(f"def f(x):\n")
    f.write(f"    return {lambdify(x, expr, 'numpy')}\n")
```

### Pickle SymPy 对象

```python
import pickle
from sympy import symbols, sin

x = symbols('x')
expr = sin(x)**2 + x

# Save
with open('expr.pkl', 'wb') as f:
    pickle.dump(expr, f)

# Load
with open('expr.pkl', 'rb') as f:
    loaded_expr = pickle.load(f)
```

## 数值评估和精度

### 高精度评估

```python
from sympy import symbols, pi, sqrt, E, exp, sin
from mpmath import mp

x = symbols('x')

# Standard precision
pi.evalf()  # 3.14159265358979

# High precision (1000 digits)
pi.evalf(1000)

# Set global precision with mpmath
mp.dps = 50  # 50 decimal places
expr = exp(pi * sqrt(163))
float(expr.evalf())

# For expressions
result = (sqrt(2) + sqrt(3)).evalf(100)
```

### 数字替换

```python
from sympy import symbols, sin, cos

x, y = symbols('x y')
expr = sin(x) + cos(y)

# Numerical evaluation
result = expr.evalf(subs={x: 1.5, y: 2.3})

# With units
from sympy.physics.units import meter, second
distance = 100 * meter
time = 10 * second
speed = distance / time
speed.evalf()
```

## 常见模式

### 模式 1：生成并执行代码

```python
from sympy import symbols, lambdify
import numpy as np

# 1. Define symbolic expression
x, y = symbols('x y')
expr = x**2 + y**2

# 2. Generate function
f = lambdify((x, y), expr, 'numpy')

# 3. Execute with numerical data
data_x = np.random.rand(1000)
data_y = np.random.rand(1000)
results = f(data_x, data_y)
```

### 模式 2：创建 LaTeX 文档

```python
from sympy import symbols, Integral, latex
from sympy.abc import x

# Define mathematical content
expr = Integral(x**2, (x, 0, 1))
result = expr.doit()

# Generate LaTeX document
latex_doc = f"""
\\documentclass{{article}}
\\usepackage{{amsmath}}
\\begin{{document}}

We compute the integral:
\\begin{{equation}}
{latex(expr)} = {latex(result)}
\\end{{equation}}

\\end{{document}}
"""

with open('document.tex', 'w') as f:
    f.write(latex_doc)
```

### 模式3：交互式计算

```python
from sympy import symbols, simplify, expand
from sympy.parsing.sympy_parser import parse_expr

x, y = symbols('x y')

# Interactive input
user_input = input("Enter expression: ")
expr = parse_expr(user_input)

# Process
simplified = simplify(expr)
expanded = expand(expr)

# Display
print(f"Simplified: {simplified}")
print(f"Expanded: {expanded}")
print(f"LaTeX: {latex(expr)}")
```

### 模式 4：批量代码生成

```python
from sympy import symbols, lambdify
from sympy.utilities.codegen import codegen

# Multiple functions
x = symbols('x')
functions = {
    'f1': x**2,
    'f2': x**3,
    'f3': x**4
}

# Generate C code for all
for name, expr in functions.items():
    [(c_name, c_code), _] = codegen((name, expr), 'C')
    with open(f'{name}.c', 'w') as f:
        f.write(c_code)
```

### 模式 5：性能优化

```python
from sympy import symbols, sin, cos, cse
import numpy as np

x, y = symbols('x y')

# Complex expression with repeated subexpressions
expr = sin(x + y)**2 + cos(x + y)**2 + sin(x + y)

# Common subexpression elimination
replacements, reduced = cse(expr)
# replacements: [(x0, sin(x + y)), (x1, cos(x + y))]
# reduced: [x0**2 + x1**2 + x0]

# Generate optimized code
for var, subexpr in replacements:
    print(f"{var} = {subexpr}")
print(f"result = {reduced[0]}")
```

## 重要提示

1. **NumPy 兼容性：** 将 `lambdify` 与 NumPy 一起使用时，请确保您的表达式使用 NumPy 中可用的函数。

2. **性能：** 对于数值计算，请始终在循环中使用 `lambdify` 或代码生成，而不是 `subs()` + `evalf()`。

3. **精度：** 在需要时使用 `mpmath` 进行任意精度算术。

4. **代码生成注意事项：** 生成的代码可能无法处理所有边缘情况。彻底测试。

5. **编译：** `autowrap` 和 `ufuncify` 需要 C/Fortran 编译器，并且可能需要在您的系统上进行配置。

6. **解析：** 解析用户输入时，进行验证和清理以避免代码注入漏洞。

7. **Jupyter：** 为了在 Jupyter 笔记本中获得最佳结果，请在会话开始时调用 `init_printing()`。