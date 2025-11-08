# Analytical Solutions Generator

本目录包含用于自动生成制造解（manufactured solutions）的工具。

## 文件说明

- `stokes.h` - Stokes 方程的制造解头文件（自动生成，**请勿手动编辑**）
- `stokes.py` - Python 脚本，用于从符号定义自动推导并生成 `stokes.h`

## 使用方法

### 生成 stokes.h

使用默认参数（Taylor–Green 涡，ν=1.0）：

```bash
cd /path/to/cross_platform
python src/analytical/stokes.py
```

或使用虚拟环境中的 Python（如果已配置）：

```bash
.venv/bin/python src/analytical/stokes.py
```

### 自定义参数

指定粘性系数：

```bash
python src/analytical/stokes.py --nu 0.25
```

指定输出路径：

```bash
python src/analytical/stokes.py --out src/analytical/custom_stokes.h
```

查看所有选项：

```bash
python src/analytical/stokes.py --help
```

## 工作原理

脚本使用 SymPy 符号计算库：

1. **定义制造解**：指定速度场 (u, v) 和压力场 p 的符号表达式
2. **自动推导强迫项**：计算 fx = ∂u/∂t - ν∇²u + ∂p/∂x，fy = ∂v/∂t - ν∇²v + ∂p/∂y
3. **生成 C++ 代码**：将符号表达式转换为 C++ 函数

## 依赖

需要安装 SymPy：

```bash
pip install sympy
```

或使用项目虚拟环境：

```bash
python -m venv .venv
source .venv/bin/activate  # macOS/Linux
# 或 .venv\Scripts\activate  # Windows
pip install sympy
```

## 当前支持的解

### Taylor–Green 涡（2D，时间依赖）

制造解定义：
- u(x,y,t) = -cos(πx) sin(πy) exp(-2π²νt)
- v(x,y,t) = sin(πx) cos(πy) exp(-2π²νt)
- p(x,y,t) = -¼(cos(2πx) + cos(2πy)) exp(-4π²νt)

自动推导的强迫项：
- fx(x,y,t) = ½π exp(-4π²t) sin(2πx)
- fy(x,y,t) = ½π exp(-4π²t) sin(2πy)

特性：
- 满足不可压条件：∇·u = 0
- 解析时间衰减
- 周期边界条件友好

## 扩展指南

要添加新的制造解，编辑 `stokes.py`：

1. 在 `taylor_green_2d()` 函数附近添加新函数，定义 u, v, p
2. 在 `main()` 中添加对应的命令行选项
3. `derive_forcing()` 会自动计算 fx, fy

示例（添加自定义解）：

```python
def custom_solution_2d(nu: float) -> Tuple[sp.Expr, sp.Expr, sp.Expr]:
    x, y, t = sp.symbols('x y t', real=True)
    pi = sp.pi
    
    u = # 你的 u(x,y,t) 表达式
    v = # 你的 v(x,y,t) 表达式
    p = # 你的 p(x,y,t) 表达式
    
    return (u, v, p)
```

然后在 `main()` 中添加：

```python
ap.add_argument('--solution', choices=['taylor_green_2d', 'custom_solution_2d'])
```

## 编译验证

生成头文件后，测试编译：

```bash
cd build
make stokes_DMStag_2D
./stokes_DMStag_2D -nx 16 -ny 16 -dt 1e-6 -T 1e-6 -heat_check_error
```

## 注意事项

- ⚠️ `stokes.h` 是自动生成的，**请勿手动编辑**
- ⚠️ 修改制造解时，请编辑 `stokes.py` 并重新运行脚本
- ✅ 推导过程是符号化的，数学上精确无误
- ✅ 生成的 C++ 代码使用 `std::` 前缀的数学函数
