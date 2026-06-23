# Cross-Platform High-Performance PDE Solver Framework

> 基于 PETSc + Kokkos 的跨平台结构网格 PDE 求解器，目标平台为 Apple Silicon (M2) macOS，未来可扩展到 GPU 加速。

---

## 项目概述

本项目构建了一套跨平台高性能科学计算软件栈，以 **Kokkos** 作为并行计算后端，**PETSc** 作为分布式计算基础设施，最终目标是实现一个通用的、基于交错网格（Staggered Grid/DMStag）的 PDE 求解框架。

当前阶段主要完成以下 PDE 的数值求解与收敛性验证：

| 方程类型 | 状态 | 说明 |
|---------|------|------|
| **Poisson 方程** | ✅ 完成 | 稳态椭圆型，二阶中心差分，已验证二阶收敛 |
| **热传导方程** | ✅ 完成 | 非稳态抛物型，隐式 Euler + 五点差分，二阶收敛 |
| **Stokes 方程** | ✅ 完成 | 稳态/非稳态鞍点问题，含投影法 |
| **Navier-Stokes** | 🔄 开发中 | 非线性对流项 |
| **非线性系统** | ✅ 完成 | JFNK (Jacobian-Free Newton-Krylov) 框架 + 弹簧质点链动力学 |

---

## 技术栈

```
┌─────────────────────────────────────────┐
│          应用层: PDE Solvers            │
├─────────────────────────────────────────┤
│  PETSc (KSP/SNES/DMStag)               │
├─────────────────────────────────────────┤
│  Kokkos + OpenMP  │  MPI (OpenMPI)      │
├─────────────────────────────────────────┤
│  GCC 14.2.0       │  OpenBLAS 0.3.27   │
├─────────────────────────────────────────┤
│  Apple Silicon (M2) + macOS             │
└─────────────────────────────────────────┘
```

### 依赖管理 (Spack)

```yaml
spack:
  specs:
  - openmpi@5.0.8
  - openblas@0.3.27
  - kokkos@4.6.01 +openmp
  - petsc@3.23.4 +kokkos +openmp
  - py-pip
  view: true
  concretizer:
    unify: true
```

### 编译器选择

- **主编译器**: GCC 14.2.0（通过 Spack 安装）
- **宿主编译器**: Homebrew 安装的基础 GCC（用于自举编译）
- macOS 自带 Clang 不支持 Fortran，无法编译 PETSc

### 并行策略

- **节点间**: MPI 分布式并行（多台 Mac 设备）
- **节点内**: OpenMP 共享内存并行（通过 Kokkos 后端）
- **未来**: Kokkos Metal 后端（等待上游支持）

---

## 项目结构

```
cross_platform/
├── CMakeLists.txt              # 构建系统
├── spack.yaml                  # Spack 环境配置
├── readme.md                   # 原始开发笔记
│
├── src/
│   ├── analytical/
│   │   └── unsteady.h          # 解析解/制造解（Taylor-Green 涡等）
│   │
│   ├── poisson/
│   │   ├── poisson_DMStag_2D.cpp       # Poisson 方程求解器
│   │   ├── empty.cpp                    # 库桩文件
│   │   ├── tests/
│   │   │   ├── poisson_tests_cli.cpp        # Poisson 单元测试
│   │   │   └── poisson_dmstag_tests_cli.cpp # DMStag 边界测试
│   │   └── rubbish/                 # 历史遗留代码
│   │
│   ├── heat/
│   │   ├── heat_DMStag_2D.cpp       # 热方程 (Dirichlet BC)
│   │   ├── heat_DMStag_CN_2D.cpp    # 热方程 (Crank-Nicolson)
│   │   └── heat_DMStag_edge_2D.cpp  # 热方程 (Neumann BC)
│   │
│   ├── stokes/
│   │   ├── stokes_DMStag_2D.cpp            # Stokes 稳态
│   │   ├── stokes_DMStag_projection_2D.cpp # 投影法
│   │   └── stokes_DMStag_unsteady_2D.cpp   # Stokes 非稳态
│   │
│   ├── nonlinear/
│   │   ├── test_JFNK.cpp       # JFNK 方法测试
│   │   └── beam_system.cpp     # 弹簧-质点链动力学
│   │
│   ├── unified_pde_solver_2D.cpp  # 统一 PDE 求解器入口
│   ├── DMStag_boundary_helpers.h   # DMStag 边界辅助函数
│   ├── DMStag_boundary_helpers.cpp
│   └── io/
│       ├── vti.cpp / vti.hpp       # VTK ImageData 输出
│       └── a.h
│
├── tests/unit/
│   ├── poisson/                # Poisson 收敛性测试脚本
│   ├── heat/                   # 热方程收敛性测试脚本
│   └── stokes/                 # Stokes 收敛性测试脚本
│
├── tools/
│   └── unsteady.py             # 解析解代码生成器
│
└── docs/
    ├── heat_equation.md               # 热方程数学描述
    ├── heat_convergence_results.md    # 热方程收敛结果
    ├── poisson_convergence_results.md # Poisson 收敛结果
    ├── poisson_compatibility_condition.md
    ├── poisson_neumann_bc.md
    └── installation.md
```

---

## 物理模型与数值方法

### 1. Poisson 方程

**问题**: $-\Delta u = f$，齐次 Dirichlet 边界条件

**数值方法**:
- 空间离散: DMStag 单元中心五点差分（二阶精度）
- 边界处理: 虚拟单元法 ($u_{\text{ghost}} = -u_{\text{interior}}$)
- 线性求解: KSP (GMRES / 红黑 Gauss-Seidel)

**验证结果**: 使用制造解 $u = \sin(\pi x)\sin(\pi y)$，确认达到二阶收敛。

### 2. 热传导方程

**问题**: $\frac{\partial u}{\partial t} = \alpha \Delta u + f$

**数值方法**:
- 时间离散: 隐式 Euler（一阶）/ Crank-Nicolson（二阶）
- 空间离散: 五点差分（二阶）
- 每个时间步求解 Helmholtz 方程 $(I - \alpha\Delta t\Delta)u^{n+1} = u^n + \Delta t f^{n+1}$

**验证结果**: 使用 $\Delta t \propto h^2$ 的耦合收敛测试，达到二阶收敛。

### 3. Stokes 方程

**问题**: 
$$\begin{cases} -\Delta \mathbf{u} + \nabla p = \mathbf{f} \\ \nabla \cdot \mathbf{u} = 0 \end{cases}$$

**数值方法**:
- 交错网格 (MAC grid): 速度在面心，压力在单元中心
- 时间推进: 投影法 (Projection Method)
- 压力 Poisson 方程使用 KSP 求解

### 4. 非线性系统 (JFNK)

**Jacobian-Free Newton-Krylov** 方法:
- 外层 Newton 迭代求解 $F(u) = 0$
- 内层 GMRES 求解线性化系统，无需显式构造 Jacobian 矩阵
- 使用有限差分近似 Jacobian-vector 乘积

### 5. 弹簧-质点链动力学

模拟竖直悬挂的 N 个弹簧-质点系统的二维运动，用于验证非线性动力学求解能力。

---

## 关键设计模式

### DMStag 边界辅助 (`DMStag_boundary_helpers.h`)

提供安全的邻居访问函数，自动处理 Dirichlet 边界条件：

```cpp
// 左邻居：若 ex==0 则返回 -u_center（反射边界实现 u=0）
PetscScalar DMStag_GetLeft(arr, ex, ey, slot, Nx);
PetscScalar DMStag_GetRight(arr, ex, ey, slot, Nx);
PetscScalar DMStag_GetDown(arr, ex, ey, slot, Ny);
PetscScalar DMStag_GetUp(arr, ex, ey, slot, Ny);
```

### 解析解生成器 (`tools/unsteady.py`)

通过 Python 脚本自动生成 C++ 解析解头文件 (`unsteady.h`)，包含：
- Taylor-Green 2D 涡的精确解
- 各方程的制造解 (Manufactured Solution) 源项
- 支持 Poisson、Heat、Stokes、Navier-Stokes 方程

---

## 编译与运行

### 环境准备

```bash
# 1. 安装宿主编译器
brew install gcc

# 2. 通过 Spack 安装完整环境
spack install

# 3. 生成解析解代码
python tools/unsteady.py
```

### 编译

```bash
mkdir -p build && cd build
cmake ..
make -j$(nproc)
```

### 运行示例

```bash
# Poisson 方程
./poisson_DMStag_2D -poisson_check_error

# 热方程
./heat_DMStag_2D -heat_check_error

# Stokes 方程
./stokes_DMStag_2D -stokes_check_error

# 多进程 MPI 运行
mpirun -np 2 ./poisson_DMStag_2D -ksp_monitor_short -ksp_converged_reason
```

---

## 收敛性验证

所有求解器均已通过 MMS (Method of Manufactured Solutions) 验证：

| 方程 | 空间精度 | 时间精度 | 验证状态 |
|------|---------|---------|---------|
| Poisson | 二阶 | — | ✅ |
| Heat (隐式 Euler) | 二阶 | 一阶 | ✅ |
| Heat (Crank-Nicolson) | 二阶 | 二阶 | ✅ |
| Stokes (稳态) | 二阶 | — | ✅ |

详见 `docs/` 目录下的各收敛性结果文档。

---

## 未来规划

1. **GPU 加速**: 等待 Kokkos 支持 Metal 后端，或通过 MoltenVK 间接支持
2. **Navier-Stokes**: 添加非线性对流项，构建完整的不可压缩流动求解器
3. **三维扩展**: 将求解器扩展到 3D 结构网格
4. **自适应网格**: 集成 PETSc 的 DMForest (p4est) 实现 AMR
5. **多物理场耦合**: 热-流-固耦合求解
