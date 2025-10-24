# Unified PDE Solver for Poisson and Heat Equations

这是一个统一的 PDE 求解器，可以通过命令行参数选择求解 **Poisson 方程**或 **Heat 方程**。

## 编译

```bash
cd build
cmake ..
make unified_pde_solver_2D
```

## 使用方法

### 1. Poisson 方程（稳态问题）

求解：`-∇²u = f` in Ω = (0,1)×(0,1)，边界条件 u = 0

**基本用法：**
```bash
./unified_pde_solver_2D -problem_type poisson -nx 32 -ny 32 -check_error
```

**参数说明：**
- `-problem_type poisson`: 选择 Poisson 方程
- `-nx, -ny`: 网格分辨率（默认 32×32）
- `-check_error`: 计算并输出 L2 误差
- `-convergence_test`: 输出可解析的收敛测试数据

**示例：**
```bash
# 不同网格分辨率测试
./unified_pde_solver_2D -problem_type poisson -nx 16 -ny 16 -check_error
./unified_pde_solver_2D -problem_type poisson -nx 32 -ny 32 -check_error
./unified_pde_solver_2D -problem_type poisson -nx 64 -ny 64 -check_error
./unified_pde_solver_2D -problem_type poisson -nx 128 -ny 128 -check_error
```

### 2. Heat 方程（时间依赖问题）

求解：`∂u/∂t - α∇²u = f` in Ω × (0, T]，边界条件 u = 0

**基本用法：**
```bash
./unified_pde_solver_2D -problem_type heat -nx 32 -ny 32 -dt 0.001 -T 0.01 -check_error
```

**参数说明：**
- `-problem_type heat`: 选择 Heat 方程
- `-nx, -ny`: 网格分辨率（默认 32×32）
- `-alpha`: 热扩散系数（默认 0.1）
- `-dt`: 时间步长（默认 0.001）
- `-T`: 最终时间（默认 0.1）
- `-check_error`: 在最终时间计算 L2 误差
- `-convergence_test`: 输出可解析的收敛测试数据

**示例：**
```bash
# 空间收敛测试
./unified_pde_solver_2D -problem_type heat -nx 16 -ny 16 -dt 0.001 -T 0.01 -check_error
./unified_pde_solver_2D -problem_type heat -nx 32 -ny 32 -dt 0.001 -T 0.01 -check_error
./unified_pde_solver_2D -problem_type heat -nx 64 -ny 64 -dt 0.001 -T 0.01 -check_error

# 时间收敛测试
./unified_pde_solver_2D -problem_type heat -nx 64 -ny 64 -dt 0.01 -T 0.1 -check_error
./unified_pde_solver_2D -problem_type heat -nx 64 -ny 64 -dt 0.005 -T 0.1 -check_error
./unified_pde_solver_2D -problem_type heat -nx 64 -ny 64 -dt 0.0025 -T 0.1 -check_error
```

## 实现细节

### 共享基础设施

1. **DMStag 网格管理**：两种方程都使用 PETSc 的 DMStag 进行网格管理
2. **Red-Black Gauss-Seidel 迭代**：共享相同的迭代求解框架
3. **Manufactured Solutions**：从 `manufactured_solutions.h` 库中获取精确解

### 区别

| 特性 | Poisson | Heat |
|------|---------|------|
| DOF 配置 | `(0,0,1)` 单元中心 | `(0,1,0)` 边界中心 |
| 时间依赖性 | 稳态 | 时间依赖 |
| 精确解 | `POISSON::TRIG_2D` | `HEAT::DECAY_2D` |
| 边界条件 | Ghost cell BC | 边界 DOF = 0 |

### 收敛率

- **Poisson**: 空间 O(h²) = 2.0
- **Heat**: 空间 O(h²) = 2.0，时间 O(dt) = 1.0 (隐式 Euler)

## 并行执行

两种模式都支持 MPI 并行：

```bash
mpirun -np 4 ./unified_pde_solver_2D -problem_type poisson -nx 128 -ny 128 -check_error
mpirun -np 4 ./unified_pde_solver_2D -problem_type heat -nx 128 -ny 128 -dt 0.001 -T 0.1 -check_error
```

## 设计思想

这个统一求解器展示了如何：
1. **共享基础设施**：不同 PDE 可以共享相同的数值框架
2. **模块化设计**：通过参数选择不同的问题类型
3. **可扩展性**：未来可以轻松添加更多方程类型（如 Stokes、Navier-Stokes）

通过这种设计，我们避免了代码重复，同时保持了灵活性。
