# Poisson 方程的相容性条件

## 问题描述

对于带有**纯 Neumann 边界条件**的 Poisson 方程：

$$
\begin{cases}
\nabla^2 p = f & \text{in } \Omega \\
\frac{\partial p}{\partial n} = 0 & \text{on } \partial\Omega
\end{cases}
$$

该问题有解**当且仅当**右端项 $f$ 满足以下**相容性条件（compatibility condition）**：

$$
\int_\Omega f \, dV = 0
$$

## 数学推导

使用散度定理（divergence theorem）：

$$
\int_\Omega \nabla^2 p \, dV = \int_{\partial\Omega} \frac{\partial p}{\partial n} \, dS
$$

由于边界条件 $\frac{\partial p}{\partial n} = 0$，右边为零：

$$
\int_\Omega \nabla^2 p \, dV = 0
$$

将 Poisson 方程 $\nabla^2 p = f$ 代入：

$$
\int_\Omega f \, dV = 0
$$

**结论**：如果 $\int_\Omega f \, dV \neq 0$，方程**无解**！

## 物理意义

在流体力学中，Poisson 方程通常用于压力修正：

$$
\nabla^2 p = \nabla \cdot \mathbf{u}^*
$$

其中 $\mathbf{u}^*$ 是中间速度场。相容性条件意味着：

- **不可压缩流动**：$\nabla \cdot \mathbf{u} = 0$，所以 $\int \nabla \cdot \mathbf{u}^* \, dV = 0$
- 通过求解 Poisson 方程得到压力场，然后修正速度使其满足无散度条件

## 解的唯一性

对于纯 Neumann 边界条件：
- 解**不是唯一的**：如果 $p$ 是解，那么 $p + C$（任意常数）也是解
- 解只在**常数意义下唯一**

为了得到唯一解，需要**额外约束**：
1. 固定一点的值：$p(x_0, y_0) = 0$
2. 强制零均值：$\int_\Omega p \, dV = 0$

## 数值求解问题

### 问题诊断

当使用迭代法（如 Gauss-Seidel）求解纯 Neumann Poisson 方程时：

**现象**：
- 残差不收敛，卡在某个较大的值
- 每次迭代残差几乎不变

**原因**：
1. RHS 不满足相容性条件：$\int f \, dV \neq 0$
2. 常数零空间未被移除

### 解决方案

#### 1. 强制 RHS 满足相容性条件

在每个时间步求解 Poisson 方程前：

```cpp
// 计算 RHS 的积分
PetscScalar divSum;
PetscCall(VecSum(div, &divSum));

// 计算并移除平均值
const PetscScalar divMean = divSum / (Nx * Ny);
PetscCall(VecShift(div, -divMean));
```

#### 2. 移除解的常数零空间

在每次迭代后强制解的零均值：

```cpp
PetscErrorCode PoissonZero(const DM &dm, Vec &u, Vec &uLocal, 
                           PetscInt Nx, PetscInt Ny) {
  // 1. 计算所有 element center DOFs 的平均值
  PetscScalar localSum = 0.0;
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      localSum += aU[ey][ex][icenter];
    }
  }
  
  // 2. 全局归约
  PetscScalar globalSum;
  MPI_Allreduce(&localSum, &globalSum, 1, MPIU_SCALAR, MPI_SUM, comm);
  const PetscScalar average = globalSum / (Nx * Ny);
  
  // 3. 从所有值中减去平均值
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      aU[ey][ex][icenter] -= average;
    }
  }
}
```

#### 3. 迭代过程

```cpp
for (int iter = 1; iter <= maxIter; ++iter) {
  // Gauss-Seidel sweep
  GaussSeidelSweep_Poisson(...);
  
  // 移除常数零空间（关键！）
  PoissonZero(dm, u, uLocal, Nx, Ny);
  
  // 计算残差
  CalculateResidual_Poisson(..., &residualNorm);
  
  if (residualNorm < tolerance) break;
}
```

## 示例：Taylor-Green Vortex

Taylor-Green vortex 是一个满足相容性条件的经典解析解：

### 速度场

$$
\begin{align}
u(x,y,t) &= -\cos(\pi x) \sin(\pi y) e^{-2\pi^2 \nu t} \\
v(x,y,t) &= \sin(\pi x) \cos(\pi y) e^{-2\pi^2 \nu t}
\end{align}
$$

### 散度

$$
\nabla \cdot \mathbf{u} = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0
$$

**散度恒为零**，自动满足相容性条件！

### 压力场

$$
p(x,y,t) = -\frac{1}{4}[\cos(2\pi x) + \cos(2\pi y)] e^{-4\pi^2 \nu t}
$$

压力满足零均值：$\int_\Omega p \, dV = 0$

## 常见错误示例

### 错误 1：使用不满足相容性的 RHS

```cpp
// ❌ 错误：cos(x)cos(y) 的积分不为零
aDiv[ey][ex][ip] = std::cos(x) * std::cos(y);

// 积分：∫∫ cos(x)cos(y) dx dy = [sin(1)]^2 ≈ 0.708 ≠ 0
```

**结果**：迭代不收敛，残差停滞

### 错误 2：忘记调用 PoissonZero

```cpp
for (int iter = 1; iter <= maxIter; ++iter) {
  GaussSeidelSweep_Poisson(...);
  // ❌ 忘记移除零空间
  CalculateResidual_Poisson(...);
}
```

**结果**：即使 RHS 满足相容性条件，迭代仍可能不收敛

### 错误 3：PoissonZero 实现错误

```cpp
// ❌ 错误：使用局部 nx*ny 而非全局 Nx*Ny
PetscInt totalCells = nx * ny;  // 局部大小
PetscScalar average = globalSum / totalCells;  // 错误的平均值
```

**结果**：平均值计算错误，零空间未被正确移除

## 验证方法

### 1. 检查 RHS 的积分

```cpp
PetscScalar sum;
VecSum(rhs, &sum);
if (std::abs(sum) > 1e-12) {
  printf("Warning: RHS integral = %e (should be ~0)\n", sum);
}
```

### 2. 检查解的均值

```cpp
PetscScalar sum;
// 只对 element center DOFs 求和
// ...
PetscScalar average = sum / (Nx * Ny);
if (std::abs(average) > 1e-12) {
  printf("Warning: Solution average = %e (should be ~0)\n", average);
}
```

### 3. 监控残差收敛

```cpp
for (int iter = 1; iter <= maxIter; ++iter) {
  // ... Gauss-Seidel sweep ...
  
  if (iter % 100 == 0) {
    printf("Iteration %d, residual = %e\n", iter, residualNorm);
  }
}
```

**正常收敛**：残差单调递减
**不收敛**：残差停滞或振荡

## 参考文献

1. Griebel, M., Dornseifer, T., & Neunhoeffer, T. (1998). *Numerical Simulation in Fluid Dynamics: A Practical Introduction*. SIAM.

2. Chorin, A. J. (1968). Numerical solution of the Navier-Stokes equations. *Mathematics of Computation*, 22(104), 745-762.

3. Quarteroni, A., & Valli, A. (1994). *Numerical Approximation of Partial Differential Equations*. Springer.

## 总结

- **纯 Neumann Poisson 方程**需要 RHS 满足 $\int f dV = 0$
- 解只在常数意义下唯一，需要额外约束（如零均值）
- 数值求解时必须：
  1. 强制 RHS 零均值
  2. 每次迭代后移除解的零空间
- Taylor-Green vortex 是理想的测试用例（自然满足相容性条件）
