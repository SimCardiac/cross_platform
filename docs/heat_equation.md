# Heat Equation Solver

## 问题描述

求解二维热传导方程：

$$\frac{\partial u}{\partial t} = \alpha \Delta u + f$$

其中：
- $u(x, y, t)$ 是温度场
- $\alpha$ 是热扩散系数
- $f(x, y, t)$ 是源项
- 边界条件：$u = 0$ on $\partial\Omega$ (齐次 Dirichlet 边界条件)
- 初始条件：$u(x, y, 0) = \sin(\pi x) \sin(\pi y)$

## 精确解

使用分离变量法，在无源项 ($f = 0$) 的情况下，精确解为：

$$u_{\text{exact}}(x, y, t) = e^{-2\pi^2\alpha t} \sin(\pi x) \sin(\pi y)$$

## 数值方法

### 时间离散化

使用**隐式欧拉法**（向后欧拉）：

$$\frac{u^{n+1} - u^n}{\Delta t} = \alpha \Delta u^{n+1} + f^{n+1}$$

重新整理为 Helmholtz 方程：

$$(I - \alpha \Delta t \Delta) u^{n+1} = u^n + \Delta t f^{n+1}$$

或者：

$$u^{n+1} - \alpha \Delta t \Delta u^{n+1} = u^n + \Delta t f^{n+1}$$

### 空间离散化

- 使用 **DMStag** 的单元中心网格（cell-centered grid）
- 五点差分格式离散化拉普拉斯算子
- 虚拟单元法处理 Dirichlet 边界条件：$u_{\text{ghost}} = -u_{\text{interior}}$

### 线性系统求解

每个时间步需要求解 Helmholtz 方程：

$$(I + \alpha \Delta t L) u^{n+1} = u^n + \Delta t f$$

其中 $L$ 是离散拉普拉斯算子矩阵。

使用 **红黑 Gauss-Seidel 迭代法**进行 matrix-free 求解：

$$u_{\text{new}} = \frac{u^n + \Delta t f + \alpha \Delta t \cdot \frac{1}{h^2}(u_l + u_r + u_d + u_u)}{1 + 4\alpha \Delta t / h^2}$$

## 编译和运行

### 编译

```bash
cd build
cmake ..
make ex_heat_stagger
```

### 运行示例

基本运行（默认参数）：
```bash
mpirun -np 4 ./ex_heat_stagger
```

自定义参数：
```bash
mpirun -np 4 ./ex_heat_stagger \
  -nx 64 -ny 64 \        # 网格大小
  -dt 0.001 \             # 时间步长
  -T 0.1 \                # 最终时间
  -heat_check_error       # 计算误差
```

### 命令行参数

- `-nx N`: x方向网格点数（默认：64）
- `-ny N`: y方向网格点数（默认：64）
- `-dt DT`: 时间步长（默认：0.001）
- `-T T_final`: 最终时间（默认：0.1）
- `-heat_check_error`: 计算与精确解的 L2 误差

## 收敛性测试

### 空间-时间收敛性

```bash
# 测试不同网格大小的收敛性
for n in 16 32 64 128; do
  dt=$(echo "scale=6; 0.1/$n/$n" | bc)
  echo "Grid: ${n}x${n}, dt=$dt"
  mpirun -np 4 ./ex_heat_stagger -nx $n -ny $n -dt $dt -T 0.05 -heat_check_error
done
```

### 预期结果

对于二阶精度方法（空间和时间都是二阶），误差应该满足：

$$\|u - u_{\text{exact}}\|_{L^2} = O(h^2 + \Delta t)$$

当 $\Delta t \propto h^2$ 时，误差为 $O(h^2)$，即网格加密2倍，误差减小4倍。

### 实际测试结果

```
n=16,  dt=.000390, ||u - u_exact||_L2 = 0.000161477
n=32,  dt=.000097, ||u - u_exact||_L2 = 4.02136e-05    (误差比: 4.01)
n=64,  dt=.000024, ||u - u_exact||_L2 = 1.00404e-05    (误差比: 4.00)
n=128, dt=.000006, ||u - u_exact||_L2 = 2.5111e-06     (误差比: 4.00)
```

**结论**：达到了预期的二阶收敛性！✅

## 代码结构

### 主要函数

1. **`SetupInitialCondition`**: 设置初始条件 $u(x,y,0) = \sin(\pi x)\sin(\pi y)$

2. **`SetupRHS`**: 设置源项 $f(x,y,t)$（本例中为0）

3. **`GaussSeidelSweep`**: 执行一次红黑 Gauss-Seidel 扫描
   - 求解 $(I + \alpha \Delta t L) u^{n+1} = u^n + \Delta t f$
   - 使用虚拟单元法处理边界条件

4. **`ComputeResidualNorm`**: 计算残差范数
   - 残差：$r = (u^n + \Delta t f) - (I + \alpha \Delta t L) u^{n+1}$

5. **`ComputeL2Error`**: 计算与精确解的 L2 误差
   - $\|u - u_{\text{exact}}\|_{L^2} = \sqrt{h_x h_y} \|u - u_{\text{exact}}\|_{\ell^2}$

### 时间演化循环

```cpp
while (t < T_final) {
  t += dt;
  VecCopy(u, uOld);                    // 保存上一步解
  SetupRHS(dm, f, fLocal, t);          // 更新源项
  
  // Gauss-Seidel 迭代求解 Helmholtz 方程
  for (its = 1; its <= maxIts; ++its) {
    GaussSeidelSweep(...);
    ComputeResidualNorm(...);
    if (resNorm <= tol) break;
  }
}
```

## 稳定性分析

隐式欧拉法是**无条件稳定**的，不受 CFL 条件限制。

对于显式方法（如向前欧拉），稳定性条件为：

$$\Delta t \leq \frac{h^2}{4\alpha}$$

而隐式方法允许使用更大的时间步长，但每步需要求解线性系统。

## 扩展方向

1. **添加源项**：修改 `rhs_f` 函数实现非齐次问题
2. **不同边界条件**：Neumann、Robin 边界条件
3. **非均匀扩散系数**：$\alpha = \alpha(x, y)$
4. **输出可视化**：集成 VTI 输出功能
5. **预条件器**：使用 PETSc KSP 求解器替代 Gauss-Seidel
6. **自适应时间步长**：根据误差估计调整 $\Delta t$

## 参考文献

1. LeVeque, R. J. (2007). *Finite Difference Methods for Ordinary and Partial Differential Equations*.
2. PETSc Documentation: https://petsc.org/release/docs/
3. DMStag Guide: https://petsc.org/release/docs/manualpages/DMSTAG/
