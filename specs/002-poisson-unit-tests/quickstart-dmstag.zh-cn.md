# 快速入门：泊松 DMStag 单元测试

## 概述

DMStag（Staggered Grid）版本的泊松方程测试使用元素中心（cell-centered）离散化，采用 Gauss-Seidel 迭代求解器。与 DMDA 版本相比，DMStag 提供更好的数值精度和边界条件处理。

## 构建

1) 确保环境中有 PETSc 和 MPI
2) 源您的 shell 初始化，使环境变量存在：
   ```bash
   source ~/.bashrc
   ```
3) 配置和构建（源外）：
   ```bash
   mkdir -p build && cd build
   cmake ..
   make -j
   ```

## 运行

### 直接 CLI 执行

- 32x32 网格上的单对：
  ```bash
  ./poisson_dmstag_tests_cli --grid 32 32 --pair sinpi
  ```

- 多对，自定义容差：
  ```bash
  ./poisson_dmstag_tests_cli --grid 16 16 --pair sinpi --pair poly2 --tol 1e-10 1e-8
  ```

- 混合边界条件：
  ```bash
  ./poisson_dmstag_tests_cli --grid 32 32 --pair cospi --bc Dirichlet Neumann Dirichlet Neumann
  ```

- 带细化的收敛性测试：
  ```bash
  ./poisson_dmstag_tests_cli --grid 16 16 --pair sinpi --refine 2
  ```

### 通过 CTest 测试

从构建目录：

- 运行所有泊松测试（包括 DMStag）：
  ```bash
  ctest -R poisson
  ```

- 仅运行 DMStag 测试：
  ```bash
  ctest -L dmstag
  ```

- 详细输出运行：
  ```bash
  ctest -L dmstag -V
  ```

- 并行运行测试：
  ```bash
  ctest -L dmstag -j4
  ```

- 运行特定测试：
  ```bash
  ctest -R poisson_dmstag_smoke_poly2
  ```

## 预期输出

CLI 运行器打印每对指标：
```
=== 泊松 DMStag 测试运行器 ===
网格：32 x 32
域：[0,1] x [0,1]
对：sinpi 

--- 测试：sinpi ---
  L2 误差：   2.51e-05
  Linf 误差： 3.14e-05
  残差：      6.11e-08
  状态：通过

=== 摘要 ===
总计：1 | 通过：1 | 失败：0
```

CTest 输出：
```
Test project /path/to/build
    Start 1: poisson_dmstag_smoke_poly2
1/3 Test #1: poisson_dmstag_smoke_poly2 .......   Passed    0.45 sec
    Start 2: poisson_dmstag_smoke_sinpi
2/3 Test #2: poisson_dmstag_smoke_sinpi .......   Passed    0.58 sec
    Start 3: poisson_dmstag_multiple_pairs
3/3 Test #3: poisson_dmstag_multiple_pairs ....   Passed    1.12 sec

100% tests passed, 0 tests failed out of 3
```

## DMStag vs DMDA

### DMStag 优势
- **更好的边界条件处理**：元素中心离散化避免边界处的特殊处理
- **更高的数值精度**：对于平滑解，通常达到更好的误差阶数
- **更简洁的实现**：Gauss-Seidel 迭代无需复杂的边界修正

### DMDA 优势
- **KSP 求解器支持**：可以使用 PETSc 的高级求解器（GMRES、CG 等）
- **更灵活的 BC 选项**：支持更多类型的边界条件
- **成熟的生态系统**：更多示例和文档

### 选择建议
- 对于简单的椭圆问题，**推荐 DMStag**
- 需要高级求解器或复杂 BC 时，使用 DMDA
- 性能关键应用，测试两者并选择更快的

## 制造解（MMS）对

### poly2 - 多项式解
- 精确解：`u(x,y) = x² + y²`
- 右侧：`f(x,y) = -4`
- 特点：常数右侧，适合测试基本算子正确性

### sinpi - 三角函数解
- 精确解：`u(x,y) = sin(πx)sin(πy)`
- 右侧：`f(x,y) = 2π²sin(πx)sin(πy)`
- 特点：零边界条件，经典测试案例

### cospi - 余弦解
- 精确解：`u(x,y) = cos(πx)cos(πy)`
- 右侧：`f(x,y) = 2π²cos(πx)cos(πy)`
- 特点：非零边界值，测试 BC 实现

## 注意事项

- 所有参数必须通过命令行提供；无配置文件
- 退出代码：如果所有测试通过则为 0，如果有任何失败则为非零
- 测试是确定性的；重新运行在容差内产生相同的结果
- DMStag 使用 Gauss-Seidel 迭代，收敛速度对网格大小敏感
- 对于大型网格（>128x128），考虑增加 `max_its` 或降低 `tol`

## 故障排除

### 测试失败：L2 误差过大
- 检查网格是否足够细（至少 16x16）
- 验证容差设置合理（默认 1e-10 for L2）
- 确认 Gauss-Seidel 收敛（残差应 < 1e-8）

### 测试超时
- 减小网格大小或增加超时限制
- DMStag 的 Gauss-Seidel 在大网格上可能较慢
- 考虑使用 MPI 并行：`mpirun -np 4 ./poisson_dmstag_tests_cli ...`

### 残差不收敛
- 检查右侧函数是否正确
- 验证边界条件与 MMS 对一致
- 对于纯 Neumann BC，确保兼容性条件（∫f = 0）
