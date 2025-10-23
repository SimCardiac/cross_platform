# 快速入门：泊松单元测试

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
  ./poisson_tests_cli --grid 32 32 --pair sinpi
  ```

- 多对，自定义容差：
  ```bash
  ./poisson_tests_cli --grid 16 16 --pair sinpi --pair poly2 --tol 1e-10 1e-8
  ```

- 混合边界条件：
  ```bash
  ./poisson_tests_cli --grid 32 32 --pair cospi --bc Dirichlet Neumann Dirichlet Neumann
  ```

- 带细化的收敛性测试：
  ```bash
  ./poisson_tests_cli --grid 16 16 --pair sinpi --refine 2
  ```

### 通过 CTest 测试

从构建目录：

- 运行所有泊松测试：
  ```bash
  ctest -R poisson
  ```

- 详细输出运行：
  ```bash
  ctest -R poisson -V
  ```

- 并行运行测试：
  ```bash
  ctest -R poisson -j4
  ```

- 运行特定测试：
  ```bash
  ctest -R poisson_smoke_poly2
  ```

## 预期输出

CLI 运行器打印每对指标：
```
=== 泊松测试运行器 ===
网格：32 x 32
域：[0,1] x [0,1]
对：sinpi 

--- 测试：sinpi ---
  L2 误差：   1.23e-05
  Linf 误差： 4.56e-06
  残差：      2.34e-11
  状态：通过

=== 摘要 ===
总计：1 | 通过：1 | 失败：0
```

CTest 输出：
```
Test project /path/to/build
    Start 1: poisson_smoke_poly2
1/4 Test #1: poisson_smoke_poly2 ..............   Passed    0.52 sec
    Start 2: poisson_smoke_sinpi
2/4 Test #2: poisson_smoke_sinpi ..............   Passed    0.61 sec
...
100% tests passed, 0 tests failed out of 4
```

## 注意事项

- 所有参数必须通过命令行提供；无配置文件
- 退出代码：如果所有测试通过则为 0，如果有任何失败则为非零
- 测试是确定性的；重新运行在容差内产生相同的结果
