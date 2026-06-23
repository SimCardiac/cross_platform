# 在 Apple Silicon macOS 上搭建 Kokkos 与 PETSc 异构计算环境

Apple 芯片应该是目前最强CPU芯片，GPU芯片也在快速进步，缺点是核心数量太少，内存太小，但是未来肯定会改善的。
我正在测试一套跨平台高性能科学计算软件栈，该方案以 Kokkos 作为并行计算后端，PETSc 作为分布式计算工具，最终目标是实现一个基于结构网格的求解器。目前，我正基于搭载 M2 芯片的 macOS 设备进行环境搭建与测试。

1. 关于编译器选择
在实际编译过程中遇到了一些问题。最初考虑使用 macOS 自带的 Clang 编译器，但由于其不支持 Fortran 代码编译，导致无法安装 PETSc。目前是统一使用 GCC 工具链， 通过 spack 安装gcc@14.2.0 本身，然后在进行编译。

2. 是否可以直接通过 Spack 编译安装 GCC？
不可行。我目前的作法是先通过 brew install gcc 安装一个基础版本的 GCC，再将其作为宿主编译器，进一步编译安装所需版本的 GCC。

3. 安装 PETSc 时常见的问题
安装过程中确实会遇到依赖项问题，尤其是 OpenBLAS 经常安装失败。目前我所使用的稳定版本是 openblas@0.3.27。

4. 安装 Kokkos 是否会出现问题？
截至目前，在安装 Kokkos 过程中尚未遇到明显问题。

5. 在 macOS 上如何实现 Kokkos + PETSc 的并行计算？
多台 Mac 设备可通过 MPI 结合 OpenMP 实现跨节点的分布式并行与节点内的多层次并行。需要注意的是，必须确保 Kokkos 与 PETSc 在编译时均启用对 OpenMP 的支持，否则无法安装支持 Kokkos 的 PETSc。

6. 如何利用 GPU 进行计算？
正在开发，等待 Kokkos 支持 metal 加速。


see [installation.md](docs/installation.md) for instructions.

# 交错网格

```bash
time mpirun -np 2 ./ex2 -ksp_monitor_short -ksp_converged_reason 
```

```bash
make ex_poisson_stagger  && ./ex_poisson_stagger -poisson_check_error
```

# 输出数据

以 vti 格式 ( vtk 格式的一种) 输出结果，具体格式参考:

https://www.pengfeima.cn/cdn/other/vtk_file_format.pdf

这里需要注意的是 cell data 还是 用 point date, ghost layer 应该怎么加，这些问题我还不是很明白。


# 测试
```bash
ctest -R heat_dmstag_edge_convergence -V
ctest -R heat_dmstag_convergence -V
ctest -R heat_dmstag_CN_temporal_convergence -V
ctest -R heat_dmstag_temporal_convergence -V
```

Implement the following solvers in order. Every solver must converge to a steady state with verified second-order accuracy.

1. Vertex-centered Poisson
2. Cell-centered Poisson
3. Staggered Poisson
4. Vertex-centered heat equation
5. Cell-centered heat equation
6. Staggered heat equation
7. Projection method for unsteady stokes equations
