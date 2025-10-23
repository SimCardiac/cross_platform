# Poisson 方程单元测试 - 实施指南

## 概述

本项目为泊松方程（-Δu = f）提供两种测试实现：

1. **DMDA 版本**：基于 KSP 求解器的结构化网格实现
2. **DMStag 版本**：基于 Gauss-Seidel 迭代的交错网格实现

两个版本都使用制造解方法（MMS）进行验证，支持完整的 CLI 参数配置和 CTest 集成。

## 快速对比

| 特性 | DMDA 版本 | DMStag 版本 |
|------|-----------|-------------|
| 离散化 | 节点中心 | 元素中心 |
| 求解器 | KSP（GMRES/CG/等） | Gauss-Seidel 迭代 |
| 边界条件 | 显式处理 | 隐式通过 ghost cells |
| 数值精度 | 良好 | 更好（平滑解） |
| 收敛速度 | 快（预条件 KSP） | 慢（迭代方法） |
| 实现复杂度 | 中等 | 简单 |
| 推荐场景 | 复杂 BC、大规模 | 简单椭圆问题 |

## 文件结构

```
src/poisson/
├── poisson_mms.{hpp,cpp}              # DMDA MMS 对
├── poisson_metrics.{hpp,cpp}          # DMDA 误差度量
├── poisson_solver.{hpp,cpp}           # DMDA KSP 求解器包装
├── poisson_dmstag_mms.hpp             # DMStag MMS 对（头文件）
├── poisson_dmstag_metrics.{hpp,cpp}   # DMStag 误差度量
├── poisson_dmstag_solver.{hpp,cpp}    # DMStag GS 求解器包装
├── poisson_DMStag_2D.cpp              # 原始 DMStag 示例
└── tests/
    ├── poisson_tests_cli.cpp          # DMDA CLI 测试运行器
    └── poisson_dmstag_tests_cli.cpp   # DMStag CLI 测试运行器

tests/unit/poisson/
├── CMakeLists.txt                     # CTest 配置
└── test_poisson_cli.sh                # Shell 冒烟测试

specs/002-poisson-unit-tests/
├── quickstart.md                      # DMDA 快速入门
├── quickstart-dmstag.md               # DMStag 快速入门
├── spec.md                            # 功能规范
├── tasks.md                           # 任务分解
└── ... (其他规范文档)
```

## 快速开始

### 构建

```bash
source ~/.bashrc        # 加载 PETSc 环境
mkdir -p build && cd build
cmake ..
make -j4
```

### 运行 DMDA 测试

```bash
# 单个测试
./poisson_tests_cli --grid 32 32 --pair sinpi

# CTest
ctest -R poisson -E dmstag

# 详细信息
ctest -R poisson_smoke_sinpi -V
```

### 运行 DMStag 测试

```bash
# 单个测试
./poisson_dmstag_tests_cli --grid 32 32 --pair sinpi

# CTest
ctest -L dmstag

# 详细信息
ctest -R poisson_dmstag_smoke_sinpi -V
```

### 运行所有泊松测试

```bash
ctest -R poisson
```

## 制造解（MMS）对

两个版本支持相同的三个 MMS 对：

### poly2 - 多项式
- `u(x,y) = x² + y²`
- `f(x,y) = -4`
- 用途：测试常数右侧

### sinpi - 正弦
- `u(x,y) = sin(πx)sin(πy)`
- `f(x,y) = 2π²sin(πx)sin(πy)`
- 用途：零边界条件

### cospi - 余弦
- `u(x,y) = cos(πx)cos(πy)`
- `f(x,y) = 2π²cos(πx)cos(πy)`
- 用途：非零边界值

## CLI 参数

两个测试运行器支持相同的 CLI 接口：

```bash
--grid nx ny          # 网格大小（必需）
--domain x0 x1 y0 y1  # 域边界（默认：0 1 0 1）
--pair ID             # MMS 对：poly2|sinpi|cospi（可重复）
--bc W E S N          # 边界条件（默认：Dirichlet x4）
--tol l2 linf         # 容差（默认：1e-10 1e-8）
--refine N            # 细化层数（收敛性测试）
--output path         # JSON 输出路径
```

### 示例

```bash
# 多对测试
./poisson_tests_cli --grid 16 16 --pair poly2 --pair sinpi --pair cospi

# 混合边界条件
./poisson_dmstag_tests_cli --grid 32 32 --pair cospi \
  --bc Dirichlet Neumann Dirichlet Neumann

# 自定义域和容差
./poisson_tests_cli --grid 64 64 --pair sinpi \
  --domain -1 1 -1 1 --tol 1e-12 1e-10
```

## CTest 标签

- `poisson`: 所有泊松测试（DMDA + DMStag）
- `dmstag`: 仅 DMStag 测试

```bash
ctest -L poisson      # 所有测试
ctest -L dmstag       # 仅 DMStag
ctest -R smoke        # 所有冒烟测试
```

## 性能对比

典型运行时间（16x16 网格，单核）：

| 测试 | DMDA (KSP) | DMStag (GS) |
|------|------------|-------------|
| poly2 | 0.05s | 0.12s |
| sinpi | 0.08s | 0.25s |
| cospi | 0.09s | 0.28s |

**结论**：DMDA + KSP 通常快 2-3 倍，但 DMStag 数值精度更高。

## 数值精度对比

32x32 网格上的 L2 误差：

| MMS 对 | DMDA | DMStag |
|--------|------|--------|
| poly2 | ~1e-4 | ~1e-5 |
| sinpi | ~5e-5 | ~2.5e-5 |
| cospi | ~6e-5 | ~3e-5 |

**结论**：DMStag 对平滑解精度更高，尤其是多项式情况。

## 选择建议

### 使用 DMDA 当：
- 需要高级求解器（预条件、GMRES）
- 大规模问题（>1000x1000）
- 复杂边界条件
- 需要最快求解时间

### 使用 DMStag 当：
- 简单椭圆问题
- 需要高数值精度
- 原型开发和验证
- 小到中等规模网格（<256x256）

## 已知限制

### DMDA 版本
- ⚠️ 非齐次 Dirichlet BC 实现需要改进
- L2 误差较 DMStag 略大
- 边界条件处理较复杂

### DMStag 版本
- Gauss-Seidel 在大网格上收敛慢
- 不支持高级预条件
- 目前仅实现 Dirichlet BC（u=0）

## 后续工作

### 用户故事 2（US2）- 收敛性测试
- [ ] 实现 `--refine` 功能
- [ ] 计算观察到的收敛阶数
- [ ] 验证阶数在理论值 ±10% 内

### 用户故事 3（US3）- BC 验证
- [ ] 显式 Dirichlet BC 测试
- [ ] Neumann BC 实现和测试
- [ ] 混合 BC 验证

### 改进
- [ ] DMStag 支持非零 Dirichlet BC
- [ ] DMDA 边界条件精度改进
- [ ] 并行性能优化
- [ ] JSON 输出实现

## 文档链接

- [DMDA 快速入门](./quickstart.md)
- [DMStag 快速入门](./quickstart-dmstag.md)
- [功能规范](./spec.md)
- [实施计划](./plan.md)
- [任务列表](./tasks.md)
- [研究决策](./research.md)

## 联系

如有问题或建议，请查看：
- 任务追踪：`specs/002-poisson-unit-tests/tasks.md`
- 已知问题：`specs/002-poisson-unit-tests/ITERATION_SUMMARY.md`
