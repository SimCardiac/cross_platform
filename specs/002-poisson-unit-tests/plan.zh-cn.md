# 实施计划：泊松方程单元测试

**分支**: `002-poisson-unit-tests` | **日期**: 2025-10-23 | **规范**: /Users/pengfei/Documents/GitHub/cross_platform/specs/002-poisson-unit-tests/spec.md
**输入**: 来自 `/specs/002-poisson-unit-tests/spec.md` 的功能规范

**注意**: 此模板由 `/speckit.plan` 命令填写。有关执行工作流程，请参见 `.specify/templates/commands/plan.md`。

## 摘要

基于现有实现（`src/poisson_DMStag_2D.cpp`）交付一个专注的泊松方程单元测试套件。
测试通过制造解验证离散算子恒等式、边界条件执行、求解器残差和网格细化收敛性。
支持多个右侧/精确解对，并通过命令行标志接受所有参数。
构建首先通过源 `~/.bashrc` 运行，然后使用 CMake + Make（源外 `build/`）。

## 技术背景

<!--
  需要操作：用项目的技术细节替换本节内容。
  此处的结构以咨询能力呈现，以指导迭代过程。
-->

**语言/版本**: C++17（如果工具链允许则 C++20）  
**主要依赖**: PETSc（Vec/Mat/DM/KSP）、CMake、Make  
**存储**: 不适用  
**测试**: 需要澄清（CTest + 轻量级 CLI 测试 vs. GoogleTest）  
**目标平台**: macOS、Linux  
**项目类型**: 单项目（CLI + 测试）  
**性能目标**: 泊松测试子集在本地完成 ≤ 60 秒  
**约束**: 确定性测试；无网络；仅 CLI 参数化；构建前源 `~/.bashrc`  
**规模/范围**: 仅单元测试级别（算子、边界条件、残差、细化收敛性）

## 章程检查

*门槛：必须在阶段 0 研究前通过。在阶段 1 设计后重新检查。*

- PETSc API 合规性：通过（泊松测试利用现有的基于 PETSc 的实现）
- 现代 C++ 标准：通过（C++17 最低；鼓励 RAII）
- 模板编程策略：通过/不适用（泊松测试不排除模板；保持可扩展设计）
- 数值正确性：通过（基于 MMS 的收敛、残差阈值、边界条件验证）
- 分阶段交付和覆盖：通过（阶段 1 是泊松；与章程路线图一致）
- 离散化/时间：通过（椭圆；验证离散拉普拉斯恒等式）
- 求解器配置：通过（可配置残差容差和迭代限制）
- 边界条件：通过（Dirichlet/Neumann 验证）
- 质量保证和文档：通过（单元 + 验证测试；添加快速入门/文档）

## 项目结构

### 文档（此功能）

```text
specs/[###-feature]/
├── plan.md              # 此文件（/speckit.plan 命令输出）
├── research.md          # 阶段 0 输出（/speckit.plan 命令）
├── data-model.md        # 阶段 1 输出（/speckit.plan 命令）
├── quickstart.md        # 阶段 1 输出（/speckit.plan 命令）
├── contracts/           # 阶段 1 输出（/speckit.plan 命令）
└── tasks.md             # 阶段 2 输出（/speckit.tasks 命令 - 非 /speckit.plan 创建）
```

### 源代码（仓库根目录）
<!--
  需要操作：用此功能的具体布局替换下面的占位符树。
  删除未使用的选项并使用真实路径扩展选择的结构
  （例如，apps/admin、packages/something）。
  交付的计划不得包含选项标签。
-->

```text
src/
├── poisson_DMStag_2D.cpp           # 现有实现（参考）
└── ...

tests/
└── unit/
  └── poisson/
    ├── test_poisson_cli.sh     # 可选的 CLI 冒烟测试
    └── (CTest 或 gtest 源待定)

specs/002-poisson-unit-tests/
├── spec.md
├── plan.md
├── research.md
├── data-model.md
├── quickstart.md
└── contracts/
  └── cli.md
```

**结构决策**: 单项目结构；测试保存在 `tests/unit/poisson/` 下。
合约表达为 CLI 参数规范（`specs/.../contracts/cli.md`）。

## 复杂性跟踪

> **仅当章程检查有必须证明的违规时填写**

| 违规 | 为何需要 | 拒绝更简单的替代方案的原因 |
|-----------|------------|-------------------------------------|
| [例如，第 4 个项目] | [当前需求] | [为什么 3 个项目不足] |
| [例如，存储库模式] | [具体问题] | [为什么直接数据库访问不足] |

---

## 阶段 0：大纲和研究

要解决的未知问题：
- 测试框架选择：仅 CTest vs. GoogleTest（需要澄清）
- 多个右侧/精确对的确切 CLI 参数架构
- 默认容差值和预期阶数

研究任务：
- 研究 C++ 单元测试的 CTest 集成和 CLI 包装
- 定义 2-3 个制造解，具有泊松的解析 u 和相应的 f
- 建立 CLI 标志：网格大小、案例选择、容差、边界条件类型

输出：research.md 整合决策、理由和替代方案。

## 阶段 1：设计和合约

要生成的工件：
- data-model.md：实体（TestCase、Boundary、Tolerances、Pair）
- contracts/cli.md：参数和多对处理的 CLI 架构
- quickstart.md：构建/运行说明，包括源 `~/.bashrc`、CMake + Make、示例命令
- 代理上下文更新（已执行）

阶段 1 后重新评估章程检查；所有门槛预期保持通过。
