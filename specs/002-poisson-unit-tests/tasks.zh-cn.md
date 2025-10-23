---

description: "泊松方程单元测试功能的任务列表"
---

# 任务：泊松方程单元测试

**输入**: 来自 `/specs/002-poisson-unit-tests/` 的设计文档
**先决条件**: plan.md（必需）、spec.md（用户故事必需）、research.md、data-model.md、contracts/

**测试**: 此功能明确请求单元测试；每个用户故事包含测试任务。

**组织**: 任务按用户故事分组，以实现每个故事的独立实现和测试。

## 格式：`[ID] [P?] [Story] 描述`

- **[P]**：可以并行运行（不同文件，无依赖）
- **[Story]**：此任务属于哪个用户故事（例如，US1、US2、US3）
- 描述中包含确切的文件路径

## 路径约定

- 单项目：仓库根目录的 `src/`、`tests/`
- 测试组织在 `tests/unit/poisson/` 下

---

## 阶段 1：设置（共享基础设施）

目的：准备所有故事所需的测试基础设施和项目脚手架。

- [x] T001 在 `CMakeLists.txt` 的顶级 CMake 中启用 CTest（添加 `enable_testing()` 和基本配置）
- [x] T002 [P] 创建测试目录结构 `tests/unit/poisson/` 并添加 `.gitkeep`
- [x] T003 [P] 添加 `tests/unit/poisson/CMakeLists.txt` 以通过 `add_test(...)` 注册基于 CLI 的测试
- [x] T004 从根 `CMakeLists.txt` 添加 `tests/unit/poisson/CMakeLists.txt` 的顶级包含

---

## 阶段 2：基础（阻塞先决条件）

目的：在任何用户故事工作之前必须完成的核心构建块。

- [x] T005 创建制造解模块 `src/poisson_mms.hpp` 和 `src/poisson_mms.cpp`（对：poly2、sinpi、cospi）
- [x] T006 [P] 创建误差范数实用程序 `src/poisson_metrics.hpp` 和 `src/poisson_metrics.cpp`（L2/Linf、残差）
- [x] T007 根据 `specs/002-poisson-unit-tests/contracts/cli.md` 实现 CLI 测试运行器 `src/poisson_tests_cli.cpp`
- [x] T008 在运行器 `src/poisson_tests_cli.cpp` 中连接 PETSc/KSP 设置（容差、最大迭代、残差捕获）
- [x] T009 重构或包装 `src/poisson_DMStag_2D.cpp` 使用：公开可调用的求解路径 `run_poisson_case(params)`（新函数）
- [x] T010 [P] 在运行器中为纯 Neumann 添加零空间处理（PETSc `MatNullSpace`）在 `src/poisson_tests_cli.cpp`
- [ ] T011 [P] 在 `src/poisson_tests_cli.cpp` 中添加可选的 JSON 摘要输出 `--output`
- [x] T012 在 CTest `tests/unit/poisson/CMakeLists.txt` 中注册冒烟测试，使用 `--grid 8 8 --pair poly2` 调用运行器

检查点：基础就绪 — 用户故事可以继续。

---

## 阶段 3：用户故事 1 - 本地运行泊松测试（优先级：P1）🎯 MVP

目标：单个命令仅运行泊松测试并打印简洁摘要和正确的退出代码。

独立测试：从 `build/`，`ctest -R poisson -j` 仅执行泊松测试并显示通过/失败计数；退出代码反映失败。

### 用户故事 1 的测试（已请求）

- [x] T013 [P] [US1] 在 `tests/unit/poisson/CMakeLists.txt` 中添加 CTest 标签 `poisson` 和测试定义
- [x] T014 [P] [US1] 创建 CLI 冒烟测试脚本 `tests/unit/poisson/test_poisson_cli.sh`（使用一对调用运行器）

### 用户故事 1 的实现

- [x] T015 [US1] 确保运行器在 `src/poisson_tests_cli.cpp` 中打印每个案例的摘要和最终聚合状态
- [x] T016 [US1] 在 `specs/002-poisson-unit-tests/quickstart.md` 中记录运行命令（ctest 和直接运行器）
- [x] T017 [US1] 添加 CTest 条目以运行 shell 冒烟测试 `tests/unit/poisson/test_poisson_cli.sh`

检查点：US1 独立可运行和可验证。

---

## 阶段 4：用户故事 2 - 验证算子一致性和收敛性（优先级：P1）

目标：使用制造解验证离散算子恒等式和细化时的预期收敛性。

独立测试：在多个网格上运行收敛测试产生观察阶数在预期的 10% 以内；恒等式测试满足误差容差。

### 用户故事 2 的测试（已请求）

- [ ] T018 [P] [US2] 在 `tests/unit/poisson/CMakeLists.txt` 中添加恒等式测试 CTest 条目（例如，常数右侧）调用运行器
- [ ] T019 [P] [US2] 在 `tests/unit/poisson/CMakeLists.txt` 中添加细化测试 CTest 条目（例如，sinpi 对，网格 16/32/64）

### 用户故事 2 的实现

- [ ] T020 [P] [US2] 在 `src/poisson_tests_cli.cpp` 中实现恒等式检查路径（报告相对误差 ≤ 1e-8）
- [ ] T021 [US2] 在 `src/poisson_tests_cli.cpp` 中实现细化循环 `--refine N` 并计算 log(error) vs log(h) 的斜率
- [ ] T022 [US2] 将 `src/poisson_mms.cpp` 中的 MMS 对使用集成到求解/评估管道中
- [ ] T023 [US2] 在 `src/poisson_tests_cli.cpp` 中确保残差阈值检查（≤ 1e-10）和快速失败行为

检查点：US2 独立验证一致性和收敛性。

---

## 阶段 5：用户故事 3 - 边界条件验证（优先级：P2）

目标：确认 Dirichlet/Neumann 边界条件在面和角上正确执行。

独立测试：仅边界的夹具在容差内将值/通量与预期进行比较；处理纯 Neumann 唯一性。

### 用户故事 3 的测试（已请求）

- [ ] T024 [P] [US3] 在 `tests/unit/poisson/CMakeLists.txt` 中添加所有边 Dirichlet 测试，使用 `--bc Dirichlet Dirichlet Dirichlet Dirichlet` 调用运行器
- [ ] T025 [P] [US3] 在 `tests/unit/poisson/CMakeLists.txt` 中添加混合 Neumann/Dirichlet 测试，使用 `--bc Dirichlet Neumann Dirichlet Neumann` 调用运行器

### 用户故事 3 的实现

- [ ] T026 [P] [US3] 在 `src/poisson_tests_cli.cpp` 中添加边界条件规范解析和应用
- [ ] T027 [US3] 在 `src/poisson_tests_cli.cpp` 中实现边界值/通量检查与 MMS（≤ 1e-8 相对）
- [ ] T028 [US3] 在 `src/poisson_tests_cli.cpp` 中确保 Neumann 案例的零空间感知验证（比较到常数）

检查点：US3 独立验证边界条件执行。

---

## 阶段 N：润色和跨领域关注点

目的：跨故事的收尾、文档和健壮性改进。

- [x] T029 [P] 在 `specs/002-poisson-unit-tests/contracts/cli.md` 中为 CLI 标志添加开发人员文档（确保与实现同步）
- [ ] T030 在 `src/poisson/tests/poisson_tests_cli.cpp`、`src/poisson/poisson_mms.*`、`src/poisson/poisson_metrics.*` 中进行代码清理和注释
- [ ] T031 [P] 在 `tests/unit/poisson/CMakeLists.txt` 中添加并行 CTest 运行配置和超时
- [ ] T032 通过全新构建和运行验证 `specs/002-poisson-unit-tests/quickstart.md`

**代码组织（已完成）：**
- [x] 重组 `src/` 结构：创建 `src/poisson/` 模块
- [x] 将测试基础设施移至 `src/poisson/tests/`
- [x] 更新 CMakeLists.txt 路径
- [x] 创建 `src/README.md` 和 `PROJECT_STRUCTURE.md` 文档

---

## 依赖关系和执行顺序

### 阶段依赖关系

- 设置（阶段 1）：无依赖 — 可以立即开始
- 基础（阶段 2）：依赖设置完成 — 阻塞所有用户故事
- 用户故事（阶段 3+）：依赖基础完成；US1 和 US2 可以在阶段 2 后并行进行；US3 也可以并行运行，但受益于 US2 实用程序
- 润色（最后阶段）：依赖选定的用户故事完成

### 用户故事依赖关系

- 用户故事 1（P1）：阶段 2 后开始；独立
- 用户故事 2（P1）：阶段 2 后开始；独立于 US1（共享基础）
- 用户故事 3（P2）：阶段 2 后开始；独立于 US1；可能重用 US2 实用程序但不是严格要求

### 每个用户故事内

- 测试（包括）应在完整实现之前编写并针对预期失败进行验证
- 实现 CLI 和计算路径，然后连接到 CTest
- 确保每个故事的标准可以独立验证

### 并行机会

- 阶段 1 中的 T002、T003 可以并行运行
- 在阶段 2 中：T006、T010、T011 可以并行运行；T012 可以在 T007 准备好后进行
- 用户故事：US1、US2、US3 可以在阶段 2 完成后由不同的贡献者并行开发

---

## 并行示例：用户故事 2

```
# 在运行器存在后并行化恒等式和细化测试连接：
任务："添加恒等式测试 CTest 条目 ..."（T018）
任务："添加细化测试 CTest 条目 ..."（T019）

# 同时实现计算路径：
任务："实现恒等式检查路径 ..."（T020）
任务："实现细化循环 ..."（T021）
```

---

## 实施策略

### MVP 优先（仅用户故事 1）

1. 完成阶段 1：设置
2. 完成阶段 2：基础（阻塞所有故事）
3. 完成阶段 3：用户故事 1
4. 停止并验证：`ctest -R poisson` 返回摘要和正确的退出代码

### 增量交付

1. 设置 + 基础 → 基础就绪
2. 添加 US1 → 独立测试 → 演示
3. 添加 US2 → 独立测试 → 演示
4. 添加 US3 → 独立测试 → 演示

---

# 报告

- 生成的文件：`/Users/pengfei/Documents/GitHub/cross_platform/specs/002-poisson-unit-tests/tasks.md`
- 总任务数：32
- 每个用户故事的任务数：
  - US1：3 个实现 + 2 个测试 = 5（T013–T017）
  - US2：4 个实现 + 2 个测试 = 6（T018–T023）
  - US3：3 个实现 + 2 个测试 = 5（T024–T028）
  - 设置 + 基础 + 润色：16（T001–T012、T029–T032）
- 确定的并行机会：阶段 1（T002、T003）、阶段 2（T006、T010、T011）、US2 测试（T018–T019）、US1 测试（T013–T014）、US3 测试（T024–T025）
- 独立测试标准：
  - US1：`ctest -R poisson` 打印摘要；退出状态反映失败
  - US2：恒等式误差阈值；收敛率在预期的 10% 以内
  - US3：边界值/通量在容差内；处理 Neumann 唯一性
- 建议的 MVP 范围：仅用户故事 1（US1）
- 格式验证：所有任务遵循所需的检查列表格式（复选框 + TaskID + 可选 [P] + 可选 [US#] + 文件路径）
