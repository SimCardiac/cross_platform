# CLI 合约：泊松测试运行器

我们将使用单个二进制文件（例如，ex_poisson_center 或小型包装器）通过 CLI 标志驱动测试。通过可重复标志在一次运行中支持多个 (u_exact, f) 对。

## 标志

- --grid nx ny
  - 必需。示例：--grid 32 32
- --domain x0 x1 y0 y1
  - 可选。默认值：0 1 0 1
- --pair ID
  - 可重复。选择制造解对。接受：poly2 | sinpi | cospi
  - 示例：--pair sinpi --pair poly2
- --bc west east south north
  - 可选。默认值：Dirichlet Dirichlet Dirichlet Dirichlet
  - 每个是 Dirichlet|Neumann
- --tol l2 linf
  - 可选。默认值：1e-10 1e-8
- --refine N
  - 可选。如果提供，执行 N 次连续均匀细化并报告观察到的阶数
- --output path
  - 可选。写入摘要 JSON 的位置（每对残差、误差和通过/失败）

## 行为

- 对于每个 --pair 出现，使用指定的网格和边界条件运行完整求解
- 报告每对指标：残差 L2、误差 L2/Linf 和状态
- 如果设置了 --refine，在 2 倍细化上重复求解并计算 log(error) vs log(h) 的斜率
- 退出代码：如果所有对通过则为 0；如果有任何失败则为非零
