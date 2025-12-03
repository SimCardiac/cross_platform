正在使用 "替换文件中的字符串"


```

```
在 PETSc 中，设置参数有两种主要方式：
1.  **命令行参数（推荐）**：灵活性最高，无需重新编译。
2.  **代码硬编码**：用于设置默认值或强制特定配置。

我已经修改了 test_JFNK.cpp，在代码中显式添加了参数设置部分。

### 1. 代码中设置参数 (已添加)

我在 `SNESSetFromOptions` 之前添加了以下配置：

*   **Newton 迭代参数 (`SNESSetTolerances`)**:
    *   `abstol = 1e-8`: 绝对误差容限 ($||F(x)|| < 10^{-8}$)
    *   `rtol = 1e-8`: 相对误差容限
    *   `maxit = 50`: 最大非线性迭代步数

*   **线性求解器类型 (`KSPSetType`)**:
    *   强制使用 `KSPGMRES`。

*   **GMRES 参数 (`KSPGMRESSetRestart`)**:
    *   设置重启参数为 30 (即 GMRES(30))。

*   **线性迭代参数 (`KSPSetTolerances`)**:
    *   `rtol = 1e-5`: 线性方程组求解的相对误差。在 JFNK 中，这个值决定了牛顿方向的准确度。

*   **预条件子 (`PCSetType`)**:
    *   显式设置为 `PCNONE`，因为纯 Matrix-Free 方法没有矩阵数据来构建 ILU 等预条件子。

### 2. 命令行参数 (运行时覆盖)

即使在代码中写了上述设置，您仍然可以通过命令行参数覆盖它们（因为代码最后调用了 `SNESSetFromOptions`）。

重新编译并运行：
```bash
cd build
make test_JFNK
```

**常用命令行参数示例：**

*   **改变 Newton 收敛标准：**
    ```bash
    ./test_JFNK -snes_rtol 1e-10 -snes_max_it 100
    ```

*   **改变线性求解器 (KSP) 参数：**
    ```bash
    ./test_JFNK -ksp_type gmres -ksp_gmres_restart 100 -ksp_rtol 1e-3
    ```

*   **查看详细收敛过程 (非常有用)：**
    ```bash
    ./test_JFNK -snes_monitor -ksp_monitor
    ```
    *   `-snes_monitor`: 显示每步 Newton 迭代的残差。
    *   `-ksp_monitor`: 显示每步 Newton 内部，线性求解器 (GMRES) 的收敛过程。

*   **查看实际使用的所有参数：**
    ```bash
    ./test_JFNK -snes_view
    ```
    这会打印出最终生效的所有求解器配置，方便检查参数是否设置成功。

已进行更改。