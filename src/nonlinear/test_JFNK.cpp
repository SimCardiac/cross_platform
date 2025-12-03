static char help[] = "JFNK Example: Solving a simple nonlinear system using Jacobian-Free Newton-Krylov.\n\n";

#include <petscsnes.h>

/* 
   ============================================================================
   1. 黑盒子函数 F(x) (Black Box Function)
   ============================================================================
   问题: 1D Bratu 问题 (强非线性 PDE 离散)
   方程: -u'' - lambda * exp(u) = 0,  x in [0, 1]
   边界: u(0) = u(1) = 0
   离散: (2u_i - u_{i-1} - u_{i+1})/h^2 - lambda * exp(u_i) = 0
   
   参数: lambda 控制非线性强度。
        lambda < 3.51: 有两个解
        lambda = 3.51: 临界点 (Turning point)
        lambda > 3.51: 无解 (对于实数域)
   
   这里我们设置 lambda = 3.0 (接近临界点，非线性较强)。
*/
PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *ctx)
{
  const PetscScalar *xx;
  PetscScalar       *ff;
  PetscInt          i, n;
  PetscReal         h, lambda = 3.0;

  PetscFunctionBeginUser;
  PetscCall(VecGetSize(x, &n));
  h = 1.0 / (n + 1);

  PetscCall(VecGetArrayRead(x, &xx));
  PetscCall(VecGetArray(f, &ff));

  for (i = 0; i < n; i++) {
    PetscScalar u = xx[i];
    PetscScalar u_left  = (i == 0)     ? 0.0 : xx[i-1];
    PetscScalar u_right = (i == n - 1) ? 0.0 : xx[i+1];
    
    // F_i = -u'' - lambda*exp(u)
    // 标准离散形式: (2u - u_L - u_R)/h^2 - lambda * exp(u)
    ff[i] = (2.0 * u - u_left - u_right) / (h * h) - lambda * std::exp(u);
  }

  PetscCall(VecRestoreArrayRead(x, &xx));
  PetscCall(VecRestoreArray(f, &ff));
  PetscFunctionReturn(0);
}

/*
   监控器：用于观察收敛过程
*/
PetscErrorCode Monitor(SNES snes, PetscInt its, PetscReal fnorm, void *ctx)
{
  PetscFunctionBeginUser;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "SNES Iteration %D, Residual Norm: %g\n", its, (double)fnorm));
  PetscFunctionReturn(0);
}

int main(int argc, char **argv)
{
  SNES      snes;         // 非线性求解器上下文
  Vec       x, r;         // 解向量, 残差向量
  Mat       J;            // 雅可比矩阵 (这里将是 Matrix-Free 的壳矩阵)
  PetscInt  size = 10;   // 问题规模: 100 个未知量
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, help));

  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "========================================\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   JFNK Example: 1D Bratu Problem       \n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "   Unknowns: %D, Lambda: 3.0 (Strong)   \n", size));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "========================================\n"));

  // -------------------------------------------------------------------
  // 1. 设置向量
  // -------------------------------------------------------------------
  PetscCall(VecCreate(PETSC_COMM_WORLD, &x));
  PetscCall(VecSetSizes(x, PETSC_DECIDE, size));
  PetscCall(VecSetFromOptions(x));
  PetscCall(VecDuplicate(x, &r));

  // 初始猜测: 设为 0.0 (Bratu 问题通常从 0 开始收敛到下分支解)
  PetscCall(VecSet(x, 0.0));

  // -------------------------------------------------------------------
  // 2. 设置 SNES (非线性求解器)
  // -------------------------------------------------------------------
  PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));
  
  // 设置 "黑盒子" 函数 F(x)
  PetscCall(SNESSetFunction(snes, r, FormFunction, NULL));

  // -------------------------------------------------------------------
  // 3. 配置 JFNK (Jacobian-Free Newton-Krylov)
  // -------------------------------------------------------------------
  /*
     JFNK 核心概念:
     牛顿法需要求解线性系统 J * dx = -F(x)。
     Krylov 子空间方法 (如 GMRES) 求解线性系统时，不需要显式的矩阵 J，
     只需要矩阵-向量乘积 (Matrix-Vector Product) 的操作：y = J * v。
     
     JFNK 利用有限差分来近似这个乘积:
     
     J * v ≈ [F(x + h*v) - F(x)] / h

     其中 h 是一个微小的扰动参数。
     
     在 PETSc 中，我们使用 MatCreateSNESMF 创建一个 "Shell Matrix" (壳矩阵)，
     它不存储任何矩阵元素，而是通过调用 FormFunction 来动态计算上述差分。
  */
  
  // 创建 Matrix-Free 雅可比矩阵上下文
  PetscCall(MatCreateSNESMF(snes, &J));
  
  // 将这个 MF 矩阵设置为 SNES 的雅可比矩阵
  // 参数 2 (J): 用于矩阵-向量乘积操作 (MatMult) 的矩阵 -> 这里是 MF 矩阵
  // 参数 3 (J): 用于构建预条件子 (Preconditioner) 的矩阵
  //             注意：如果我们在这里也传 MF 矩阵，通常意味着我们无法构建复杂的预条件子
  //             (因为没有矩阵元素)，只能用 None 或 User-Defined PC。
  //             在实际大规模问题中，通常会在这里传入一个近似的、稀疏的雅可比矩阵用于预条件 (Preconditioning Matrix)。
  PetscCall(SNESSetJacobian(snes, J, J, MatMFFDComputeJacobian, NULL));

  // -------------------------------------------------------------------
  // 4. 自定义 JFNK 参数 (参数解释)
  // -------------------------------------------------------------------
  /*
     差分参数 'h' (Differencing Parameter) 是 JFNK 精度和稳定性的关键。
     h 太小 -> 浮点数舍入误差 (Cancellation Error)
     h 太大 -> 截断误差 (Truncation Error, 导数近似不准)
     
     PETSc 默认使用 "wp" (Walker and Pernice) 算法来动态计算最优的 h:
     h = error_rel * sqrt(1 + ||x||) / ||v||
     
     MatMFFDSetFunctionError 设置的是 F(x) 计算的相对误差估计 (error_rel)。
     默认约为 1e-8 (双精度机器精度平方根)。如果你的 F(x) 包含噪声，应该调大这个值。
  */
  PetscCall(MatMFFDSetFunctionError(J, 1e-6)); 

  // -------------------------------------------------------------------
  // 4.5 设置求解器参数 (Newton & GMRES)
  // -------------------------------------------------------------------
  // 注意：PETSc 推荐优先使用命令行参数设置这些选项 (例如 -snes_rtol 1e-8)，
  // 但这里展示如何在代码中硬编码默认值。

  // A. 设置非线性求解器 (Newton) 参数
  // 参数: snes, abstol, rtol, stol, maxit, maxf
  // abstol: 绝对收敛误差 (F_norm < abstol)
  // rtol:   相对收敛误差 (F_norm / F_initial < rtol)
  // stol:   步长收敛误差 (dx < stol)
  // maxit:  最大非线性迭代次数
  // maxf:   最大函数评估次数
  // 注意: JFNK 中，每次线性迭代都会调用一次函数计算。
  // 因此 maxf 必须设置得足够大 (大于 maxit_newton * maxit_linear)
  PetscCall(SNESSetTolerances(snes, 1e-8, 1e-8, 1e-8, 500, 10000));

  // A.1 设置求解器类型 (Line Search vs Trust Region)
  // 默认是 SNESNEWTONLS (Newton Line Search)
  // 可选: SNESNEWTONTR (Newton Trust Region)
  // PetscCall(SNESSetType(snes, SNESNEWTONLS)); // 显式设置为线搜索
  // PetscCall(SNESSetType(snes, SNESNEWTONTR)); // 显式设置为可信域

  // A.2 配置线搜索 (如果使用 SNESNEWTONLS)
  // 获取 Line Search 上下文
  SNESLineSearch linesearch;
  PetscCall(SNESGetLineSearch(snes, &linesearch));
  
  // 设置线搜索类型:
  // SNESLINESEARCHBT (Backtracking, 默认) - 适合大多数情况
  // SNESLINESEARCHCP (Critical Point)
  // SNESLINESEARCHL2 (Secant search)
  PetscCall(SNESLineSearchSetType(linesearch, SNESLINESEARCHBT));
  
  // 设置线搜索参数 (alpha, damping, maxstep)
  // alpha: 下降充分性条件参数 (通常 1e-4)
  // maxstep: 最大步长
  // damping: 初始步长阻尼因子 (默认 1.0, 即全牛顿步)
  // PetscCall(SNESLineSearchSetDamping(linesearch, 1.0)); 

  // B. 获取线性求解器上下文 (KSP)
  // JFNK 中，Newton 的每一步都需要解一个线性系统 J*dx = -F
  KSP ksp;
  PetscCall(SNESGetKSP(snes, &ksp));

  // C. 设置线性求解器类型为 GMRES
  PetscCall(KSPSetType(ksp, KSPGMRES));

  // D. 设置 GMRES 重启参数 (Restart)
  // GMRES(30) -> 每 30 次迭代重启一次，节省内存
  PetscCall(KSPGMRESSetRestart(ksp, 30));

  // E. 设置线性求解器收敛参数
  // 参数: ksp, rtol, abstol, dtol, maxits
  // rtol:   相对误差 (相对于线性系统的初始残差) -> JFNK 中通常需要 "过解" (Oversolving) 以保证 Newton 方向准确
  //         但在非精确牛顿法 (Inexact Newton) 中，这个值可以动态调整 (Eisenstat-Walker forcing term)
  // abstol: 绝对误差
  // dtol:   发散判定
  // maxits: 最大线性迭代次数
  PetscCall(KSPSetTolerances(ksp, 1e-8, 1e-8, PETSC_DEFAULT, 100000000));

  // F. 设置 PC (预条件子) 为 None
  // 因为我们使用的是 Matrix-Free 且没有提供近似矩阵 Pmat，
  // 所以默认无法构建 ILU 等预条件子。这里显式设为 NONE。
  // 如果提供了 Pmat，可以使用 PCJACOBI, PCILU 等。
  PC pc;
  PetscCall(KSPGetPC(ksp, &pc));
  PetscCall(PCSetType(pc, PCNONE));

  // 设置监控器
  PetscCall(SNESMonitorSet(snes, Monitor, NULL, NULL));

  // 允许从命令行覆盖设置 (例如 -snes_monitor -ksp_monitor -snes_mf_operator)
  PetscCall(SNESSetFromOptions(snes));

  // -------------------------------------------------------------------
  // 5. 求解
  // -------------------------------------------------------------------
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Starting Solver...\n"));
  PetscCall(SNESSolve(snes, NULL, x));

  // -------------------------------------------------------------------
  // 6. 检查结果
  // -------------------------------------------------------------------
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\nSolution Vector:\n"));
  PetscCall(VecView(x, PETSC_VIEWER_STDOUT_WORLD));

  // 清理资源
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&r));
  PetscCall(MatDestroy(&J));
  PetscCall(SNESDestroy(&snes));
  PetscCall(PetscFinalize());
  return 0;
}


// make test_JFNK && ./test_JFNK -snes_monitor -snes_converged_reason -ksp_monitor -snes_max_funcs 200000
// make test_JFNK && ./test_JFNK -snes_monitor -snes_max_funcs 200000
// -ksp_rtol 1e-6 -ksp_atol 1e-12 -snes_rtol 1e-6