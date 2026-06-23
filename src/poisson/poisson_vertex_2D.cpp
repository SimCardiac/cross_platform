#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <petscsys.h>

#include <cmath>
#include <iostream>
#include <vector>

// ============================================================================
// Manufactured solution: u(x,y) = sin(πx) sin(πy)
// Poisson equation:    -Δu = f  with  f = 2π² sin(πx) sin(πy)
// Domain: [0,1]×[0,1], Dirichlet BC: u = 0 on ∂Ω
// ============================================================================

static inline PetscScalar u_exact(PetscScalar x, PetscScalar y) {
  return std::sin(M_PI * x) * std::sin(M_PI * y);
}

static inline PetscScalar f_rhs(PetscScalar x, PetscScalar y) {
  return 2.0 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
}

// ============================================================================
// Assemble the 5-point Laplacian matrix on a vertex-centered DMDA grid.
// Dirichlet BC: boundary rows are identity (diag=1, RHS=0).
// Interior:  -(u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1} - 4 u_{i,j}) / h²
// ============================================================================
PetscErrorCode AssembleMatrix(DM da, Mat A, PetscInt Nx, PetscInt Ny) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / (Nx - 1);
  const PetscReal hy = 1.0 / (Ny - 1);
  const PetscReal ix2 = 1.0 / (hx * hx);
  const PetscReal iy2 = 1.0 / (hy * hy);
  const PetscReal diag = 2.0 * (ix2 + iy2);

  PetscInt xs, ys, xm, ym;
  PetscCall(DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL));

  for (PetscInt j = ys; j < ys + ym; ++j) {
    for (PetscInt i = xs; i < xs + xm; ++i) {
      MatStencil row;
      row.i = i;
      row.j = j;
      row.k = 0;
      row.c = 0;

      // Boundary nodes: Dirichlet BC u = 0
      if (i == 0 || i == Nx - 1 || j == 0 || j == Ny - 1) {
        PetscScalar v = 1.0;
        PetscCall(MatSetValuesStencil(A, 1, &row, 1, &row, &v, INSERT_VALUES));
      } else {
        // Interior: 5-point stencil
        MatStencil col[5];
        PetscScalar val[5];
        PetscInt ncols = 0;

        // Center
        col[ncols] = row;
        val[ncols] = diag;
        ncols++;

        // Left neighbor
        col[ncols].i = i - 1; col[ncols].j = j; col[ncols].k = 0; col[ncols].c = 0;
        val[ncols] = -ix2;
        ncols++;

        // Right neighbor
        col[ncols].i = i + 1; col[ncols].j = j; col[ncols].k = 0; col[ncols].c = 0;
        val[ncols] = -ix2;
        ncols++;

        // Down neighbor
        col[ncols].i = i; col[ncols].j = j - 1; col[ncols].k = 0; col[ncols].c = 0;
        val[ncols] = -iy2;
        ncols++;

        // Up neighbor
        col[ncols].i = i; col[ncols].j = j + 1; col[ncols].k = 0; col[ncols].c = 0;
        val[ncols] = -iy2;
        ncols++;

        PetscCall(MatSetValuesStencil(A, 1, &row, ncols, col, val, INSERT_VALUES));
      }
    }
  }

  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  PetscFunctionReturn(0);
}

// ============================================================================
// Fill the RHS vector: f[i][j] = 2π² sin(πx_i) sin(πy_j) at interior,
//                         f[i][j] = 0 at boundaries (for Dirichlet BC).
// ============================================================================
PetscErrorCode AssembleRHS(DM da, Vec b, PetscInt Nx, PetscInt Ny) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / (Nx - 1);
  const PetscReal hy = 1.0 / (Ny - 1);

  PetscInt xs, ys, xm, ym;
  PetscCall(DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL));

  PetscScalar **arr;
  PetscCall(DMDAVecGetArray(da, b, &arr));

  for (PetscInt j = ys; j < ys + ym; ++j) {
    for (PetscInt i = xs; i < xs + xm; ++i) {
      if (i == 0 || i == Nx - 1 || j == 0 || j == Ny - 1) {
        arr[j][i] = 0.0; // Dirichlet BC: u=0
      } else {
        const PetscScalar x = i * hx;
        const PetscScalar y = j * hy;
        arr[j][i] = f_rhs(x, y);
      }
    }
  }

  PetscCall(DMDAVecRestoreArray(da, b, &arr));
  PetscFunctionReturn(0);
}

// ============================================================================
// Compute L2 error: ||u - u_exact||_L² = sqrt(hx*hy * Σ(u_ij - u_exact_ij)²)
// ============================================================================
PetscErrorCode ComputeL2Error(DM da, Vec u, PetscInt Nx, PetscInt Ny,
                               PetscReal *l2error) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / (Nx - 1);
  const PetscReal hy = 1.0 / (Ny - 1);

  // Build exact solution vector
  Vec uExact;
  PetscCall(DMCreateGlobalVector(da, &uExact));

  PetscInt xs, ys, xm, ym;
  PetscCall(DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL));

  PetscScalar **arrExact;
  PetscCall(DMDAVecGetArray(da, uExact, &arrExact));
  for (PetscInt j = ys; j < ys + ym; ++j) {
    for (PetscInt i = xs; i < xs + xm; ++i) {
      arrExact[j][i] = u_exact(i * hx, j * hy);
    }
  }
  PetscCall(DMDAVecRestoreArray(da, uExact, &arrExact));

  // Compute difference
  Vec diff;
  PetscCall(VecDuplicate(u, &diff));
  PetscCall(VecCopy(u, diff));
  PetscCall(VecAXPY(diff, -1.0, uExact));

  PetscReal nrm2;
  PetscCall(VecNorm(diff, NORM_2, &nrm2));

  // Scale by sqrt(cell area) for L² norm
  *l2error = nrm2 * std::sqrt(hx * hy);

  PetscCall(VecDestroy(&diff));
  PetscCall(VecDestroy(&uExact));
  PetscFunctionReturn(0);
}

// ============================================================================
// Main: vertex-centered Poisson solver on [0,1]² with manufactured solution
//   Usage: ./poisson_vertex_2D [options]
//   Options:
//     -nx <N>             Number of grid points in x (default: 64)
//     -ny <N>             Number of grid points in y (default: 64)
//     -poisson_check_error          Compute and print L2 error
//     -convergence_test             Output CONVERGENCE line for scripting
//     -ksp_type <type>              KSP solver type (default: gmres)
//     -pc_type <type>               Preconditioner (default: ilu)
// ============================================================================
int main(int argc, char **argv) {
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL,
                            "Vertex-centered Poisson solver on DMDA"));

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // ---- Parse options ----
  PetscInt Nx = 64, Ny = 64;
  PetscBool compute_error = PETSC_FALSE;
  PetscBool convergence_test = PETSC_FALSE;

  PetscCall(PetscOptionsGetInt(NULL, NULL, "-nx", &Nx, NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-ny", &Ny, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-poisson_check_error", &compute_error, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-convergence_test", &convergence_test, NULL));

  if (Nx < 3 || Ny < 3) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Error: Nx, Ny must be >= 3\n"));
    PetscCall(PetscFinalize());
    return 1;
  }

  // ---- Create DMDA (vertex-centered structured grid) ----
  DM da;
  PetscCall(DMDACreate2d(PETSC_COMM_WORLD,
                         DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                         DMDA_STENCIL_STAR,
                         Nx, Ny,
                         PETSC_DECIDE, PETSC_DECIDE,
                         1,       // dof per node
                         1,       // stencil width
                         NULL, NULL,
                         &da));
  PetscCall(DMSetFromOptions(da));
  PetscCall(DMSetUp(da));

  // Set uniform coordinates [0,1]×[0,1]
  PetscCall(DMDASetUniformCoordinates(da, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0));

  // ---- Create matrix, solution, and RHS ----
  Mat A;
  Vec u, b;
  PetscCall(DMCreateMatrix(da, &A));
  PetscCall(DMCreateGlobalVector(da, &u));
  PetscCall(DMCreateGlobalVector(da, &b));

  PetscCall(VecSet(u, 0.0));

  // ---- Assemble system ----
  PetscCall(AssembleMatrix(da, A, Nx, Ny));
  PetscCall(AssembleRHS(da, b, Nx, Ny));

  // ---- Set up KSP solver ----
  KSP ksp;
  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
  PetscCall(KSPSetOperators(ksp, A, A));
  PetscCall(KSPSetFromOptions(ksp));
  PetscCall(KSPSetTolerances(ksp, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, 10000));

  // ---- Solve ----
  if (rank == 0) {
    std::cout << "Vertex-centered Poisson: Nx=" << Nx << ", Ny=" << Ny
              << ", h=" << 1.0 / (Nx - 1) << std::endl;
  }
  PetscCall(KSPSolve(ksp, b, u));

  PetscInt its;
  PetscReal rnorm;
  PetscCall(KSPGetIterationNumber(ksp, &its));
  PetscCall(KSPGetResidualNorm(ksp, &rnorm));

  KSPConvergedReason reason;
  PetscCall(KSPGetConvergedReason(ksp, &reason));

  if (rank == 0) {
    std::cout << "KSP iterations: " << its
              << ", residual: " << rnorm
              << ", reason: " << reason << std::endl;
  }

  // ---- Compute L2 error ----
  PetscReal l2Error = 0.0;
  if (compute_error) {
    PetscCall(ComputeL2Error(da, u, Nx, Ny, &l2Error));
    if (rank == 0) {
      std::cout << "||u - u_exact||_L2 = " << l2Error << std::endl;
    }
  }

  // ---- Convergence test mode ----
  if (convergence_test && rank == 0) {
    const PetscReal h = 1.0 / (Nx - 1);
    std::cout << "CONVERGENCE: " << Nx << " " << h << " " << l2Error
              << " " << its << std::endl;
  }

  // ---- Cleanup ----
  PetscCall(KSPDestroy(&ksp));
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&b));
  PetscCall(MatDestroy(&A));
  PetscCall(DMDestroy(&da));
  PetscCall(PetscFinalize());
  return 0;
}
