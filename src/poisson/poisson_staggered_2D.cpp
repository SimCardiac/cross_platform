#include <petscdm.h>
#include <petscdmstag.h>
#include <petscksp.h>
#include <petscsys.h>

#include <cmath>
#include <iostream>

#include "analytical/unsteady.h"
using namespace UNSTEADY::SINPI_SCALAR_2D;

// ============================================================================
// Staggered (mixed first-order) Poisson solver on DMStag
//
// Unknowns:
//   dof0: u  at element centers       (Nx × Ny)
//   dof1: gx at left faces            ((Nx+1) × Ny)
//   dof2: gy at down faces            (Nx × (Ny+1))
//
// First-order system equivalent to -Δu = f:
//   gx - ∂u/∂x = 0   (on left faces)
//   gy - ∂u/∂y = 0   (on down faces)
//   -(∂gx/∂x + ∂gy/∂y) = f   (on elements)
//
// Manufactured solution: u = sin(πx)sin(πy), f = 2π²sin(πx)sin(πy)
// (from SINPI_SCALAR_2D in analytical/unsteady.h)
// ============================================================================

// ============================================================================
// Assembly helpers
// ============================================================================

// Get row/col indices for each DOF type in the monolithic system.
// Ordering: element DOFs first, then left-face DOFs, then down-face DOFs.
static inline PetscInt idx_u(PetscInt ex, PetscInt ey, PetscInt Nx, PetscInt Ny) {
  (void)Ny;
  return ey * Nx + ex;
}
static inline PetscInt idx_gx(PetscInt fx, PetscInt fy, PetscInt Nx, PetscInt Ny) {
  return Nx * Ny + fy * (Nx + 1) + fx;
}
static inline PetscInt idx_gy(PetscInt fx, PetscInt fy, PetscInt Nx, PetscInt Ny) {
  return Nx * Ny + (Nx + 1) * Ny + fx * (Ny + 1) + fy;
}
static inline PetscInt totalDOFs(PetscInt Nx, PetscInt Ny) {
  return Nx * Ny + (Nx + 1) * Ny + Nx * (Ny + 1);
}

// ============================================================================
// Assemble the monolithic mixed system matrix
//   [  I     0    -Gx ] [gx]   [0]
//   [  0     I    -Gy ] [gy] = [0]
//   [ -Dx  -Dy    0   ] [u ]   [f]
// ============================================================================
PetscErrorCode AssembleMatrix(Mat A, PetscInt Nx, PetscInt Ny, PetscInt Nlocal) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  const PetscReal ihx = 1.0 / hx;   // for interior gradient (distance = h)
  const PetscReal ihx2 = 2.0 / hx;  // for boundary gradient (distance = h/2)
  const PetscReal ihy = 1.0 / hy;
  const PetscReal ihy2 = 2.0 / hy;

  PetscInt rstart, rend;
  PetscCall(MatGetOwnershipRange(A, &rstart, &rend));

  // ---- gx equations: gx - Gx u = 0 (on left faces) ----
  for (PetscInt fy = 0; fy < Ny; ++fy) {
    for (PetscInt fx = 0; fx <= Nx; ++fx) {
      PetscInt row = idx_gx(fx, fy, Nx, Ny);
      if (row < rstart || row >= rend) continue;

      // Diagonal: identity for gx
      PetscCall(MatSetValue(A, row, row, 1.0, INSERT_VALUES));

      if (fx == 0) {
        // x=0 boundary: distance to first cell center = h/2
        // gx = (u[0,fy] - 0) / (h/2) = 2*u[0,fy]/h
        PetscInt col = idx_u(0, fy, Nx, Ny);
        PetscCall(MatSetValue(A, row, col, -ihx2, INSERT_VALUES));
      } else if (fx == Nx) {
        // x=1 boundary: distance from last cell center = h/2
        // gx = (0 - u[Nx-1,fy]) / (h/2) = -2*u[Nx-1,fy]/h
        PetscInt col = idx_u(Nx - 1, fy, Nx, Ny);
        PetscCall(MatSetValue(A, row, col, ihx2, INSERT_VALUES));
      } else {
        // Interior: gx = (u[fx,fy] - u[fx-1,fy])/h
        PetscInt colR = idx_u(fx, fy, Nx, Ny);
        PetscInt colL = idx_u(fx - 1, fy, Nx, Ny);
        PetscCall(MatSetValue(A, row, colR, -ihx, INSERT_VALUES));
        PetscCall(MatSetValue(A, row, colL, ihx, INSERT_VALUES));
      }
    }
  }

  // ---- gy equations: gy - Gy u = 0 (on down faces) ----
  for (PetscInt fx = 0; fx < Nx; ++fx) {
    for (PetscInt fy = 0; fy <= Ny; ++fy) {
      PetscInt row = idx_gy(fx, fy, Nx, Ny);
      if (row < rstart || row >= rend) continue;

      // Diagonal: identity for gy
      PetscCall(MatSetValue(A, row, row, 1.0, INSERT_VALUES));

      if (fy == 0) {
        // y=0 boundary
        PetscInt col = idx_u(fx, 0, Nx, Ny);
        PetscCall(MatSetValue(A, row, col, -ihy2, INSERT_VALUES));
      } else if (fy == Ny) {
        // y=1 boundary
        PetscInt col = idx_u(fx, Ny - 1, Nx, Ny);
        PetscCall(MatSetValue(A, row, col, ihy2, INSERT_VALUES));
      } else {
        // Interior
        PetscInt colU = idx_u(fx, fy, Nx, Ny);
        PetscInt colD = idx_u(fx, fy - 1, Nx, Ny);
        PetscCall(MatSetValue(A, row, colU, -ihy, INSERT_VALUES));
        PetscCall(MatSetValue(A, row, colD, ihy, INSERT_VALUES));
      }
    }
  }

  // ---- Divergence equation: -(Dx gx + Dy gy) = f (on elements) ----
  for (PetscInt ey = 0; ey < Ny; ++ey) {
    for (PetscInt ex = 0; ex < Nx; ++ex) {
      PetscInt row = idx_u(ex, ey, Nx, Ny);
      if (row < rstart || row >= rend) continue;

      // -(gx[ex+1,ey] - gx[ex,ey])/h - (gy[ex,ey+1] - gy[ex,ey])/h = f
      PetscInt col_gxR = idx_gx(ex + 1, ey, Nx, Ny);
      PetscInt col_gxL = idx_gx(ex, ey, Nx, Ny);
      PetscInt col_gyU = idx_gy(ex, ey + 1, Nx, Ny);
      PetscInt col_gyD = idx_gy(ex, ey, Nx, Ny);

      PetscCall(MatSetValue(A, row, col_gxR, -ihx, INSERT_VALUES));
      PetscCall(MatSetValue(A, row, col_gxL, ihx, INSERT_VALUES));
      PetscCall(MatSetValue(A, row, col_gyU, -ihy, INSERT_VALUES));
      PetscCall(MatSetValue(A, row, col_gyD, ihy, INSERT_VALUES));
    }
  }

  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  PetscFunctionReturn(0);
}

// ============================================================================
// Fill the RHS vector [0; 0; f]
// ============================================================================
PetscErrorCode AssembleRHS(Vec b, PetscInt Nx, PetscInt Ny) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;

  PetscCall(VecSet(b, 0.0));

  PetscInt rstart, rend;
  PetscCall(VecGetOwnershipRange(b, &rstart, &rend));

  for (PetscInt ey = 0; ey < Ny; ++ey) {
    for (PetscInt ex = 0; ex < Nx; ++ex) {
      PetscInt row = idx_u(ex, ey, Nx, Ny);
      if (row < rstart || row >= rend) continue;
      const PetscScalar x = (ex + 0.5) * hx;
      const PetscScalar y = (ey + 0.5) * hy;
      PetscCall(VecSetValue(b, row, f_poisson(x, y), INSERT_VALUES));
    }
  }

  PetscCall(VecAssemblyBegin(b));
  PetscCall(VecAssemblyEnd(b));
  PetscFunctionReturn(0);
}

// ============================================================================
// Compute L2 error for u (element DOFs only)
// ============================================================================
PetscErrorCode ComputeL2Error(Vec x, PetscInt Nx, PetscInt Ny, PetscReal *l2err) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;

  PetscReal localSum = 0.0;
  PetscInt rstart, rend;
  PetscCall(VecGetOwnershipRange(x, &rstart, &rend));

  for (PetscInt ey = 0; ey < Ny; ++ey) {
    for (PetscInt ex = 0; ex < Nx; ++ex) {
      PetscInt row = idx_u(ex, ey, Nx, Ny);
      if (row < rstart || row >= rend) continue;

      PetscScalar val;
      PetscCall(VecGetValues(x, 1, &row, &val));

      const PetscScalar xc = (ex + 0.5) * hx;
      const PetscScalar yc = (ey + 0.5) * hy;
      const PetscScalar diff = val - u_steady(xc, yc);
      localSum += PetscRealPart(diff * diff);
    }
  }

  PetscReal globalSum;
  MPI_Allreduce(&localSum, &globalSum, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
  *l2err = std::sqrt(globalSum * hx * hy);
  PetscFunctionReturn(0);
}

// ============================================================================
// Main
// ============================================================================
int main(int argc, char **argv) {
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL,
                            "Staggered (mixed) Poisson solver on DMStag"));

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscInt Nx = 64, Ny = 64;
  PetscBool compute_error = PETSC_FALSE;
  PetscBool convergence_test = PETSC_FALSE;

  PetscCall(PetscOptionsGetInt(NULL, NULL, "-nx", &Nx, NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-ny", &Ny, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-poisson_check_error", &compute_error, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-convergence_test", &convergence_test, NULL));

  if (Nx < 2 || Ny < 2) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Error: Nx, Ny must be >= 2\n"));
    PetscCall(PetscFinalize());
    return 1;
  }

  // ---- Global size (for monolithic matrix) ----
  PetscInt Ntotal = totalDOFs(Nx, Ny);

  // ---- Create monolithic matrix and vectors ----
  Mat A;
  PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
  PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, Ntotal, Ntotal));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatSetUp(A));

  // Pre-allocate: gx/gy rows ≤3, u (divergence) rows ≤4
  PetscCall(MatSeqAIJSetPreallocation(A, 5, NULL));
  PetscCall(MatMPIAIJSetPreallocation(A, 5, NULL, 5, NULL));

  PetscCall(AssembleMatrix(A, Nx, Ny, Ntotal));

  Vec b, x;
  PetscCall(MatCreateVecs(A, &x, &b));
  PetscCall(VecSet(x, 0.0));
  PetscCall(AssembleRHS(b, Nx, Ny));

  // ---- Solve ----
  KSP ksp;
  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
  PetscCall(KSPSetOperators(ksp, A, A));
  PetscCall(KSPSetType(ksp, KSPGMRES));
  PetscCall(KSPSetTolerances(ksp, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, 10000));
  // Default: no preconditioner (saddle-point — u-block has zero diagonal)
  {
    PC pc;
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, PCNONE));
  }
  // Allow command-line options to override everything above
  PetscCall(KSPSetFromOptions(ksp));

  if (rank == 0) {
    std::cout << "Staggered (mixed) Poisson: Nx=" << Nx << ", Ny=" << Ny
              << ", h=" << 1.0 / Nx << ", total DOFs=" << Ntotal << std::endl;
  }

  PetscCall(KSPSolve(ksp, b, x));

  PetscInt its;
  PetscReal rnorm;
  PetscCall(KSPGetIterationNumber(ksp, &its));
  PetscCall(KSPGetResidualNorm(ksp, &rnorm));

  if (rank == 0) {
    KSPConvergedReason reason;
    PetscCall(KSPGetConvergedReason(ksp, &reason));
    std::cout << "KSP iterations: " << its
              << ", residual: " << rnorm
              << ", reason: " << reason << std::endl;
  }

  // ---- Compute L2 error for u ----
  PetscReal l2Error = 0.0;
  if (compute_error) {
    PetscCall(ComputeL2Error(x, Nx, Ny, &l2Error));
    if (rank == 0) {
      std::cout << "||u - u_exact||_L2 = " << l2Error << std::endl;
    }
  }

  if (convergence_test && rank == 0) {
    const PetscReal h = 1.0 / Nx;
    std::cout << "CONVERGENCE: " << Nx << " " << h << " " << l2Error
              << " " << its << std::endl;
  }

  // ---- Cleanup ----
  PetscCall(KSPDestroy(&ksp));
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&b));
  PetscCall(MatDestroy(&A));
  PetscCall(PetscFinalize());
  return 0;
}
