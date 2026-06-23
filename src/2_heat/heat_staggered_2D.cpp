#include <petscdm.h>
#include <petscdmstag.h>
#include <petscksp.h>
#include <petscsys.h>

#include <cmath>
#include <iostream>

#include "analytical/unsteady.h"
using namespace UNSTEADY::SINPI_SCALAR_2D;

// ============================================================================
// Staggered (mixed first-order) heat equation on DMStag (Item 6)
//
// First-order system:
//   gx - ∂u/∂x = 0       (left faces)
//   gy - ∂u/∂y = 0       (down faces)
//   ∂u/∂t - α(∂gx/∂x + ∂gy/∂y) = 0   (elements)
//
// Implicit Euler: [ I     0     -Gx      ] [gx^{n+1}]   [0    ]
//                 [ 0     I     -Gy      ] [gy^{n+1}] = [0    ]
//                 [-αΔtDx -αΔtDy  I     ] [u^{n+1}]    [u^n  ]
//
// Manufactured solution (SINPI_SCALAR_2D):
//   u = exp(-2π²αt) sin(πx) sin(πy)
// ============================================================================

static const PetscScalar alpha = 0.1;

// DOF indexing (same as staggered Poisson)
static inline PetscInt idx_u(PetscInt ex, PetscInt ey, PetscInt Nx, PetscInt Ny) {
  (void)Ny; return ey * Nx + ex;
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
// Assemble the monolithic system for one time step
//   [ I     0     -Gx      ] [gx]
//   [ 0     I     -Gy      ] [gy]
//   [-c*Dx  -c*Dy    I     ] [u ]    where c = alpha*dt
// ============================================================================
PetscErrorCode AssembleMatrix(Mat A, PetscInt Nx, PetscInt Ny, PetscReal dt) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  const PetscReal c = alpha * dt;
  const PetscReal ihx = 1.0 / hx;
  const PetscReal ihx2 = 2.0 / hx;   // boundary: distance = h/2
  const PetscReal ihy = 1.0 / hy;
  const PetscReal ihy2 = 2.0 / hy;
  const PetscReal cihx = c * ihx;
  const PetscReal cihx2 = c * ihx2;
  const PetscReal cihy = c * ihy;
  const PetscReal cihy2 = c * ihy2;

  PetscInt rstart, rend;
  PetscCall(MatGetOwnershipRange(A, &rstart, &rend));

  // ---- gx equations: gx - Gx u = 0 ----
  for (PetscInt fy = 0; fy < Ny; ++fy) {
    for (PetscInt fx = 0; fx <= Nx; ++fx) {
      PetscInt row = idx_gx(fx, fy, Nx, Ny);
      if (row < rstart || row >= rend) continue;
      PetscCall(MatSetValue(A, row, row, 1.0, INSERT_VALUES));

      if (fx == 0) {
        PetscInt col = idx_u(0, fy, Nx, Ny);
        PetscCall(MatSetValue(A, row, col, -ihx2, INSERT_VALUES));
      } else if (fx == Nx) {
        PetscInt col = idx_u(Nx - 1, fy, Nx, Ny);
        PetscCall(MatSetValue(A, row, col, ihx2, INSERT_VALUES));
      } else {
        PetscCall(MatSetValue(A, row, idx_u(fx, fy, Nx, Ny), -ihx, INSERT_VALUES));
        PetscCall(MatSetValue(A, row, idx_u(fx - 1, fy, Nx, Ny), ihx, INSERT_VALUES));
      }
    }
  }

  // ---- gy equations: gy - Gy u = 0 ----
  for (PetscInt fx = 0; fx < Nx; ++fx) {
    for (PetscInt fy = 0; fy <= Ny; ++fy) {
      PetscInt row = idx_gy(fx, fy, Nx, Ny);
      if (row < rstart || row >= rend) continue;
      PetscCall(MatSetValue(A, row, row, 1.0, INSERT_VALUES));

      if (fy == 0) {
        PetscCall(MatSetValue(A, row, idx_u(fx, 0, Nx, Ny), -ihy2, INSERT_VALUES));
      } else if (fy == Ny) {
        PetscCall(MatSetValue(A, row, idx_u(fx, Ny - 1, Nx, Ny), ihy2, INSERT_VALUES));
      } else {
        PetscCall(MatSetValue(A, row, idx_u(fx, fy, Nx, Ny), -ihy, INSERT_VALUES));
        PetscCall(MatSetValue(A, row, idx_u(fx, fy - 1, Nx, Ny), ihy, INSERT_VALUES));
      }
    }
  }

  // ---- Heat equation on elements: u - c*(Dx gx + Dy gy) = u_old ----
  for (PetscInt ey = 0; ey < Ny; ++ey) {
    for (PetscInt ex = 0; ex < Nx; ++ex) {
      PetscInt row = idx_u(ex, ey, Nx, Ny);
      if (row < rstart || row >= rend) continue;

      // Diagonal: identity for u
      PetscCall(MatSetValue(A, row, row, 1.0, INSERT_VALUES));

      // -c*Dx: -(gx[ex+1,ey] - gx[ex,ey]) * c/h
      PetscCall(MatSetValue(A, row, idx_gx(ex + 1, ey, Nx, Ny), -cihx, INSERT_VALUES));
      PetscCall(MatSetValue(A, row, idx_gx(ex, ey, Nx, Ny), cihx, INSERT_VALUES));

      // -c*Dy: -(gy[ex,ey+1] - gy[ex,ey]) * c/h
      PetscCall(MatSetValue(A, row, idx_gy(ex, ey + 1, Nx, Ny), -cihy, INSERT_VALUES));
      PetscCall(MatSetValue(A, row, idx_gy(ex, ey, Nx, Ny), cihy, INSERT_VALUES));
    }
  }

  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  PetscFunctionReturn(0);
}

// ============================================================================
// Set initial condition: u at elements, gx/gy derived from exact gradients
// ============================================================================
PetscErrorCode SetInitialCondition(Vec x, PetscInt Nx, PetscInt Ny) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;

  PetscInt rstart, rend;
  PetscCall(VecGetOwnershipRange(x, &rstart, &rend));

  // u at element centers
  for (PetscInt ey = 0; ey < Ny; ++ey) {
    for (PetscInt ex = 0; ex < Nx; ++ex) {
      PetscInt row = idx_u(ex, ey, Nx, Ny);
      if (row < rstart || row >= rend) continue;
      const PetscScalar xc = (ex + 0.5) * hx;
      const PetscScalar yc = (ey + 0.5) * hy;
      PetscCall(VecSetValue(x, row, u_steady(xc, yc), INSERT_VALUES));
    }
  }

  // gx at left faces (initial: exact gradient)
  for (PetscInt fy = 0; fy < Ny; ++fy) {
    for (PetscInt fx = 0; fx <= Nx; ++fx) {
      PetscInt row = idx_gx(fx, fy, Nx, Ny);
      if (row < rstart || row >= rend) continue;
      const PetscScalar xf = fx * hx;
      const PetscScalar yc = (fy + 0.5) * hy;
      PetscCall(VecSetValue(x, row, gx_exact(xf, yc), INSERT_VALUES));
    }
  }

  // gy at down faces
  for (PetscInt fx = 0; fx < Nx; ++fx) {
    for (PetscInt fy = 0; fy <= Ny; ++fy) {
      PetscInt row = idx_gy(fx, fy, Nx, Ny);
      if (row < rstart || row >= rend) continue;
      const PetscScalar xc = (fx + 0.5) * hx;
      const PetscScalar yf = fy * hy;
      PetscCall(VecSetValue(x, row, gy_exact(xc, yf), INSERT_VALUES));
    }
  }

  PetscCall(VecAssemblyBegin(x));
  PetscCall(VecAssemblyEnd(x));
  PetscFunctionReturn(0);
}

// ============================================================================
// Build RHS for step: copy u^n into the RHS, zeros for gx/gy equations
// ============================================================================
PetscErrorCode BuildRHS(Vec b, const Vec xOld, PetscInt Nx, PetscInt Ny) {
  PetscFunctionBeginUser;
  PetscCall(VecSet(b, 0.0));

  PetscInt rstart, rend;
  PetscCall(VecGetOwnershipRange(b, &rstart, &rend));

  // Copy u_old values into RHS (for u rows only)
  for (PetscInt ey = 0; ey < Ny; ++ey) {
    for (PetscInt ex = 0; ex < Nx; ++ex) {
      PetscInt row = idx_u(ex, ey, Nx, Ny);
      if (row < rstart || row >= rend) continue;
      PetscScalar val;
      PetscCall(VecGetValues(xOld, 1, &row, &val));
      PetscCall(VecSetValue(b, row, val, INSERT_VALUES));
    }
  }

  PetscCall(VecAssemblyBegin(b));
  PetscCall(VecAssemblyEnd(b));
  PetscFunctionReturn(0);
}

// ============================================================================
// Compute L2 error for u (element DOFs)
// ============================================================================
PetscErrorCode ComputeL2Error(Vec x, PetscInt Nx, PetscInt Ny,
                               PetscReal t, PetscReal *l2err) {
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
      const PetscScalar diff = val - u_exact(xc, yc, t, alpha);
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
                            "Staggered (mixed) heat equation on DMStag"));

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscInt Nx = 64, Ny = 64;
  PetscReal T_final = 0.05;
  PetscBool compute_error = PETSC_FALSE;
  PetscBool convergence_test = PETSC_FALSE;

  PetscCall(PetscOptionsGetInt(NULL, NULL, "-nx", &Nx, NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-ny", &Ny, NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-T", &T_final, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-heat_check_error", &compute_error, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-convergence_test", &convergence_test, NULL));

  // Time step: dt ∝ h² for convergence test
  const PetscReal h = 1.0 / Nx;
  const PetscReal dt = h * h;
  const PetscInt Nsteps = static_cast<PetscInt>(std::ceil(T_final / dt));
  const PetscReal dtActual = T_final / Nsteps;

  PetscInt Ntotal = totalDOFs(Nx, Ny);

  // ---- Create monolithic matrix ----
  Mat A;
  PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
  PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, Ntotal, Ntotal));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatSetUp(A));
  PetscCall(MatSeqAIJSetPreallocation(A, 5, NULL));
  PetscCall(MatMPIAIJSetPreallocation(A, 5, NULL, 5, NULL));

  PetscCall(AssembleMatrix(A, Nx, Ny, dtActual));

  Vec x, xOld, b;
  PetscCall(MatCreateVecs(A, &x, &b));
  PetscCall(VecDuplicate(x, &xOld));

  // Initial condition
  PetscCall(SetInitialCondition(x, Nx, Ny));

  // ---- KSP solver ----
  KSP ksp;
  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
  PetscCall(KSPSetOperators(ksp, A, A));
  PetscCall(KSPSetType(ksp, KSPGMRES));
  PetscCall(KSPSetTolerances(ksp, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT, 2000));
  PetscCall(KSPSetFromOptions(ksp));

  if (rank == 0) {
    std::cout << "Staggered heat: N=" << Nx << "x" << Ny
              << ", h=" << h << ", dt=" << dtActual
              << ", steps=" << Nsteps << ", DOFs=" << Ntotal << std::endl;
  }

  // ---- Time stepping ----
  PetscReal t = 0.0;
  for (PetscInt step = 1; step <= Nsteps; ++step) {
    t += dtActual;
    PetscCall(VecCopy(x, xOld));
    PetscCall(BuildRHS(b, xOld, Nx, Ny));
    PetscCall(KSPSolve(ksp, b, x));

    if (rank == 0 && (step % 100 == 0 || step == Nsteps)) {
      PetscInt its;
      PetscCall(KSPGetIterationNumber(ksp, &its));
      std::cout << "  step " << step << "/" << Nsteps
                << ", KSP its=" << its << std::endl;
    }
  }

  // ---- Compute error ----
  PetscReal l2Error = 0.0;
  if (compute_error) {
    PetscCall(ComputeL2Error(x, Nx, Ny, t, &l2Error));
    if (rank == 0) {
      std::cout << "||u(T) - u_exact(T)||_L2 = " << l2Error << std::endl;
    }
  }

  if (convergence_test && rank == 0) {
    std::cout << "CONVERGENCE: " << Nx << " " << h << " " << l2Error
              << " " << Nsteps << std::endl;
  }

  // ---- Cleanup ----
  PetscCall(KSPDestroy(&ksp));
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&xOld));
  PetscCall(VecDestroy(&b));
  PetscCall(MatDestroy(&A));
  PetscCall(PetscFinalize());
  return 0;
}
