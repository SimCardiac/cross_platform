#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <petscsys.h>

#include <cmath>
#include <iostream>

// ============================================================================
// Vertex-centered heat equation on DMDA
//   ∂u/∂t = α Δu + f   on [0,1]²,  u=0 on boundary
//   Implicit Euler: (I - αΔt Δ_h) u^{n+1} = u^n + Δt f^{n+1}
//
// Manufactured solution: u(x,y,t) = e^{-2π²αt} sin(πx) sin(πy)
//   → f = 0,  u(x,y,0) = sin(πx) sin(πy)
// ============================================================================

static const PetscScalar alpha = 0.1;  // thermal diffusivity

static inline PetscScalar u_exact(PetscScalar x, PetscScalar y, PetscScalar t) {
  return std::exp(-2.0 * M_PI * M_PI * alpha * t) *
         std::sin(M_PI * x) * std::sin(M_PI * y);
}

static inline PetscScalar u_init(PetscScalar x, PetscScalar y) {
  return std::sin(M_PI * x) * std::sin(M_PI * y);
}

// ============================================================================
// Assemble Helmholtz matrix: A = I - αΔt Δ_h  (including boundary rows)
// ============================================================================
PetscErrorCode AssembleHelmholtz(Mat A, PetscInt Nx, PetscInt Ny,
                                  PetscReal dt) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / (Nx - 1);
  const PetscReal hy = 1.0 / (Ny - 1);
  const PetscReal c = alpha * dt;
  const PetscReal ix2 = c / (hx * hx);
  const PetscReal iy2 = c / (hy * hy);
  const PetscReal diag = 1.0 + 2.0 * (ix2 + iy2);

  PetscInt rstart, rend;
  PetscCall(MatGetOwnershipRange(A, &rstart, &rend));

  for (PetscInt j = 0; j < Ny; ++j) {
    for (PetscInt i = 0; i < Nx; ++i) {
      PetscInt row = j * Nx + i;
      if (row < rstart || row >= rend) continue;

      if (i == 0 || i == Nx - 1 || j == 0 || j == Ny - 1) {
        // Dirichlet BC: identity row
        PetscCall(MatSetValue(A, row, row, 1.0, INSERT_VALUES));
      } else {
        // Interior: 5-point Helmholtz stencil
        PetscCall(MatSetValue(A, row, row,                          diag, INSERT_VALUES));
        PetscCall(MatSetValue(A, row, j * Nx + (i - 1),             -ix2, INSERT_VALUES));
        PetscCall(MatSetValue(A, row, j * Nx + (i + 1),             -ix2, INSERT_VALUES));
        PetscCall(MatSetValue(A, row, (j - 1) * Nx + i,             -iy2, INSERT_VALUES));
        PetscCall(MatSetValue(A, row, (j + 1) * Nx + i,             -iy2, INSERT_VALUES));
      }
    }
  }

  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  PetscFunctionReturn(0);
}

// ============================================================================
// Fill initial condition
// ============================================================================
PetscErrorCode SetInitialCondition(DM da, Vec u, PetscInt Nx, PetscInt Ny) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / (Nx - 1);
  const PetscReal hy = 1.0 / (Ny - 1);

  PetscInt xs, ys, xm, ym;
  PetscCall(DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL));

  PetscScalar **arr;
  PetscCall(DMDAVecGetArray(da, u, &arr));
  for (PetscInt j = ys; j < ys + ym; ++j) {
    for (PetscInt i = xs; i < xs + xm; ++i) {
      arr[j][i] = u_init(i * hx, j * hy);
    }
  }
  PetscCall(DMDAVecRestoreArray(da, u, &arr));
  PetscFunctionReturn(0);
}

// ============================================================================
// Compute L2 error at final time
// ============================================================================
PetscErrorCode ComputeL2Error(DM da, Vec u, PetscInt Nx, PetscInt Ny,
                               PetscScalar t, PetscReal *l2err) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / (Nx - 1);
  const PetscReal hy = 1.0 / (Ny - 1);

  Vec uExact;
  PetscCall(DMCreateGlobalVector(da, &uExact));

  PetscInt xs, ys, xm, ym;
  PetscCall(DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL));

  PetscScalar **arrE;
  PetscCall(DMDAVecGetArray(da, uExact, &arrE));
  for (PetscInt j = ys; j < ys + ym; ++j) {
    for (PetscInt i = xs; i < xs + xm; ++i) {
      arrE[j][i] = u_exact(i * hx, j * hy, t);
    }
  }
  PetscCall(DMDAVecRestoreArray(da, uExact, &arrE));

  Vec diff;
  PetscCall(VecDuplicate(u, &diff));
  PetscCall(VecCopy(u, diff));
  PetscCall(VecAXPY(diff, -1.0, uExact));

  PetscReal nrm2;
  PetscCall(VecNorm(diff, NORM_2, &nrm2));
  *l2err = nrm2 * std::sqrt(hx * hy);

  PetscCall(VecDestroy(&diff));
  PetscCall(VecDestroy(&uExact));
  PetscFunctionReturn(0);
}

// ============================================================================
// Main
// ============================================================================
int main(int argc, char **argv) {
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL,
                            "Vertex-centered heat equation on DMDA"));

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscInt Nx = 64, Ny = 64;
  PetscReal Tfinal = 0.05;
  PetscBool compute_error = PETSC_FALSE;
  PetscBool convergence_test = PETSC_FALSE;

  PetscCall(PetscOptionsGetInt(NULL, NULL, "-nx", &Nx, NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-ny", &Ny, NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-T", &Tfinal, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-heat_check_error", &compute_error, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-convergence_test", &convergence_test, NULL));

  if (Nx < 3 || Ny < 3) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Error: Nx, Ny >= 3\n"));
    PetscCall(PetscFinalize());
    return 1;
  }

  // Time step: for convergence test, use Δt ∝ h² (coupled refinement)
  const PetscReal h = 1.0 / (Nx - 1);
  const PetscReal dt = h * h;  // Δt ∝ h² for second-order spatial-temporal coupling
  const PetscInt Nsteps = static_cast<PetscInt>(std::ceil(Tfinal / dt));
  const PetscReal dtActual = Tfinal / Nsteps;

  // ---- Create DMDA ----
  DM da;
  PetscCall(DMDACreate2d(PETSC_COMM_WORLD,
                         DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                         DMDA_STENCIL_STAR,
                         Nx, Ny, PETSC_DECIDE, PETSC_DECIDE,
                         1, 1, NULL, NULL, &da));
  PetscCall(DMSetFromOptions(da));
  PetscCall(DMSetUp(da));
  PetscCall(DMDASetUniformCoordinates(da, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0));

  // ---- Create matrix, vectors ----
  Mat A;
  PetscCall(DMCreateMatrix(da, &A));
  PetscCall(AssembleHelmholtz(A, Nx, Ny, dtActual));

  Vec u, uOld, b;
  PetscCall(DMCreateGlobalVector(da, &u));
  PetscCall(DMCreateGlobalVector(da, &uOld));
  PetscCall(DMCreateGlobalVector(da, &b));

  PetscCall(SetInitialCondition(da, u, Nx, Ny));

  // ---- KSP solver ----
  KSP ksp;
  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
  PetscCall(KSPSetOperators(ksp, A, A));
  PetscCall(KSPSetType(ksp, KSPGMRES));
  PetscCall(KSPSetTolerances(ksp, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT, 2000));
  PetscCall(KSPSetFromOptions(ksp));

  if (rank == 0) {
    std::cout << "Vertex-centered heat: N=" << Nx << "x" << Ny
              << ", h=" << h << ", dt=" << dtActual
              << ", steps=" << Nsteps << ", T=" << Tfinal << std::endl;
  }

  // ---- Time stepping ----
  for (PetscInt step = 1; step <= Nsteps; ++step) {
    PetscCall(VecCopy(u, uOld));          // uOld = u^n
    PetscCall(VecCopy(uOld, b));          // RHS = u^n (since f=0)
    PetscCall(KSPSolve(ksp, b, u));       // solve Helmholtz

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
    PetscCall(ComputeL2Error(da, u, Nx, Ny, Tfinal, &l2Error));
    if (rank == 0) {
      std::cout << "||u - u_exact||_L2 = " << l2Error << std::endl;
    }
  }

  if (convergence_test && rank == 0) {
    std::cout << "CONVERGENCE: " << Nx << " " << h << " " << l2Error
              << " " << Nsteps << std::endl;
  }

  // ---- Cleanup ----
  PetscCall(KSPDestroy(&ksp));
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&uOld));
  PetscCall(VecDestroy(&b));
  PetscCall(MatDestroy(&A));
  PetscCall(DMDestroy(&da));
  PetscCall(PetscFinalize());
  return 0;
}
