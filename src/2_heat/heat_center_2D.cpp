#include <petscdm.h>
#include <petscdmstag.h>
#include <petscksp.h>
#include <petscsys.h>

#include <cmath>
#include <iostream>

#include "analytical/unsteady.h"

// ============================================================================
// Inline DMStag boundary helpers (ghost-cell reflection for Dirichlet BC)
// ============================================================================
static inline PetscScalar DMStag_GetLeft(PetscScalar ***arr, PetscInt ex, PetscInt ey,
                                          PetscInt slot, PetscInt Nx) {
  if (ex == 0) return -arr[ey][ex][slot];
  return arr[ey][ex - 1][slot];
}
static inline PetscScalar DMStag_GetRight(PetscScalar ***arr, PetscInt ex, PetscInt ey,
                                           PetscInt slot, PetscInt Nx) {
  if (ex == Nx - 1) return -arr[ey][ex][slot];
  return arr[ey][ex + 1][slot];
}
static inline PetscScalar DMStag_GetDown(PetscScalar ***arr, PetscInt ex, PetscInt ey,
                                          PetscInt slot, PetscInt Ny) {
  if (ey == 0) return -arr[ey][ex][slot];
  return arr[ey - 1][ex][slot];
}
static inline PetscScalar DMStag_GetUp(PetscScalar ***arr, PetscInt ex, PetscInt ey,
                                        PetscInt slot, PetscInt Ny) {
  if (ey == Ny - 1) return -arr[ey][ex][slot];
  return arr[ey + 1][ex][slot];
}

// ============================================================================
// Cell-centered heat equation on DMStag (Item 5)
//   ∂u/∂t = α Δu   on [0,1]²,  u=0 on boundary
//   Implicit Euler: (I - αΔt Δ_h) u^{n+1} = u^n
//
// Manufactured solution (SINPI_SCALAR_2D):
//   u(x,y,t) = exp(-2π²αt) sin(πx) sin(πy)
// ============================================================================

using namespace UNSTEADY::SINPI_SCALAR_2D;

static const PetscScalar alpha = 0.1;  // thermal diffusivity

// ============================================================================
// Set initial condition on element DOFs
// ============================================================================
PetscErrorCode SetInitialCondition(const DM &dm, Vec &u, Vec &uLocal) {
  PetscFunctionBeginUser;
  PetscScalar ***aU;
  PetscScalar **cX, **cY;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt icenter, ip;

  PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
  PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                             &nEx[0], &nEx[1], NULL));
  PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));

  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      const PetscScalar x = cX[ex][icenter];
      const PetscScalar y = cY[ey][icenter];
      aU[ey][ex][ip] = u_steady(x, y);
    }
  }

  PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
  PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(DMLocalToGlobal(dm, uLocal, INSERT_VALUES, u));
  PetscFunctionReturn(0);
}

// ============================================================================
// Compute residual norm for implicit Euler: r = u^n - (I - αΔt Δ_h) u^{n+1}
// ============================================================================
PetscErrorCode ComputeResidualNorm(const DM &dm,
                                   const Vec &u, const Vec &uLocal,
                                   const Vec &uOld, const Vec &uOldLocal,
                                   PetscInt Nx, PetscInt Ny,
                                   PetscReal dt, PetscReal *residualNorm) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  const PetscReal ix2 = 1.0 / (hx * hx);
  const PetscReal iy2 = 1.0 / (hy * hy);
  const PetscReal coef = alpha * dt;
  const PetscReal diag = 1.0 + 2.0 * coef * (ix2 + iy2);

  PetscScalar ***aU, ***aUold;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt icenter;

  PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMGlobalToLocalBegin(dm, uOld, INSERT_VALUES, uOldLocal));
  PetscCall(DMGlobalToLocalEnd(dm, uOld, INSERT_VALUES, uOldLocal));
  PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
  PetscCall(DMStagVecGetArray(dm, uOldLocal, &aUold));
  PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                             &nEx[0], &nEx[1], NULL));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &icenter));

  PetscReal localSum = 0.0;
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      const PetscScalar uc = aU[ey][ex][icenter];
      const PetscScalar uc_old = aUold[ey][ex][icenter];

      // Helmholtz stencil: (I - αΔt Δ_h) u^{n+1}
      const PetscScalar ul = DMStag_GetLeft(aU, ex, ey, icenter, Nx);
      const PetscScalar ur = DMStag_GetRight(aU, ex, ey, icenter, Nx);
      const PetscScalar ud = DMStag_GetDown(aU, ex, ey, icenter, Ny);
      const PetscScalar uu = DMStag_GetUp(aU, ex, ey, icenter, Ny);

      const PetscScalar lhs = diag * uc - coef * (ix2 * (ul + ur) + iy2 * (ud + uu));
      const PetscScalar r = uc_old - lhs;  // residual = u^n - (I - αΔt Δ_h) u^{n+1}
      localSum += PetscRealPart(r * r);
    }
  }

  PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
  PetscCall(DMStagVecRestoreArray(dm, uOldLocal, &aUold));

  PetscReal outNorm;
  MPI_Allreduce(&localSum, &outNorm, 1, MPIU_REAL, MPIU_SUM, PETSC_COMM_WORLD);
  *residualNorm = std::sqrt(outNorm);
  PetscFunctionReturn(0);
}

// ============================================================================
// Red-black Gauss-Seidel sweep for: (I - αΔt Δ_h) u^{n+1} = u^n
// ============================================================================
PetscErrorCode GaussSeidelSweep(const DM &dm, Vec &u, Vec &uLocal,
                                const Vec &uOld, const Vec &uOldLocal,
                                PetscInt Nx, PetscInt Ny, PetscReal dt) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  const PetscReal ix2 = 1.0 / (hx * hx);
  const PetscReal iy2 = 1.0 / (hy * hy);
  const PetscReal coef = alpha * dt;
  const PetscReal diag = 1.0 + 2.0 * coef * (ix2 + iy2);

  for (int color = 0; color < 2; ++color) {
    PetscScalar ***aU, ***aUold;
    PetscInt startx, starty, nx, ny, nEx[2];
    PetscInt icenter;

    PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
    PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
    PetscCall(DMGlobalToLocalBegin(dm, uOld, INSERT_VALUES, uOldLocal));
    PetscCall(DMGlobalToLocalEnd(dm, uOld, INSERT_VALUES, uOldLocal));
    PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
    PetscCall(DMStagVecGetArray(dm, uOldLocal, &aUold));
    PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                               &nEx[0], &nEx[1], NULL));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &icenter));

    for (PetscInt ey = starty; ey < starty + ny; ++ey) {
      for (PetscInt ex = startx; ex < startx + nx; ++ex) {
        if (((ex + ey) & 1) != color) continue;

        const PetscScalar uc_old = aUold[ey][ex][icenter];

        const PetscScalar ul = DMStag_GetLeft(aU, ex, ey, icenter, Nx);
        const PetscScalar ur = DMStag_GetRight(aU, ex, ey, icenter, Nx);
        const PetscScalar ud = DMStag_GetDown(aU, ex, ey, icenter, Ny);
        const PetscScalar uu = DMStag_GetUp(aU, ex, ey, icenter, Ny);

        // unew = (u_old + αΔt*(neighbor_sum)) / diag
        const PetscScalar unew = (uc_old + coef * (ix2 * (ul + ur) + iy2 * (ud + uu))) / diag;
        aU[ey][ex][icenter] = unew;
      }
    }

    PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
    PetscCall(DMStagVecRestoreArray(dm, uOldLocal, &aUold));
    PetscCall(DMLocalToGlobal(dm, uLocal, INSERT_VALUES, u));
  }
  PetscFunctionReturn(0);
}

// ============================================================================
// Compute L2 error at final time
// ============================================================================
PetscErrorCode ComputeError(const DM &dm, const Vec &u, Vec &uLocal,
                               PetscInt Nx, PetscInt Ny,
                               PetscReal t, PetscReal *error) {
  PetscFunctionBeginUser;
  Vec uExact, diff;
  PetscCall(DMCreateGlobalVector(dm, &uExact));
  PetscCall(DMGetLocalVector(dm, &uLocal));

  PetscScalar ***aUe;
  PetscScalar **cX, **cY;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt icenter, ip;

  PetscCall(DMGlobalToLocalBegin(dm, uExact, INSERT_VALUES, uLocal));
  PetscCall(DMGlobalToLocalEnd(dm, uExact, INSERT_VALUES, uLocal));
  PetscCall(DMStagVecGetArray(dm, uLocal, &aUe));
  PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                             &nEx[0], &nEx[1], NULL));
  PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));

  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      const PetscScalar x = cX[ex][icenter];
      const PetscScalar y = cY[ey][icenter];
      aUe[ey][ex][ip] = u_exact(x, y, t, alpha);
    }
  }

  PetscCall(DMStagVecRestoreArray(dm, uLocal, &aUe));
  PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(DMLocalToGlobal(dm, uLocal, INSERT_VALUES, uExact));

  PetscCall(VecDuplicate(u, &diff));
  PetscCall(VecCopy(u, diff));
  PetscCall(VecAXPY(diff, -1.0, uExact));

  PetscReal nrm2;
  PetscCall(VecNorm(diff, NORM_2, &nrm2));

  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  *error = nrm2 * std::sqrt(hx * hy);

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
                            "Cell-centered heat equation on DMStag"));

  DM dm;
  Vec u, uLocal, uOld, uOldLocal;
  PetscInt Nx = 64, Ny = 64;
  const PetscInt dof0 = 1, dof1 = 0, dof2 = 0;  // element DOFs only (scalar)
  const PetscInt stencilWidth = 1;
  PetscReal dt = 0.001;
  PetscReal T_final = 0.05;
  const PetscReal tol = 1e-8;
  const PetscInt maxIts = 20000;
  int rank;
  PetscBool compute_error = PETSC_FALSE;
  PetscBool convergence_test = PETSC_FALSE;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscCall(PetscOptionsGetInt(NULL, NULL, "-nx", &Nx, NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-ny", &Ny, NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-T", &T_final, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-heat_check_error", &compute_error, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-convergence_test", &convergence_test, NULL));

  // Time step: for convergence test, dt ∝ h²
  if (!convergence_test) {
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-dt", &dt, NULL));
  }
  if (convergence_test) {
    const PetscReal h = 1.0 / Nx;
    dt = h * h;
  }

  // Create DMStag (element-centered only)
  PetscCall(DMStagCreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                           Nx, Ny, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2,
                           DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, &dm));
  PetscCall(DMSetFromOptions(dm));
  PetscCall(DMSetUp(dm));
  PetscCall(DMStagSetUniformCoordinatesProduct(dm, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0));

  PetscCall(DMCreateGlobalVector(dm, &u));
  PetscCall(DMCreateGlobalVector(dm, &uOld));
  PetscCall(DMGetLocalVector(dm, &uLocal));
  PetscCall(DMGetLocalVector(dm, &uOldLocal));

  PetscInt Nglob[2];
  PetscCall(DMStagGetGlobalSizes(dm, &Nglob[0], &Nglob[1], NULL));

  const PetscInt Nsteps = static_cast<PetscInt>(std::ceil(T_final / dt));
  const PetscReal dtActual = T_final / Nsteps;

  if (rank == 0) {
    std::cout << "Cell-centered heat: N=" << Nglob[0] << "x" << Nglob[1]
              << ", h=" << 1.0 / Nglob[0]
              << ", dt=" << dtActual << ", steps=" << Nsteps
              << ", T=" << T_final << std::endl;
  }

  // Set initial condition
  PetscCall(SetInitialCondition(dm, u, uLocal));

  // Time stepping
  PetscReal t = 0.0;
  for (PetscInt step = 1; step <= Nsteps; ++step) {
    t += dtActual;

    // uOld = u^n
    PetscCall(VecCopy(u, uOld));

    // Gauss-Seidel solve: (I - αΔt Δ_h) u^{n+1} = u^n
    PetscReal resNorm = 0.0;
    PetscInt its = 0;
    for (its = 1; its <= maxIts; ++its) {
      PetscCall(GaussSeidelSweep(dm, u, uLocal, uOld, uOldLocal,
                                 Nglob[0], Nglob[1], dtActual));
      PetscCall(ComputeResidualNorm(dm, u, uLocal, uOld, uOldLocal,
                                    Nglob[0], Nglob[1], dtActual, &resNorm));
      if (resNorm <= tol) break;
    }

    if (rank == 0 && (step % 100 == 0 || step == Nsteps)) {
      std::cout << "  step " << step << "/" << Nsteps
                << ", t=" << t << ", GS its=" << its
                << ", ||r||=" << resNorm << std::endl;
    }
  }

  // Compute error
  PetscReal error = 0.0;
  if (compute_error) {
    PetscCall(ComputeError(dm, u, uLocal, Nglob[0], Nglob[1], t, &error));
    if (rank == 0) {
      std::cout << "||u(T) - u_exact(T)||_L2 = " << error << std::endl;
    }
  }

  if (convergence_test && rank == 0) {
    const PetscReal h = 1.0 / Nglob[0];
    std::cout << "CONVERGENCE: " << Nglob[0] << " " << h << " " << error
              << " " << Nsteps << std::endl;
  }

  // Cleanup
  PetscCall(DMRestoreLocalVector(dm, &uLocal));
  PetscCall(DMRestoreLocalVector(dm, &uOldLocal));
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&uOld));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFinalize());
  return 0;
}
