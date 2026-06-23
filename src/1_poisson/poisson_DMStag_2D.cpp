#include <petscdm.h>
#include <petscdmstag.h>
#include <petscksp.h>
#include <petscsys.h>

#include <cmath>
#include <iostream>
#include "analytical/unsteady.h"

static inline PetscScalar u_exact(PetscScalar x, PetscScalar y) {
  return std::sin(M_PI * x) * std::sin(M_PI * y)+x;
  // return UNSTEADY::TAYLOR_GREEN_2D::p_exact(x, y, 0.0);
}

static inline PetscScalar rhs_f(PetscScalar x, PetscScalar y) {
  return 2.0 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
  // return UNSTEADY::TAYLOR_GREEN_2D::f_poisson(x, y, 0.0);
}

// ============================================================================
// Setup right-hand side vector
// ============================================================================
PetscErrorCode SetupRHS(const DM &dm, Vec &f, Vec &fLocal,
                        PetscScalar (*rhs_func)(PetscScalar, PetscScalar)) {
  PetscFunctionBeginUser;
  // Fill RHS f at element centers
  {
    PetscScalar ***aF;
    PetscScalar **cX, **cY;
    PetscInt startx, starty, nx, ny, nEx[2];
    PetscInt icenter, ip;
    PetscCall(DMGlobalToLocalBegin(dm, f, INSERT_VALUES, fLocal));
    PetscCall(DMGlobalToLocalEnd(dm, f, INSERT_VALUES, fLocal));
    PetscCall(DMStagVecGetArray(dm, fLocal, &aF));
    PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                               &nEx[0], &nEx[1], NULL));
    PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
    PetscCall(
        DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));
    for (PetscInt ey = starty; ey < starty + ny; ++ey) {
      for (PetscInt ex = startx; ex < startx + nx; ++ex) {
        const PetscScalar x = cX[ex][icenter];
        const PetscScalar y = cY[ey][icenter];
        aF[ey][ex][ip] = rhs_f(x, y);
      }
    }
    PetscCall(DMStagVecRestoreArray(dm, fLocal, &aF));
    PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
    PetscCall(DMLocalToGlobal(dm, fLocal, INSERT_VALUES, f));
  }
  PetscFunctionReturn(0);
}

// ============================================================================
// Compute residual norm
// ============================================================================
PetscErrorCode ComputeResidualNorm(const DM &dm, const Vec &u,
                                   const Vec &uLocal, const Vec &f,
                                   const Vec &fLocal, PetscInt Nx, PetscInt Ny,
                                   PetscReal *residualNorm) {
  PetscFunctionBeginUser;
  PetscScalar ***aU, ***aF;
  PetscScalar **cX, **cY;
  PetscScalar outNorm;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt icenter, ip;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  const PetscReal ix2 = 1.0 / (hx * hx);
  const PetscReal iy2 = 1.0 / (hy * hy);
  const PetscReal diag = 2.0 * ix2 + 2.0 * iy2;

  PetscCallAbort(PETSC_COMM_WORLD,
                 DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
  PetscCallAbort(PETSC_COMM_WORLD,
                 DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
  PetscCallAbort(PETSC_COMM_WORLD,
                 DMGlobalToLocalBegin(dm, f, INSERT_VALUES, fLocal));
  PetscCallAbort(PETSC_COMM_WORLD,
                 DMGlobalToLocalEnd(dm, f, INSERT_VALUES, fLocal));
  PetscCallAbort(PETSC_COMM_WORLD, DMStagVecGetArray(dm, uLocal, &aU));
  PetscCallAbort(PETSC_COMM_WORLD, DMStagVecGetArray(dm, fLocal, &aF));
  PetscCallAbort(PETSC_COMM_WORLD,
                 DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                                  &nEx[0], &nEx[1], NULL));
  PetscCallAbort(PETSC_COMM_WORLD,
                 DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
  PetscCallAbort(PETSC_COMM_WORLD,
                 DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));
  PetscCallAbort(PETSC_COMM_WORLD, DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));

  PetscReal localSum = 0.0;
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      // For cell-centered grid, all cells are interior
      const PetscScalar uc = aU[ey][ex][ip];
      const PetscScalar x = cX[ex][icenter];
      const PetscScalar y = cY[ey][icenter];

      // Get neighbor values with ghost cell BC
      PetscScalar ul, ur, ud, uu;
      if (ex == 0) {
        PetscScalar u_bnd = u_exact(0.0, y);
        ul = 2.0 * u_bnd - uc;  // This is correct for cell-centered
      } else {
        ul = aU[ey][ex - 1][ip];
      }

      if (ex == Nx - 1) {
        PetscScalar u_bnd = u_exact(1.0, y);
        ur = 2.0 * u_bnd - uc;
      } else {
        ur = aU[ey][ex + 1][ip];
      }

      if (ey == 0) {
        PetscScalar u_bnd = u_exact(x, 0.0);
        ud = 2.0 * u_bnd - uc;
      } else {
        ud = aU[ey - 1][ex][ip];
      }

      if (ey == Ny - 1) {
        PetscScalar u_bnd = u_exact(x, 1.0);
        uu = 2.0 * u_bnd - uc;
      } else {
        uu = aU[ey + 1][ex][ip];
      }

      const PetscScalar Au = diag * uc - ix2 * (ul + ur) - iy2 * (ud + uu);
      const PetscScalar ff = aF[ey][ex][ip];
      const PetscScalar r = ff - Au; // residual for -Laplace(u)=f
      localSum += PetscRealPart(r * r);
    }
  }
  PetscCallAbort(PETSC_COMM_WORLD, DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCallAbort(PETSC_COMM_WORLD, DMStagVecRestoreArray(dm, uLocal, &aU));
  PetscCallAbort(PETSC_COMM_WORLD, DMStagVecRestoreArray(dm, fLocal, &aF));
  MPI_Allreduce(&localSum, &outNorm, 1, MPIU_REAL, MPIU_SUM, PETSC_COMM_WORLD);
  *residualNorm = std::sqrt(outNorm);
  PetscFunctionReturn(0);
}

// ============================================================================
// Perform one red-black Gauss-Seidel sweep
// ============================================================================
PetscErrorCode GaussSeidelSweep(const DM &dm, Vec &u, Vec &uLocal, const Vec &f,
                                const Vec &fLocal, PetscInt Nx, PetscInt Ny) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  const PetscReal ix2 = 1.0 / (hx * hx);
  const PetscReal iy2 = 1.0 / (hy * hy);
  const PetscReal diag = 2.0 * ix2 + 2.0 * iy2;
  for (int color = 0; color < 2; ++color) {
    PetscScalar ***aU, ***aF;
    PetscScalar **cX, **cY;
    PetscInt startx, starty, nx, ny, nEx[2];
    PetscInt icenter, ip;
    PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
    PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
    PetscCall(DMGlobalToLocalBegin(dm, f, INSERT_VALUES, fLocal));
    PetscCall(DMGlobalToLocalEnd(dm, f, INSERT_VALUES, fLocal));
    PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
    PetscCall(DMStagVecGetArray(dm, fLocal, &aF));
    PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                               &nEx[0], &nEx[1], NULL));
    PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));
    PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));

    for (PetscInt ey = starty; ey < starty + ny; ++ey) {
      for (PetscInt ex = startx; ex < startx + nx; ++ex) {
        if (((ex + ey) & 1) != color)
          continue;

        PetscScalar ul, ur, ud, uu;
        const PetscScalar x = cX[ex][icenter];
        const PetscScalar y = cY[ey][icenter];
        if (ex == 0) {
        PetscScalar u_bnd = u_exact(0.0, y);
          ul = 2*u_bnd-aU[ey][ex][ip]; // Ghost cell: u_{-1} = -u_0 for Dirichlet BC u=0
        } else {
          ul = aU[ey][ex - 1][ip];
        }

        if (ex == Nx - 1) {
        PetscScalar u_bnd = u_exact(1.0, y);
          ur = 2*u_bnd-aU[ey][ex][ip]; // Ghost cell: u_{N} = -u_{N-1}
        } else {
          ur = aU[ey][ex + 1][ip];
        }

        if (ey == 0) {
        PetscScalar u_bnd = u_exact(x, 0.0);
          ud = 2*u_bnd-aU[ey][ex][ip]; // Ghost cell: u_{-1} = -u_0
        } else {
          ud = aU[ey - 1][ex][ip];
        }

        if (ey == Ny - 1) {
        PetscScalar u_bnd = u_exact(x, 1.0);
          uu = 2*u_bnd-aU[ey][ex][ip]; // Ghost cell: u_{N} = -u_{N-1}
        } else {
          uu = aU[ey + 1][ex][ip];
        }

        const PetscScalar ff = aF[ey][ex][ip];
        const PetscScalar unew =
            (ix2 * (ul + ur) + iy2 * (ud + uu) + ff) / diag;
        aU[ey][ex][ip] = unew;
      }
    }

    PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
    PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
    PetscCall(DMStagVecRestoreArray(dm, fLocal, &aF));
    // scatter updated owned values back to global, then refresh ghosts for
    // next color
    PetscCall(DMLocalToGlobal(dm, uLocal, INSERT_VALUES, u));
  }
  PetscFunctionReturn(0);
}

// ============================================================================
// Compute L2 error against exact solution
// ============================================================================
PetscErrorCode ComputeL2Error(const DM &dm, const Vec &u, Vec &uLocal,
                              PetscScalar (*exact_func)(PetscScalar,
                                                        PetscScalar),
                              PetscInt Nx, PetscInt Ny, PetscReal *l2Error) {
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
  PetscCall(
      DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));

  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      const PetscScalar x = cX[ex][icenter];
      const PetscScalar y = cY[ey][icenter];
      aUe[ey][ex][ip] = exact_func(x, y);
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
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  *l2Error = nrm2 * std::sqrt(hx * hy);
  if (rank == 0)
    std::cout << "||u - u_exact||_2 = " << nrm2 * std::sqrt(hx * hy)
              << std::endl;
  PetscCall(VecDestroy(&diff));
  PetscCall(VecDestroy(&uExact));
  PetscFunctionReturn(0);
}

int main(int argc, char **argv) {
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, "Matrix-free Gauss-Seidel for Poisson on DMStag (element-centered)"));

  DM dm;
  Vec u, uLocal, f, fLocal;
  PetscInt Nx = 128, Ny = 128;
  const PetscInt dof0 = 0, dof1 = 0, dof2 = 1;
  const PetscInt stencilWidth = 1;
  const PetscReal tol = 1e-8;
  const PetscInt maxIts = 40000;
  PetscInt its = 0;
  PetscReal resNorm = 0.0;
  PetscReal res0 = -1.0;
  const PetscReal atol = tol;
  int rank;
  PetscBool compute_error = PETSC_FALSE;
  PetscBool convergence_test = PETSC_FALSE;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  // Get options from command line
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-nx", &Nx, NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-ny", &Ny, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-poisson_check_error", &compute_error, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-convergence_test", &convergence_test, NULL));

  PetscCall(DMStagCreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                           Nx, Ny, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2,
                           DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, &dm));
  PetscCall(DMSetFromOptions(dm));
  PetscCall(DMSetUp(dm));
  PetscCall(DMStagSetUniformCoordinatesProduct(dm, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0));
  PetscCall(DMCreateGlobalVector(dm, &u));
  PetscCall(DMCreateGlobalVector(dm, &f));
  PetscCall(DMGetLocalVector(dm, &uLocal));
  PetscCall(DMGetLocalVector(dm, &fLocal));
  PetscCall(VecSet(u, 0.0));

  // 1. Setup RHS
  PetscCall(SetupRHS(dm, f, fLocal, rhs_f));

  PetscInt Nglob[2];
  PetscCall(DMStagGetGlobalSizes(dm, &Nglob[0], &Nglob[1], NULL));

  // compute_residual_norm(resNorm);
  ComputeResidualNorm(dm, u, uLocal, f, fLocal, Nglob[0], Nglob[1], &resNorm);
  res0 = resNorm;
  if (rank == 0)
    std::cout << "GS start: ||r|| = " << resNorm << std::endl;

  // Red-black Gauss-Seidel iterations
  for (its = 1; its <= maxIts; ++its) {
    PetscCall(GaussSeidelSweep(dm, u, uLocal, f, fLocal, Nglob[0], Nglob[1]));
    PetscCall(ComputeResidualNorm(dm, u, uLocal, f, fLocal, Nglob[0], Nglob[1],
                                  &resNorm));
    if (rank == 0 && (its % 5000 == 0 || resNorm <= atol)) {
      std::cout << "it=" << its << ", ||r||=" << resNorm << " (rel "
                << resNorm / (res0 > 0 ? res0 : 1.0) << ")" << std::endl;
    }
    if (resNorm <= atol)
      break;
  }

  if (rank == 0) {
    std::cout << "Done. it=" << its << ", ||r||=" << resNorm << std::endl;
  }

  PetscReal l2Error = 0.0;
  if (compute_error) {
    PetscCall(ComputeL2Error(dm, u, uLocal, u_exact, Nglob[0], Nglob[1], &l2Error));
  }
  
  // Convergence test mode: print error in parseable format
  if (convergence_test && rank == 0) {
    const PetscReal hx = 1.0 / Nglob[0];
    const PetscReal hy = 1.0 / Nglob[1];
    const PetscReal h = std::sqrt(hx * hy);
    std::cout << "CONVERGENCE: " << Nglob[0] << " " << h << " " << l2Error << " " << its << std::endl;
  }

  PetscCall(DMRestoreLocalVector(dm, &uLocal));
  PetscCall(DMRestoreLocalVector(dm, &fLocal));
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&f));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFinalize());
  return 0;
}

// time mpirun -np 4 ./ex_poisson_stagger -poisson_check_error
// GS start: ||r|| = 1263.31
// it=5000, ||r||=87.8774 (rel 0.0695613)
// it=10000, ||r||=4.32179 (rel 0.00342101)
// it=15000, ||r||=0.212544 (rel 0.000168244)
// it=20000, ||r||=0.0104529 (rel 8.2742e-06)
// it=25000, ||r||=0.00051407 (rel 4.06923e-07)
// it=30000, ||r||=2.52818e-05 (rel 2.00124e-08)
// it=35000, ||r||=1.24335e-06 (rel 9.84201e-10)
// it=40000, ||r||=6.11473e-08 (rel 4.84025e-11)
// Done. it=40001, ||r||=6.11473e-08
// ||u - u_exact||_2 = 2.51004e-05
// mpirun -np 4 ./ex_poisson_stagger -poisson_check_error  21.04s user 0.37s
// system 395% cpu 5.417 total

// GS start: ||r|| = 1263.31
// it=5000, ||r||=87.8774 (rel 0.0695613)
// it=10000, ||r||=4.32179 (rel 0.00342101)
// it=15000, ||r||=0.212544 (rel 0.000168244)
// it=20000, ||r||=0.0104529 (rel 8.2742e-06)
// it=25000, ||r||=0.00051407 (rel 4.06923e-07)
// it=30000, ||r||=2.52818e-05 (rel 2.00124e-08)
// it=35000, ||r||=1.24335e-06 (rel 9.84201e-10)
// it=40000, ||r||=6.11473e-08 (rel 4.84025e-11)
// Done. it=40001, ||r||=6.11473e-08
// ||u - u_exact||_2 = 2.51004e-05
// mpirun -np 1 ./ex_poisson_stagger -poisson_check_error  16.33s user 0.20s
// system 94% cpu 17.433 total