#include <cmath>
#include <iostream>
#include <petscdm.h>
#include <petscdmstag.h>
#include <petscksp.h>
#include <petscsys.h>
#include <string>

#include "manufactured_solutions.h"

// ============================================================================
// Compute velocity divergence from edge DOFs and store in element centers
// This is used as RHS for pressure Poisson equation
// ============================================================================
PetscErrorCode ComputeDivergence(const DM &dm, const Vec &u, const Vec &uLocal,
                                 Vec &div, Vec &divLocal, PetscInt Nx,
                                 PetscInt Ny) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;

  PetscScalar ***aU, ***aDiv;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt iux, iuy, ip;

  // Update ghost cells
  PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMGlobalToLocalBegin(dm, div, INSERT_VALUES, divLocal));
  PetscCall(DMGlobalToLocalEnd(dm, div, INSERT_VALUES, divLocal));

  PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
  PetscCall(DMStagVecGetArray(dm, divLocal, &aDiv));
  PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                             &nEx[0], &nEx[1], NULL));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));

  // Compute divergence at each cell center
  // div(u) = du/dx + dv/dy
  // Using finite differences:
  // du/dx ≈ (u_right - u_left) / hx
  // dv/dy ≈ (v_up - v_down) / hy
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      // u component on LEFT edges: ex and ex+1
      const PetscScalar u_left = aU[ey][ex][iux];      // LEFT edge at ex
      const PetscScalar u_right = aU[ey][ex + 1][iux]; // LEFT edge at ex+1
      const PetscScalar dudx = (u_right - u_left) / hx;

      // v component on DOWN edges: ey and ey+1
      const PetscScalar v_down = aU[ey][ex][iuy];   // DOWN edge at ey
      const PetscScalar v_up = aU[ey + 1][ex][iuy]; // DOWN edge at ey+1
      const PetscScalar dvdy = (v_up - v_down) / hy;

      // Divergence at cell center
      aDiv[ey][ex][ip] = dudx + dvdy;
    }
  }

  PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
  PetscCall(DMStagVecRestoreArray(dm, divLocal, &aDiv));
  PetscCall(DMLocalToGlobal(dm, divLocal, INSERT_VALUES, div));
  PetscFunctionReturn(0);
}

// ============================================================================
// Velocity update using pressure gradient (projection step)
// u^{n+1} = u* - dt * grad(p)
// - Updates edge-centered velocities (LEFT -> u_x, DOWN -> u_y)
// - Reads pressure from element-centered slot
// - Re-enforces boundary edges (remain unchanged as we skip boundaries)
// ============================================================================
PetscErrorCode UpdateVelocity(const DM &dm, Vec &u, Vec &uLocal, PetscReal dt) {
  PetscFunctionBeginUser;
  // Grid spacing from global sizes
  PetscInt Nx, Ny;
  PetscCall(DMStagGetGlobalSizes(dm, &Nx, &Ny, NULL));
  const PetscReal hx = 1.0 / static_cast<PetscReal>(Nx);
  const PetscReal hy = 1.0 / static_cast<PetscReal>(Ny);

  // Access arrays and slots
  PetscScalar ***aU;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt iux, iuy, ip;

  // Refresh ghosts for u
  PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));

  PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
  PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                             &nEx[0], &nEx[1], NULL));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));

  // Correct u_y on DOWN edges (interior only)
  for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
    if (ey == 0 || ey == Ny) continue; // skip boundary edges
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      // Pressure in cells below and above the edge
      const PetscScalar pD = aU[ey - 1][ex][ip];
      const PetscScalar pU = aU[ey][ex][ip];
      const PetscScalar dpdy = (pU - pD) / hy;
      aU[ey][ex][iuy] -= dt * dpdy;
    }
  }

  // Correct u_x on LEFT edges (interior only)
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx + nEx[0]; ++ex) {
      if (ex == 0 || ex == Nx) continue; // skip boundary edges
      // Pressure in cells left and right of the edge
      const PetscScalar pL = aU[ey][ex - 1][ip];
      const PetscScalar pR = aU[ey][ex][ip];
      const PetscScalar dpdx = (pR - pL) / hx;
      aU[ey][ex][iux] -= dt * dpdx;
    }
  }

  PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
  // Scatter updated velocities back to global
  PetscCall(DMLocalToGlobal(dm, uLocal, INSERT_VALUES, u));
  PetscFunctionReturn(0);
}

// ============================================================================
// Setup RHS for Poisson using velocity divergence
// ============================================================================
PetscErrorCode SetupRHS_Poisson_Divergence(const DM &dm, Vec &f, Vec &fLocal,
                                           const Vec &u, const Vec &uLocal,
                                           PetscInt Nx, PetscInt Ny) {
  PetscFunctionBeginUser;
  PetscCall(ComputeDivergence(dm, u, uLocal, f, fLocal, Nx, Ny));
  PetscFunctionReturn(0);
}

// ============================================================================
// Setup RHS for Heat (edge-centered, time-dependent)
// ============================================================================
PetscErrorCode SetupRHS_Heat(const DM &dm, Vec &f, Vec &fLocal, PetscReal t) {
  PetscFunctionBeginUser;
  using namespace HEAT::DECAY_2D;

  PetscScalar ***aF;
  PetscScalar **cX, **cY;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt iprev, icenter, iux, iuy;

  PetscCall(DMGlobalToLocalBegin(dm, f, INSERT_VALUES, fLocal));
  PetscCall(DMGlobalToLocalEnd(dm, f, INSERT_VALUES, fLocal));
  PetscCall(DMStagVecGetArray(dm, fLocal, &aF));
  PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                             &nEx[0], &nEx[1], NULL));
  PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(
      DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_LEFT, &iprev));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));

  for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx + nEx[0]; ++ex) {
      aF[ey][ex][iuy] = rhs_f(cX[ex][icenter], cY[ey][iprev], t);
      aF[ey][ex][iux] = rhs_f(cX[ex][iprev], cY[ey][icenter], t);
    }
  }

  PetscCall(DMStagVecRestoreArray(dm, fLocal, &aF));
  PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(DMLocalToGlobal(dm, fLocal, INSERT_VALUES, f));
  PetscFunctionReturn(0);
}

// ============================================================================
// Setup initial condition for Heat
// ============================================================================
PetscErrorCode SetupInitialCondition_Heat(const DM &dm, Vec &u, Vec &uLocal) {
  PetscFunctionBeginUser;
  using namespace HEAT::DECAY_2D;

  PetscScalar ***aU;
  PetscScalar **cX, **cY;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt iprev, icenter, iux, iuy;

  PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
  PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                             &nEx[0], &nEx[1], NULL));
  PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(
      DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_LEFT, &iprev));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));

  for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx + nEx[0]; ++ex) {
      aU[ey][ex][iuy] = u_initial(cX[ex][icenter], cY[ey][iprev]);
      aU[ey][ex][iux] = u_initial(cX[ex][iprev], cY[ey][icenter]);
    }
  }

  PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
  PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(DMLocalToGlobal(dm, uLocal, INSERT_VALUES, u));
  PetscFunctionReturn(0);
}

// ============================================================================
// Gauss-Seidel sweep for Poisson (element-centered)
// ============================================================================
PetscErrorCode GaussSeidelSweep_Poisson(const DM &dm, Vec &u, Vec &uLocal,
                                        const Vec &f, const Vec &fLocal,
                                        PetscInt Nx, PetscInt Ny) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  const PetscReal ix2 = 1.0 / (hx * hx);
  const PetscReal iy2 = 1.0 / (hy * hy);
  const PetscReal diag = 2.0 * (ix2 + iy2);

  for (int color = 0; color < 2; ++color) {
    PetscScalar ***aU, ***aF;
    PetscInt startx, starty, nx, ny, nEx[2];
    PetscInt icenter;

    PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
    PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
    PetscCall(DMGlobalToLocalBegin(dm, f, INSERT_VALUES, fLocal));
    PetscCall(DMGlobalToLocalEnd(dm, f, INSERT_VALUES, fLocal));
    PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
    PetscCall(DMStagVecGetArray(dm, fLocal, &aF));
    PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                               &nEx[0], &nEx[1], NULL));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &icenter));

    for (PetscInt ey = starty; ey < starty + ny; ++ey) {
      for (PetscInt ex = startx; ex < startx + nx; ++ex) {
        if (((ex + ey) & 1) != color)
          continue;

        const PetscScalar uc = aU[ey][ex][icenter];
        PetscScalar ul = (ex == 0) ? -uc : aU[ey][ex - 1][icenter];
        PetscScalar ur = (ex == Nx - 1) ? -uc : aU[ey][ex + 1][icenter];
        PetscScalar ud = (ey == 0) ? -uc : aU[ey - 1][ex][icenter];
        PetscScalar uu = (ey == Ny - 1) ? -uc : aU[ey + 1][ex][icenter];

        const PetscScalar ff = aF[ey][ex][icenter];
        aU[ey][ex][icenter] = (ix2 * (ul + ur) + iy2 * (ud + uu) + ff) / diag;
      }
    }

    PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
    PetscCall(DMStagVecRestoreArray(dm, fLocal, &aF));
    PetscCall(DMLocalToGlobal(dm, uLocal, INSERT_VALUES, u));
  }
  PetscFunctionReturn(0);
}

// ============================================================================
// Gauss-Seidel sweep for Heat (edge-centered)
// ============================================================================
PetscErrorCode GaussSeidelSweep_Heat(const DM &dm, Vec &u, Vec &uLocal,
                                     const Vec &uOld, const Vec &uOldLocal,
                                     const Vec &f, const Vec &fLocal,
                                     PetscInt Nx, PetscInt Ny, PetscReal alpha,
                                     PetscReal dt) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  const PetscReal ix2 = 1.0 / (hx * hx);
  const PetscReal iy2 = 1.0 / (hy * hy);
  const PetscReal coef = alpha * dt;
  const PetscReal diag = 1.0 + 2.0 * coef * (ix2 + iy2);

  for (int color = 0; color < 2; ++color) {
    PetscScalar ***aU, ***aUold, ***aF;
    PetscInt startx, starty, nx, ny, nEx[2];
    PetscInt iux, iuy;

    PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
    PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
    PetscCall(DMGlobalToLocalBegin(dm, uOld, INSERT_VALUES, uOldLocal));
    PetscCall(DMGlobalToLocalEnd(dm, uOld, INSERT_VALUES, uOldLocal));
    PetscCall(DMGlobalToLocalBegin(dm, f, INSERT_VALUES, fLocal));
    PetscCall(DMGlobalToLocalEnd(dm, f, INSERT_VALUES, fLocal));
    PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
    PetscCall(DMStagVecGetArray(dm, uOldLocal, &aUold));
    PetscCall(DMStagVecGetArray(dm, fLocal, &aF));
    PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                               &nEx[0], &nEx[1], NULL));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));

    // Update DOWN edges
    for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
      if (ey == 0 || ey == Ny)
        continue;

      for (PetscInt ex = startx; ex < startx + nx; ++ex) {
        if (((ex + ey) & 1) != color)
          continue;

        const PetscScalar uc = aU[ey][ex][iuy];
        PetscScalar ul = (ex == 0) ? -uc : aU[ey][ex - 1][iuy];
        PetscScalar ur = (ex == Nx - 1) ? -uc : aU[ey][ex + 1][iuy];
        PetscScalar ud = (ey == 1) ? 0.0 : aU[ey - 1][ex][iuy];
        PetscScalar uu = (ey == Ny - 1) ? 0.0 : aU[ey + 1][ex][iuy];

        const PetscScalar rhs = aUold[ey][ex][iuy] + dt * aF[ey][ex][iuy];
        aU[ey][ex][iuy] =
            (rhs + coef * (ix2 * (ul + ur) + iy2 * (ud + uu))) / diag;
      }
    }

    // Update LEFT edges
    for (PetscInt ey = starty; ey < starty + ny; ++ey) {
      for (PetscInt ex = startx; ex < startx + nx + nEx[0]; ++ex) {
        if (ex == 0 || ex == Nx)
          continue;
        if (((ex + ey) & 1) != color)
          continue;

        const PetscScalar uc = aU[ey][ex][iux];
        PetscScalar ul = (ex == 1) ? 0.0 : aU[ey][ex - 1][iux];
        PetscScalar ur = (ex == Nx - 1) ? 0.0 : aU[ey][ex + 1][iux];
        PetscScalar ud = (ey == 0) ? -uc : aU[ey - 1][ex][iux];
        PetscScalar uu = (ey == Ny - 1) ? -uc : aU[ey + 1][ex][iux];

        const PetscScalar rhs = aUold[ey][ex][iux] + dt * aF[ey][ex][iux];
        aU[ey][ex][iux] =
            (rhs + coef * (ix2 * (ul + ur) + iy2 * (ud + uu))) / diag;
      }
    }

    PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
    PetscCall(DMStagVecRestoreArray(dm, uOldLocal, &aUold));
    PetscCall(DMStagVecRestoreArray(dm, fLocal, &aF));
    PetscCall(DMLocalToGlobal(dm, uLocal, INSERT_VALUES, u));
  }
  PetscFunctionReturn(0);
}

// ============================================================================
// Compute L2 error for Heat
// ============================================================================
PetscErrorCode ComputeL2Error_Heat(const DM &dm, const Vec &u, Vec &uLocal,
                                   PetscInt Nx, PetscInt Ny, PetscReal t,
                                   PetscReal alpha, PetscReal *l2Error) {
  PetscFunctionBeginUser;
  using namespace HEAT::DECAY_2D;

  PetscScalar ***aU;
  PetscScalar **cX, **cY;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt iprev, icenter, iux, iuy, ip;

  PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
  PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                             &nEx[0], &nEx[1], NULL));
  PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(
      DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_LEFT, &iprev));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));

  PetscReal localError2 = 0.0;

  // DOWN edges (interior)
  for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
    if (ey == 0 || ey == Ny)
      continue;
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      const PetscScalar x = cX[ex][icenter];
      const PetscScalar y = cY[ey][iprev];
      const PetscScalar diff = aU[ey][ex][iuy] - u_exact(x, y, t, alpha);
      localError2 += diff * diff;
    }
  }

  // LEFT edges (interior)
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx + nEx[0]; ++ex) {
      if (ex == 0 || ex == Nx)
        continue;
      const PetscScalar x = cX[ex][iprev];
      const PetscScalar y = cY[ey][icenter];
      const PetscScalar diff = aU[ey][ex][iux] - u_exact(x, y, t, alpha);
      localError2 += diff * diff;
    }
  }

  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      const PetscScalar x = cX[ex][icenter];
      const PetscScalar y = cY[ey][icenter];
      const PetscScalar diff = aU[ey][ex][ip] - POISSON::TRIG_2D::u_exact(x, y);
      localError2 += diff * diff;
    }
  }

  PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
  PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));

  PetscReal globalError2;
  PetscCall(MPI_Allreduce(&localError2, &globalError2, 1, MPIU_REAL, MPI_SUM,
                          PetscObjectComm((PetscObject)dm)));

  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  *l2Error = std::sqrt(globalError2 * hx * hy);
  PetscFunctionReturn(0);
}

// ============================================================================
// Main
// ============================================================================
int main(int argc, char **argv) {
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(
      &argc, &argv, NULL, "Unified PDE solver for Poisson and Heat equations"));

  // Common parameters
  DM dm;
  Vec u, uLocal, f, fLocal, uOld, uOldLocal;
  PetscInt Nx = 32, Ny = 32;
  const PetscInt stencilWidth = 1;
  const PetscReal tol = 1e-8;
  const PetscInt maxIts = 40000;
  int rank;
  PetscBool compute_error = PETSC_FALSE;
  PetscBool convergence_test = PETSC_FALSE;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscCall(PetscOptionsGetInt(NULL, NULL, "-nx", &Nx, NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-ny", &Ny, NULL));
  PetscCall(
      PetscOptionsGetBool(NULL, NULL, "-check_error", &compute_error, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-convergence_test",
                                &convergence_test, NULL));

  // Problem-specific parameters
  PetscInt dof0 = 0, dof1 = 1, dof2 = 1;
  PetscReal alpha = 0.1, dt = 0.001, T_final = 0.1;

  // Create DMStag
  PetscCall(DMStagCreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                           Nx, Ny, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2,
                           DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, &dm));
  PetscCall(DMSetFromOptions(dm));
  PetscCall(DMSetUp(dm));
  PetscCall(
      DMStagSetUniformCoordinatesProduct(dm, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0));

  // Create vectors
  PetscCall(DMCreateGlobalVector(dm, &u));
  PetscCall(DMCreateGlobalVector(dm, &f));
  PetscCall(DMGetLocalVector(dm, &uLocal));
  PetscCall(DMGetLocalVector(dm, &fLocal));
  PetscCall(VecSet(u, 0.0));
  PetscCall(DMCreateGlobalVector(dm, &uOld));
  PetscCall(DMGetLocalVector(dm, &uOldLocal));

  PetscInt Nglob[2];
  PetscCall(DMStagGetGlobalSizes(dm, &Nglob[0], &Nglob[1], NULL));

  if (rank == 0) {
    std::cout << "Grid: " << Nglob[0] << " x " << Nglob[1] << std::endl;
  }

  PetscCall(SetupInitialCondition_Heat(dm, u, uLocal));
  PetscReal t = 0.0;
  PetscInt step = 0;
  const PetscInt output_interval = std::max(1, (int)(T_final / dt / 10));

  // Time-stepping loop
  while (t < T_final) {
    step++;
    t += dt;

    // Tentitive velocity update
    PetscCall(VecCopy(u, uOld));
    PetscCall(SetupRHS_Heat(dm, f, fLocal, t));

    for (PetscInt its = 1; its <= maxIts; ++its) {
      PetscCall(GaussSeidelSweep_Heat(dm, u, uLocal, uOld, uOldLocal, f, fLocal,
                                      Nglob[0], Nglob[1], alpha, dt));
    }

    if (rank == 0 && (step % output_interval == 0 || step == 1)) {
      std::cout << "Step " << step << ", t=" << t << std::endl;
    }

    // Solve Poisson equation
    PetscCall(SetupRHS_Poisson_Divergence(dm, f, fLocal, u, uLocal, Nglob[0],
                                          Nglob[1]));

    for (PetscInt its = 1; its <= maxIts; ++its) {
      PetscCall(GaussSeidelSweep_Poisson(dm, u, uLocal, f, fLocal, Nglob[0],
                                         Nglob[1]));

      if (its % 5000 == 0 && rank == 0) {
        std::cout << "Iteration " << its << std::endl;
      }
    }

    // Update Velocity using pressure gradient (projection step)
    PetscCall(UpdateVelocity(dm, u, uLocal, dt));
  }

  if (rank == 0)
    std::cout << "\nTime integration complete. Final time: " << t << std::endl;

  if (compute_error || convergence_test) {
    PetscReal l2Error;
    PetscCall(ComputeL2Error_Heat(dm, u, uLocal, Nglob[0], Nglob[1], t, alpha,
                                  &l2Error));
    if (rank == 0) {
      std::cout << "||u(T) - u_exact(T)||_L2 = " << l2Error << std::endl;
      if (convergence_test) {
        const PetscReal h = std::sqrt(1.0 / (Nglob[0] * Nglob[1]));
        std::cout << "CONVERGENCE: " << Nglob[0] << " " << h << " " << l2Error
                  << std::endl;
      }
    }
  }

  // Cleanup
  PetscCall(DMRestoreLocalVector(dm, &uLocal));
  PetscCall(DMRestoreLocalVector(dm, &fLocal));
  PetscCall(DMRestoreLocalVector(dm, &uOldLocal));
  PetscCall(VecDestroy(&uOld));
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&f));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFinalize());
  return 0;
}
