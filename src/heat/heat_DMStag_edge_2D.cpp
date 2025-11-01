#include <petscdm.h>
#include <petscdmstag.h>
#include <petscksp.h>
#include <petscsys.h>
#include <cmath>
#include <iostream>

#include "manufactured_solutions.h"

// Use manufactured solution from library with non-zero Dirichlet BC
using namespace HEAT::NONZERO_DIRICHLET_2D;
namespace FUNC = HEAT::FUNC_2D;


// ============================================================================
// Compute velocity divergence from edge DOFs and store in element centers
// This is used as RHS for pressure Poisson equation in Stokes/NS solvers
// 
// For the heat equation: Although we're solving a scalar PDE, the edge DOFs 
// (u_x on LEFT edges, u_y on DOWN edges) can be interpreted as components
// of a vector field. The divergence measures the "compressibility" of this field.
// For the exact solution, div should match: ∂²u/∂x² + ∂²u/∂y² (related to Laplacian)
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
// Setup right-hand side vector
// ============================================================================
PetscErrorCode SetupRHS(const DM &dm, Vec &f, Vec &fLocal,
                        PetscScalar t, FUNC::RHS _rhs_f) {
  PetscFunctionBeginUser;
  PetscScalar ***aF;
  PetscScalar **cX, **cY;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt iprev, icenter, ip, iux, iuy;
  PetscCall(DMGlobalToLocalBegin(dm, f, INSERT_VALUES, fLocal));
  PetscCall(DMGlobalToLocalEnd(dm, f, INSERT_VALUES, fLocal));
  PetscCall(DMStagVecGetArray(dm, fLocal, &aF));
  PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));

  PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_LEFT, &iprev));
  
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy)); 
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));

  for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx + nEx[0]; ++ex) {
      aF[ey][ex][iuy] = _rhs_f(cX[ex][icenter], cY[ey][iprev], t);
      aF[ey][ex][iux] = _rhs_f(cX[ex][iprev], cY[ey][icenter], t);
      if (ey < starty + ny && ex < startx + nx) {
        aF[ey][ex][ip] = _rhs_f(cX[ex][icenter], cY[ey][icenter], t);
      }
    }
  }
  PetscCall(DMStagVecRestoreArray(dm, fLocal, &aF));
  PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(DMLocalToGlobal(dm, fLocal, INSERT_VALUES, f));
  PetscFunctionReturn(0);
}

// ============================================================================
// Setup initial condition
// ============================================================================
PetscErrorCode SetupInitialCondition(const DM &dm, Vec &u, Vec &uLocal, FUNC::INITIAL _u_initial) {
  PetscFunctionBeginUser;
  PetscScalar ***aU;
  PetscScalar **cX, **cY;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt iprev, icenter, ip, iux, iuy;
  
  PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
  PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));

  PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_LEFT, &iprev));
  
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy)); 
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));
  
  for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx + nEx[0]; ++ex) {
      aU[ey][ex][iuy] = _u_initial(cX[ex][icenter], cY[ey][iprev]);
      aU[ey][ex][iux] = _u_initial(cX[ex][iprev], cY[ey][icenter]);
      if (ey < starty + ny && ex < startx + nx) {
        aU[ey][ex][ip] = _u_initial(cX[ex][icenter], cY[ey][icenter]);
      }
    }
  }
  
  PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
  PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(DMLocalToGlobal(dm, uLocal, INSERT_VALUES, u));
  PetscFunctionReturn(0);
}

// ============================================================================
// Compute residual norm
// ============================================================================
PetscErrorCode ComputeResidualNorm(const DM &dm, 
                                   const Vec &u, const Vec &uLocal, 
                                   const Vec &uOld, const Vec &uOldLocal, 
                                   const Vec &f, const Vec &fLocal, 
                                   PetscInt Nx, PetscInt Ny,
                                   PetscReal alpha, PetscReal dt,
                                   PetscReal *residualNorm,
                                   PetscReal t, FUNC::EXACT _u_exact) {
  PetscFunctionBeginUser;
  PetscScalar ***aU, ***aUold, ***aF;
  PetscScalar **cX, **cY;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt iprev, icenter, iux, iuy;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  const PetscReal ix2 = 1.0 / (hx * hx);
  const PetscReal iy2 = 1.0 / (hy * hy);
  const PetscReal coef = alpha * dt;
  const PetscReal diag = 1.0 + 2.0 * coef * (ix2 + iy2);

  PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMGlobalToLocalBegin(dm, uOld, INSERT_VALUES, uOldLocal));
  PetscCall(DMGlobalToLocalEnd(dm, uOld, INSERT_VALUES, uOldLocal));
  PetscCall(DMGlobalToLocalBegin(dm, f, INSERT_VALUES, fLocal));
  PetscCall(DMGlobalToLocalEnd(dm, f, INSERT_VALUES, fLocal));
  PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
  PetscCall(DMStagVecGetArray(dm, uOldLocal, &aUold));
  PetscCall(DMStagVecGetArray(dm, fLocal, &aF));
  PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));
  PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_LEFT, &iprev));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));
  
  PetscReal localSum = 0.0;
  
  // Loop over DOWN edges (horizontal edges at y-boundaries between cells)
  for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
    // Skip boundary edges where we enforce BC
    if (ey == 0 || ey == Ny) continue;
    
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      const PetscScalar uc = aU[ey][ex][iuy];
      
      // Get neighbor values (with boundary handling)
      PetscScalar ul, ur, ud, uu;
      
      // Left neighbor
      if (ex == 0) {
        const PetscScalar x = 0.0;  // Left boundary x coordinate
        const PetscScalar y = cY[ey][iprev];
        ul = 2.0 * _u_exact(x, y, t, alpha) - uc;
      } else {
        ul = aU[ey][ex-1][iuy];
      }
      
      // Right neighbor
      if (ex == Nx - 1) {
        const PetscScalar x = 1.0;  // Right boundary x coordinate
        const PetscScalar y = cY[ey][iprev];
        ur = 2.0 * _u_exact(x, y, t, alpha) - uc;
      } else {
        ur = aU[ey][ex+1][iuy];
      }
      
      // Down neighbor (ey-1)
      if (ey == 1) {
        ud = _u_exact(cX[ex][icenter], cY[0][iprev], t, alpha);
      } else {
        ud = aU[ey-1][ex][iuy];
      }
      
      // Up neighbor (ey+1)
      if (ey == Ny - 1) {
        uu = _u_exact(cX[ex][icenter], cY[Ny][iprev], t, alpha);
      } else {
        uu = aU[ey+1][ex][iuy];
      }

      const PetscScalar lhs = diag * uc - coef * (ix2 * (ul + ur) + iy2 * (ud + uu));
      const PetscScalar rhs = aUold[ey][ex][iuy] + dt * aF[ey][ex][iuy];
      const PetscScalar r = rhs - lhs;
      localSum += PetscRealPart(r * r);
    }
  }
  
  // Loop over LEFT edges (vertical edges at x-boundaries between cells)
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx + nEx[0]; ++ex) {
      // Skip boundary edges where we enforce BC
      if (ex == 0 || ex == Nx) continue;
      
      const PetscScalar uc = aU[ey][ex][iux];
      
      // Get neighbor values (with boundary handling)
      PetscScalar ul, ur, ud, uu;
      
      // Left neighbor (ex-1)
      if (ex == 1) {
        ul = _u_exact(cX[0][iprev], cY[ey][icenter], t, alpha);
      } else {
        ul = aU[ey][ex-1][iux];
      }
      
      // Right neighbor (ex+1)
      if (ex == Nx - 1) {
        ur = _u_exact(cX[Nx][iprev], cY[ey][icenter], t, alpha);
      } else {
        ur = aU[ey][ex+1][iux];
      }
      
      // Down neighbor
      if (ey == 0) {
        const PetscScalar x = cX[ex][iprev];
        const PetscScalar y = 0.0;  // Bottom boundary y coordinate
        ud = 2.0 * _u_exact(x, y, t, alpha) - uc;
      } else {
        ud = aU[ey-1][ex][iux];
      }
      
      // Up neighbor
      if (ey == Ny - 1) {
        const PetscScalar x = cX[ex][iprev];
        const PetscScalar y = 1.0;  // Top boundary y coordinate
        uu = 2.0 * _u_exact(x, y, t, alpha) - uc;
      } else {
        uu = aU[ey+1][ex][iux];
      }

      const PetscScalar lhs = diag * uc - coef * (ix2 * (ul + ur) + iy2 * (ud + uu));
      const PetscScalar rhs = aUold[ey][ex][iux] + dt * aF[ey][ex][iux];
      const PetscScalar r = rhs - lhs;
      localSum += PetscRealPart(r * r);
    }
  }
  PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
  PetscCall(DMStagVecRestoreArray(dm, uOldLocal, &aUold));
  PetscCall(DMStagVecRestoreArray(dm, fLocal, &aF));
  PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscReal outNorm;
  MPI_Allreduce(&localSum, &outNorm, 1, MPIU_REAL, MPIU_SUM, PETSC_COMM_WORLD);
  *residualNorm = std::sqrt(outNorm);
  PetscFunctionReturn(0);
}

// ============================================================================
// Perform one red-black Gauss-Seidel sweep
// ============================================================================
PetscErrorCode GaussSeidelSweep(const DM &dm, Vec &u, Vec &uLocal, 
                                const Vec &uOld, const Vec &uOldLocal,
                                const Vec &f, const Vec &fLocal, 
                                PetscInt Nx, PetscInt Ny,
                                PetscReal alpha, PetscReal dt,
                                PetscReal t, FUNC::EXACT _u_exact) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  const PetscReal ix2 = 1.0 / (hx * hx);
  const PetscReal iy2 = 1.0 / (hy * hy);
  const PetscReal coef = alpha * dt;
  const PetscReal diag = 1.0 + 2.0 * coef * (ix2 + iy2);
  
  for (int color = 0; color < 2; ++color) {
    PetscScalar ***aU, ***aUold, ***aF;
    PetscScalar **cX, **cY;
    PetscInt startx, starty, nx, ny, nEx[2];
    PetscInt iprev, icenter, iux, iuy;
    PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
    PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
    PetscCall(DMGlobalToLocalBegin(dm, uOld, INSERT_VALUES, uOldLocal));
    PetscCall(DMGlobalToLocalEnd(dm, uOld, INSERT_VALUES, uOldLocal));
    PetscCall(DMGlobalToLocalBegin(dm, f, INSERT_VALUES, fLocal));
    PetscCall(DMGlobalToLocalEnd(dm, f, INSERT_VALUES, fLocal));
    PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
    PetscCall(DMStagVecGetArray(dm, uOldLocal, &aUold));
    PetscCall(DMStagVecGetArray(dm, fLocal, &aF));
    PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));
    PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
    PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
    PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_LEFT, &iprev));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));
    
    // Update DOWN edges
    for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
      for (PetscInt ex = startx; ex < startx + nx; ++ex) {
        // Apply Dirichlet BC at boundaries
        if (ey == 0 || ey == Ny) {
          const PetscScalar x = cX[ex][icenter];
          const PetscScalar y = cY[ey][iprev];
          aU[ey][ex][iuy] = _u_exact(x, y, t, alpha);
          continue;
        }
        
        if (((ex + ey) & 1) != color)
          continue;

        const PetscScalar uc = aU[ey][ex][iuy];
        
        // Get neighbor values
        PetscScalar ul, ur, ud, uu;
        
        if (ex == 0) {
          // Left boundary: extrapolate from center to ghost
          const PetscScalar x = 0.0;  // Left boundary x coordinate
          const PetscScalar y = cY[ey][iprev];
          ul = 2.0 * _u_exact(x, y, t, alpha) - uc;
        } else {
          ul = aU[ey][ex-1][iuy];
        }
        
        if (ex == Nx - 1) {
          // Right boundary: extrapolate from center to ghost
          const PetscScalar x = 1.0;  // Right boundary x coordinate
          const PetscScalar y = cY[ey][iprev];
          ur = 2.0 * _u_exact(x, y, t, alpha) - uc;
        } else {
          ur = aU[ey][ex+1][iuy];
        }
        
        // Down neighbor (ey-1)
        if (ey == 1) {
          ud = _u_exact(cX[ex][icenter], cY[0][iprev], t, alpha);
        } else {
          ud = aU[ey-1][ex][iuy];
        }
        
        // Up neighbor (ey+1)
        if (ey == Ny - 1) {
          uu = _u_exact(cX[ex][icenter], cY[Ny][iprev], t, alpha);
        } else {
          uu = aU[ey+1][ex][iuy];
        }

        const PetscScalar rhs = aUold[ey][ex][iuy] + dt * aF[ey][ex][iuy];
        const PetscScalar unew = (rhs + coef * (ix2 * (ul + ur) + iy2 * (ud + uu))) / diag;
        aU[ey][ex][iuy] = unew;
      }
    }
    
    // Update LEFT edges
    for (PetscInt ey = starty; ey < starty + ny; ++ey) {
      for (PetscInt ex = startx; ex < startx + nx + nEx[0]; ++ex) {
        // Apply Dirichlet BC at boundaries
        if (ex == 0 || ex == Nx) {
          const PetscScalar x = cX[ex][iprev];
          const PetscScalar y = cY[ey][icenter];
          aU[ey][ex][iux] = _u_exact(x, y, t, alpha);
          continue;
        }
        
        if (((ex + ey) & 1) != color)
          continue;

        const PetscScalar uc = aU[ey][ex][iux];
        
        // Get neighbor values
        PetscScalar ul = (ex == 1) ? _u_exact(cX[0][iprev], cY[ey][icenter], t, alpha) : aU[ey][ex-1][iux];
        PetscScalar ur = (ex == Nx - 1) ? _u_exact(cX[Nx][iprev], cY[ey][icenter], t, alpha) : aU[ey][ex+1][iux];
        PetscScalar ud = (ey == 0) ? 2.0 * _u_exact(cX[ex][iprev], 0.0, t, alpha) - uc : aU[ey-1][ex][iux];
        PetscScalar uu = (ey == Ny - 1) ? 2.0 * _u_exact(cX[ex][iprev], 1.0, t, alpha) - uc : aU[ey+1][ex][iux];

        const PetscScalar rhs = aUold[ey][ex][iux] + dt * aF[ey][ex][iux];
        const PetscScalar unew = (rhs + coef * (ix2 * (ul + ur) + iy2 * (ud + uu))) / diag;
        aU[ey][ex][iux] = unew;
      }
    }

    PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
    PetscCall(DMStagVecRestoreArray(dm, uOldLocal, &aUold));
    PetscCall(DMStagVecRestoreArray(dm, fLocal, &aF));
    PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
    PetscCall(DMLocalToGlobal(dm, uLocal, INSERT_VALUES, u));
  }
  PetscFunctionReturn(0);
}

// ============================================================================
// Compute L2 error against exact solution at time t
// ============================================================================
PetscErrorCode ComputeL2Error(const DM &dm, const Vec &u, Vec &uLocal,
                              PetscInt Nx, PetscInt Ny, 
                              PetscReal t, PetscReal alpha,
                              PetscReal *l2Error, FUNC::EXACT _u_exact) {
  PetscFunctionBeginUser;
  PetscScalar ***aU;
  PetscScalar **cX, **cY;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt iprev, icenter, iux, iuy;
  
  PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
  PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));
  PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_LEFT, &iprev));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));

  // Manually compute L2 error for both DOWN and LEFT edge DOFs
  PetscReal localError2 = 0.0;
  
  // DOWN edges (including boundary edges)
  for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      const PetscScalar x = cX[ex][icenter];
      const PetscScalar y = cY[ey][iprev];
      const PetscScalar uExact = _u_exact(x, y, t, alpha);
      const PetscScalar diff = aU[ey][ex][iuy] - uExact;
      localError2 += diff * diff;
    }
  }
  
  // LEFT edges (including boundary edges)
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx + nEx[0]; ++ex) {
      const PetscScalar x = cX[ex][iprev];
      const PetscScalar y = cY[ey][icenter];
      const PetscScalar uExact = _u_exact(x, y, t, alpha);
      const PetscScalar diff = aU[ey][ex][iux] - uExact;
      localError2 += diff * diff;
    }
  }
  
  PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
  PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  
  // Sum across all processes
  PetscReal globalError2;
  PetscCall(MPI_Allreduce(&localError2, &globalError2, 1, MPIU_REAL, MPI_SUM, PetscObjectComm((PetscObject)dm)));
  
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  *l2Error = std::sqrt(globalError2 * hx * hy);
  PetscFunctionReturn(0);
}

int main(int argc, char **argv) {
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, "Heat equation solver using DMStag"));

  DM dm;
  Vec u, uLocal, uOld, uOldLocal, f, fLocal;
  Vec div, divLocal;  // For divergence computation
  PetscInt Nx = 64, Ny = 64;
  const PetscInt dof0 = 0, dof1 = 1, dof2 = 1;  // Added dof2=1 for divergence at cell center
  const PetscInt stencilWidth = 1;
  const PetscReal alpha = 0.1;  // Thermal diffusivity
  PetscReal dt = 0.001;         // Time step
  PetscReal T_final = 0.1;      // Final time
  const PetscReal tol = 1e-8;   // Tolerance for GS iterations
  const PetscInt maxIts = 10000;
  int rank;
  PetscBool compute_error = PETSC_FALSE;
  PetscBool output_divergence = PETSC_FALSE;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  // Get options from command line
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-nx", &Nx, NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-ny", &Ny, NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-dt", &dt, NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-T", &T_final, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-heat_check_error", &compute_error, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-output_divergence", &output_divergence, NULL));

  // Create DMStag
  PetscCall(DMStagCreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                           Nx, Ny, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2,
                           DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, &dm));
  PetscCall(DMSetFromOptions(dm));
  PetscCall(DMSetUp(dm));
  PetscCall(DMStagSetUniformCoordinatesProduct(dm, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0));
  
  // Create vectors
  PetscCall(DMCreateGlobalVector(dm, &u));
  PetscCall(DMCreateGlobalVector(dm, &uOld));
  PetscCall(DMCreateGlobalVector(dm, &f));
  PetscCall(DMCreateGlobalVector(dm, &div));
  PetscCall(DMGetLocalVector(dm, &uLocal));
  PetscCall(DMGetLocalVector(dm, &uOldLocal));
  PetscCall(DMGetLocalVector(dm, &fLocal));
  PetscCall(DMGetLocalVector(dm, &divLocal));
  
  // Get global sizes
  PetscInt Nglob[2];
  PetscCall(DMStagGetGlobalSizes(dm, &Nglob[0], &Nglob[1], NULL));
  
  if (rank == 0) {
    std::cout << "=== Heat Equation Solver ===" << std::endl;
    std::cout << "Grid: " << Nglob[0] << " x " << Nglob[1] << std::endl;
    std::cout << "Alpha: " << alpha << ", dt: " << dt << ", T_final: " << T_final << std::endl;
    std::cout << "Time steps: " << (int)(T_final / dt) << std::endl;
  }

  // Setup initial condition
  PetscCall(SetupInitialCondition(dm, u, uLocal, u_initial));
  
  // Time stepping loop
  PetscReal t = 0.0;
  PetscInt step = 0;
  const PetscInt output_interval = std::max(1, (int)(T_final / dt / 10));
  
  while (t < T_final) {
    step++;
    t += dt;
    
    // u_old = u^n
    PetscCall(VecCopy(u, uOld));
    
    // Setup RHS (source term at time t)
    PetscCall(SetupRHS(dm, f, fLocal, t, rhs_f));
    
    // Solve (I - dt*alpha*Laplace)*u^{n+1} = u^n + dt*f using Gauss-Seidel
    PetscReal resNorm = 0.0;
    PetscInt its = 0;
    
    for (its = 1; its <= maxIts; ++its) {
      PetscCall(GaussSeidelSweep(dm, u, uLocal, uOld, uOldLocal, f, fLocal, 
                                 Nglob[0], Nglob[1], alpha, dt, t, u_exact));
      PetscCall(ComputeResidualNorm(dm, u, uLocal, uOld, uOldLocal, f, fLocal, 
                                    Nglob[0], Nglob[1], alpha, dt, &resNorm, t, u_exact));
      if (resNorm <= tol)
        break;
    }
    
    if (rank == 0 && (step % output_interval == 0 || step == 1)) {
      std::cout << "Step " << step << ", t=" << t << ", GS its=" << its 
                << ", ||r||=" << resNorm;
      
      // Compute and output divergence statistics if requested
      if (output_divergence) {
        PetscCall(ComputeDivergence(dm, u, uLocal, div, divLocal, Nglob[0], Nglob[1]));
        PetscReal divNorm, divMax, divMin;
        PetscCall(VecNorm(div, NORM_2, &divNorm));
        PetscCall(VecMax(div, NULL, &divMax));
        PetscCall(VecMin(div, NULL, &divMin));
        
        // Normalize L2 norm by grid size for comparison
        const PetscReal hx = 1.0 / Nglob[0];
        const PetscReal hy = 1.0 / Nglob[1];
        const PetscReal normalizedDivNorm = divNorm * std::sqrt(hx * hy);
        
        std::cout << ", ||div||=" << normalizedDivNorm 
                  << ", div_max=" << divMax 
                  << ", div_min=" << divMin;
      }
      std::cout << std::endl;
    }
  }
  
  if (rank == 0) {
    std::cout << "\nTime integration complete. Final time: " << t << std::endl;
  }
  
  // Compute error if requested
  if (compute_error) {
    PetscReal l2Error = 0.0;
    PetscCall(ComputeL2Error(dm, u, uLocal, Nglob[0], Nglob[1], t, alpha, &l2Error, u_exact));
    if (rank == 0) {
      std::cout << "||u(T) - u_exact(T)||_L2 = " << l2Error << std::endl;
    }
  }

  // Cleanup
  PetscCall(DMRestoreLocalVector(dm, &uLocal));
  PetscCall(DMRestoreLocalVector(dm, &uOldLocal));
  PetscCall(DMRestoreLocalVector(dm, &fLocal));
  PetscCall(DMRestoreLocalVector(dm, &divLocal));
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&uOld));
  PetscCall(VecDestroy(&f));
  PetscCall(VecDestroy(&div));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFinalize());
  return 0;
}
