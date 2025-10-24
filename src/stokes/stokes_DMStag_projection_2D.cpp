#include <petscdm.h>
#include <petscdmstag.h>
#include <petscksp.h>
#include <petscsys.h>
#include <cmath>
#include <iostream>

#include "manufactured_solutions.h"

// Use Stokes manufactured solution from library
using namespace STOKES::TAYLOR_GREEN_STEADY_2D;

// ============================================================================
// Setup force vector
// ============================================================================
PetscErrorCode SetupForce(const DM &dm, Vec &f, Vec &fLocal) {
  PetscFunctionBeginUser;
  PetscScalar ***aF;
  PetscScalar **cX, **cY;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt iprev, icenter, iux, iuy;
  
  PetscCall(DMGlobalToLocalBegin(dm, f, INSERT_VALUES, fLocal));
  PetscCall(DMGlobalToLocalEnd(dm, f, INSERT_VALUES, fLocal));
  PetscCall(DMStagVecGetArray(dm, fLocal, &aF));
  PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));
  PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_LEFT, &iprev));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy));

  // Force at LEFT edges (u component)
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx + nEx[0]; ++ex) {
      const PetscScalar x = cX[ex][iprev];
      const PetscScalar y = cY[ey][icenter];
      aF[ey][ex][iux] = fx(x, y);
    }
  }
  
  // Force at DOWN edges (v component)
  for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      const PetscScalar x = cX[ex][icenter];
      const PetscScalar y = cY[ey][iprev];
      aF[ey][ex][iuy] = fy(x, y);
    }
  }
  
  PetscCall(DMStagVecRestoreArray(dm, fLocal, &aF));
  PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(DMLocalToGlobal(dm, fLocal, INSERT_VALUES, f));
  PetscFunctionReturn(0);
}

// ============================================================================
// Step 1: Solve momentum equation for intermediate velocity (without pressure)
// -ν∇²u* = f
// ============================================================================
PetscErrorCode SolveViscousStep(const DM &dm, Vec &ustar, Vec &ustarLocal, 
                                 const Vec &f, const Vec &fLocal,
                                 PetscInt Nx, PetscInt Ny, PetscReal nu,
                                 PetscReal tol, PetscInt maxIts) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  const PetscReal ix2 = 1.0 / (hx * hx);
  const PetscReal iy2 = 1.0 / (hy * hy);
  const PetscReal diag = 2.0 * nu * (ix2 + iy2);
  
  PetscInt iux, iuy;
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy));
  
  // Gauss-Seidel iterations for viscous solve
  for (PetscInt its = 0; its < maxIts; ++its) {
    for (int color = 0; color < 2; ++color) {
      PetscScalar ***aU, ***aF;
      PetscInt startx, starty, nx, ny, nEx[2];
      
      PetscCall(DMGlobalToLocalBegin(dm, ustar, INSERT_VALUES, ustarLocal));
      PetscCall(DMGlobalToLocalEnd(dm, ustar, INSERT_VALUES, ustarLocal));
      PetscCall(DMGlobalToLocalBegin(dm, f, INSERT_VALUES, fLocal));
      PetscCall(DMGlobalToLocalEnd(dm, f, INSERT_VALUES, fLocal));
      PetscCall(DMStagVecGetArray(dm, ustarLocal, &aU));
      PetscCall(DMStagVecGetArray(dm, fLocal, &aF));
      PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));
      
      // Update u component (LEFT edges)
      for (PetscInt ey = starty; ey < starty + ny; ++ey) {
        for (PetscInt ex = startx; ex < startx + nx + nEx[0]; ++ex) {
          if (((ex + ey) & 1) != color) continue;
          
          // Boundary conditions: u = 0 at x=0 and x=1
          if (ex == 0 || ex == Nx) {
            aU[ey][ex][iux] = 0.0;
            continue;
          }
          
          const PetscScalar uc = aU[ey][ex][iux];
          PetscScalar ul = aU[ey][ex-1][iux];
          PetscScalar ur = aU[ey][ex+1][iux];
          PetscScalar ud = (ey == 0) ? -uc : aU[ey-1][ex][iux];
          PetscScalar uu = (ey == Ny-1) ? -uc : aU[ey+1][ex][iux];
          
          const PetscScalar rhs = aF[ey][ex][iux];
          aU[ey][ex][iux] = (rhs + nu * (ix2 * (ul + ur) + iy2 * (ud + uu))) / diag;
        }
      }
      
      // Update v component (DOWN edges)
      for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
        for (PetscInt ex = startx; ex < startx + nx; ++ex) {
          if (((ex + ey) & 1) != color) continue;
          
          // Boundary conditions: v = 0 at y=0 and y=1
          if (ey == 0 || ey == Ny) {
            aU[ey][ex][iuy] = 0.0;
            continue;
          }
          
          const PetscScalar vc = aU[ey][ex][iuy];
          PetscScalar vl = (ex == 0) ? -vc : aU[ey][ex-1][iuy];
          PetscScalar vr = (ex == Nx-1) ? -vc : aU[ey][ex+1][iuy];
          PetscScalar vd = aU[ey-1][ex][iuy];
          PetscScalar vu = aU[ey+1][ex][iuy];
          
          const PetscScalar rhs = aF[ey][ex][iuy];
          aU[ey][ex][iuy] = (rhs + nu * (ix2 * (vl + vr) + iy2 * (vd + vu))) / diag;
        }
      }
      
      PetscCall(DMStagVecRestoreArray(dm, ustarLocal, &aU));
      PetscCall(DMStagVecRestoreArray(dm, fLocal, &aF));
      PetscCall(DMLocalToGlobal(dm, ustarLocal, INSERT_VALUES, ustar));
    }
  }
  
  PetscFunctionReturn(0);
}

// ============================================================================
// Step 2: Solve pressure Poisson equation
// ∇²p = ∇·u* 
// ============================================================================
PetscErrorCode SolvePressurePoisson(const DM &dm, Vec &p, Vec &pLocal,
                                     const Vec &ustar, const Vec &ustarLocal,
                                     PetscInt Nx, PetscInt Ny,
                                     PetscReal tol, PetscInt maxIts) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  const PetscReal ix = 1.0 / hx;
  const PetscReal iy = 1.0 / hy;
  const PetscReal ix2 = 1.0 / (hx * hx);
  const PetscReal iy2 = 1.0 / (hy * hy);
  const PetscReal diag = 2.0 * (ix2 + iy2);
  
  PetscInt iux, iuy, ip;
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));
  
  // Gauss-Seidel iterations for pressure Poisson
  for (PetscInt its = 0; its < maxIts; ++its) {
    for (int color = 0; color < 2; ++color) {
      PetscScalar ***aP, ***aU;
      PetscInt startx, starty, nx, ny, nEx[2];
      
      PetscCall(DMGlobalToLocalBegin(dm, p, INSERT_VALUES, pLocal));
      PetscCall(DMGlobalToLocalEnd(dm, p, INSERT_VALUES, pLocal));
      PetscCall(DMGlobalToLocalBegin(dm, ustar, INSERT_VALUES, ustarLocal));
      PetscCall(DMGlobalToLocalEnd(dm, ustar, INSERT_VALUES, ustarLocal));
      PetscCall(DMStagVecGetArray(dm, pLocal, &aP));
      PetscCall(DMStagVecGetArray(dm, ustarLocal, &aU));
      PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));
      
      for (PetscInt ey = starty; ey < starty + ny; ++ey) {
        for (PetscInt ex = startx; ex < startx + nx; ++ex) {
          if (((ex + ey) & 1) != color) continue;
          
          const PetscScalar pc = aP[ey][ex][ip];
          PetscScalar pl = (ex == 0) ? -pc : aP[ey][ex-1][ip];
          PetscScalar pr = (ex == Nx-1) ? -pc : aP[ey][ex+1][ip];
          PetscScalar pd = (ey == 0) ? -pc : aP[ey-1][ex][ip];
          PetscScalar pu = (ey == Ny-1) ? -pc : aP[ey+1][ex][ip];
          
          // RHS = ∇·u* = ∂u/∂x + ∂v/∂y
          const PetscScalar dudx = ix * (aU[ey][ex+1][iux] - aU[ey][ex][iux]);
          const PetscScalar dvdy = iy * (aU[ey+1][ex][iuy] - aU[ey][ex][iuy]);
          const PetscScalar div_ustar = dudx + dvdy;
          
          aP[ey][ex][ip] = (div_ustar + ix2 * (pl + pr) + iy2 * (pd + pu)) / diag;
        }
      }
      
      PetscCall(DMStagVecRestoreArray(dm, pLocal, &aP));
      PetscCall(DMStagVecRestoreArray(dm, ustarLocal, &aU));
      PetscCall(DMLocalToGlobal(dm, pLocal, INSERT_VALUES, p));
    }
  }
  
  PetscFunctionReturn(0);
}

// ============================================================================
// Step 3: Velocity correction
// u^{n+1} = u* - ∇p
// ============================================================================
PetscErrorCode CorrectVelocity(const DM &dm, Vec &u, Vec &uLocal,
                                const Vec &ustar, const Vec &ustarLocal,
                                const Vec &p, const Vec &pLocal,
                                PetscInt Nx, PetscInt Ny) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  const PetscReal ix = 1.0 / hx;
  const PetscReal iy = 1.0 / hy;
  
  PetscScalar ***aU, ***aUstar, ***aP;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt iux, iuy, ip;
  
  PetscCall(DMGlobalToLocalBegin(dm, ustar, INSERT_VALUES, ustarLocal));
  PetscCall(DMGlobalToLocalEnd(dm, ustar, INSERT_VALUES, ustarLocal));
  PetscCall(DMGlobalToLocalBegin(dm, p, INSERT_VALUES, pLocal));
  PetscCall(DMGlobalToLocalEnd(dm, p, INSERT_VALUES, pLocal));
  PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
  
  PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
  PetscCall(DMStagVecGetArray(dm, ustarLocal, &aUstar));
  PetscCall(DMStagVecGetArray(dm, pLocal, &aP));
  PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));
  
  // Correct u component: u = u* - ∂p/∂x
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx + nEx[0]; ++ex) {
      if (ex == 0 || ex == Nx) continue;
      const PetscScalar dpdx = ix * (aP[ey][ex][ip] - aP[ey][ex-1][ip]);
      aU[ey][ex][iux] = aUstar[ey][ex][iux] - dpdx;
    }
  }
  
  // Correct v component: v = v* - ∂p/∂y
  for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      if (ey == 0 || ey == Ny) continue;
      const PetscScalar dpdy = iy * (aP[ey][ex][ip] - aP[ey-1][ex][ip]);
      aU[ey][ex][iuy] = aUstar[ey][ex][iuy] - dpdy;
    }
  }
  
  PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
  PetscCall(DMStagVecRestoreArray(dm, ustarLocal, &aUstar));
  PetscCall(DMStagVecRestoreArray(dm, pLocal, &aP));
  PetscCall(DMLocalToGlobal(dm, uLocal, INSERT_VALUES, u));
  PetscFunctionReturn(0);
}

// ============================================================================
// Compute L2 errors
// ============================================================================
PetscErrorCode ComputeL2Error(const DM &dm, const Vec &u, Vec &uLocal,
                               const Vec &p, Vec &pLocal,
                               PetscInt Nx, PetscInt Ny,
                               PetscReal *u_error, PetscReal *v_error, PetscReal *p_error) {
  PetscFunctionBeginUser;
  PetscScalar ***aU, ***aP;
  PetscScalar **cX, **cY;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt iprev, icenter, iux, iuy, ip;
  
  PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
  PetscCall(DMGlobalToLocalBegin(dm, p, INSERT_VALUES, pLocal));
  PetscCall(DMGlobalToLocalEnd(dm, p, INSERT_VALUES, pLocal));
  
  PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
  PetscCall(DMStagVecGetArray(dm, pLocal, &aP));
  PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));
  PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_LEFT, &iprev));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));
  
  PetscReal localError_u = 0.0, localError_v = 0.0, localError_p = 0.0;
  
  // u component error
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx + nEx[0]; ++ex) {
      const PetscScalar x = cX[ex][iprev];
      const PetscScalar y = cY[ey][icenter];
      const PetscScalar diff = aU[ey][ex][iux] - u_exact(x, y);
      localError_u += diff * diff;
    }
  }
  
  // v component error
  for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      const PetscScalar x = cX[ex][icenter];
      const PetscScalar y = cY[ey][iprev];
      const PetscScalar diff = aU[ey][ex][iuy] - v_exact(x, y);
      localError_v += diff * diff;
    }
  }
  
  // p component error
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      const PetscScalar x = cX[ex][icenter];
      const PetscScalar y = cY[ey][icenter];
      const PetscScalar diff = aP[ey][ex][ip] - p_exact(x, y);
      localError_p += diff * diff;
    }
  }
  
  PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
  PetscCall(DMStagVecRestoreArray(dm, pLocal, &aP));
  PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  
  PetscReal globalError_u, globalError_v, globalError_p;
  PetscCall(MPI_Allreduce(&localError_u, &globalError_u, 1, MPIU_REAL, MPI_SUM, PetscObjectComm((PetscObject)dm)));
  PetscCall(MPI_Allreduce(&localError_v, &globalError_v, 1, MPIU_REAL, MPI_SUM, PetscObjectComm((PetscObject)dm)));
  PetscCall(MPI_Allreduce(&localError_p, &globalError_p, 1, MPIU_REAL, MPI_SUM, PetscObjectComm((PetscObject)dm)));
  
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  *u_error = std::sqrt(globalError_u * hx * hy);
  *v_error = std::sqrt(globalError_v * hx * hy);
  *p_error = std::sqrt(globalError_p * hx * hy);
  
  PetscFunctionReturn(0);
}

int main(int argc, char **argv) {
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, "Stokes equation solver using projection method"));

  DM dm;
  Vec u, uLocal, ustar, ustarLocal, p, pLocal, f, fLocal;
  PetscInt Nx = 32, Ny = 32;
  const PetscInt dof0 = 0, dof1 = 1, dof2 = 1;  // edge velocities + element pressure
  const PetscInt stencilWidth = 1;
  const PetscReal nu = 1.0;  // Viscosity
  const PetscReal tol = 1e-8;
  const PetscInt maxIts = 10000;
  int rank;
  PetscBool compute_error = PETSC_FALSE;
  PetscBool convergence_test = PETSC_FALSE;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-nx", &Nx, NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-ny", &Ny, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-stokes_check_error", &compute_error, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-convergence_test", &convergence_test, NULL));

  // Create DMStag with staggered grid
  PetscCall(DMStagCreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                           Nx, Ny, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2,
                           DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, &dm));
  PetscCall(DMSetFromOptions(dm));
  PetscCall(DMSetUp(dm));
  PetscCall(DMStagSetUniformCoordinatesProduct(dm, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0));
  
  // Create vectors
  PetscCall(DMCreateGlobalVector(dm, &u));
  PetscCall(DMCreateGlobalVector(dm, &ustar));
  PetscCall(DMCreateGlobalVector(dm, &p));
  PetscCall(DMCreateGlobalVector(dm, &f));
  PetscCall(DMGetLocalVector(dm, &uLocal));
  PetscCall(DMGetLocalVector(dm, &ustarLocal));
  PetscCall(DMGetLocalVector(dm, &pLocal));
  PetscCall(DMGetLocalVector(dm, &fLocal));
  
  PetscCall(VecSet(u, 0.0));
  PetscCall(VecSet(ustar, 0.0));
  PetscCall(VecSet(p, 0.0));
  
  PetscInt Nglob[2];
  PetscCall(DMStagGetGlobalSizes(dm, &Nglob[0], &Nglob[1], NULL));
  
  if (rank == 0) {
    std::cout << "=== Stokes Equation Solver (Projection Method) ===" << std::endl;
    std::cout << "Grid: " << Nglob[0] << " x " << Nglob[1] << std::endl;
    std::cout << "Viscosity: " << nu << std::endl;
  }

  // Setup force
  PetscCall(SetupForce(dm, f, fLocal));
  
  // Projection method: 3 steps
  if (rank == 0) std::cout << "\nStep 1: Solving viscous momentum equation..." << std::endl;
  PetscCall(SolveViscousStep(dm, ustar, ustarLocal, f, fLocal, Nglob[0], Nglob[1], nu, tol, maxIts));
  
  if (rank == 0) std::cout << "Step 2: Solving pressure Poisson equation..." << std::endl;
  PetscCall(SolvePressurePoisson(dm, p, pLocal, ustar, ustarLocal, Nglob[0], Nglob[1], tol, maxIts));
  
  if (rank == 0) std::cout << "Step 3: Correcting velocity..." << std::endl;
  PetscCall(CorrectVelocity(dm, u, uLocal, ustar, ustarLocal, p, pLocal, Nglob[0], Nglob[1]));
  
  if (rank == 0) std::cout << "\nProjection method complete." << std::endl;
  
  // Compute errors
  if (compute_error || convergence_test) {
    PetscReal u_err, v_err, p_err;
    PetscCall(ComputeL2Error(dm, u, uLocal, p, pLocal, Nglob[0], Nglob[1], &u_err, &v_err, &p_err));
    
    if (rank == 0) {
      std::cout << "\n||u - u_exact||_L2 = " << u_err << std::endl;
      std::cout << "||v - v_exact||_L2 = " << v_err << std::endl;
      std::cout << "||p - p_exact||_L2 = " << p_err << std::endl;
      
      if (convergence_test) {
        const PetscReal hx = 1.0 / Nglob[0];
        const PetscReal hy = 1.0 / Nglob[1];
        const PetscReal h = std::sqrt(hx * hy);
        std::cout << "CONVERGENCE: " << Nglob[0] << " " << h << " " 
                  << u_err << " " << v_err << " " << p_err << std::endl;
      }
    }
  }

  // Cleanup
  PetscCall(DMRestoreLocalVector(dm, &uLocal));
  PetscCall(DMRestoreLocalVector(dm, &ustarLocal));
  PetscCall(DMRestoreLocalVector(dm, &pLocal));
  PetscCall(DMRestoreLocalVector(dm, &fLocal));
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&ustar));
  PetscCall(VecDestroy(&p));
  PetscCall(VecDestroy(&f));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFinalize());
  return 0;
}
