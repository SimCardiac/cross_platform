#include <petscdm.h>
#include <petscdmstag.h>
#include <petscksp.h>
#include <petscsys.h>

#include <cmath>
#include <iostream>

#include "manufactured_solutions.h"

// Use manufactured solution from library
using namespace HEAT::DECAY_2D;
namespace FUNC = HEAT::FUNC_2D;

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
                                   PetscReal *residualNorm) {
  PetscFunctionBeginUser;
  PetscScalar ***aU, ***aUold, ***aF;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt iux, iuy;
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
        ul = -uc;  // Dirichlet BC: u=0 at boundary
      } else {
        ul = aU[ey][ex-1][iuy];
      }
      
      // Right neighbor
      if (ex == Nx - 1) {
        ur = -uc;
      } else {
        ur = aU[ey][ex+1][iuy];
      }
      
      // Down neighbor (ey-1)
      if (ey == 1) {
        ud = 0.0;  // Boundary edge at ey=0 has u=0
      } else {
        ud = aU[ey-1][ex][iuy];
      }
      
      // Up neighbor (ey+1)
      if (ey == Ny - 1) {
        uu = 0.0;  // Boundary edge at ey=Ny has u=0
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
        ul = 0.0;  // Boundary edge at ex=0 has u=0
      } else {
        ul = aU[ey][ex-1][iux];
      }
      
      // Right neighbor (ex+1)
      if (ex == Nx - 1) {
        ur = 0.0;  // Boundary edge at ex=Nx has u=0
      } else {
        ur = aU[ey][ex+1][iux];
      }
      
      // Down neighbor
      if (ey == 0) {
        ud = -uc;  // Dirichlet BC: u=0 at boundary
      } else {
        ud = aU[ey-1][ex][iux];
      }
      
      // Up neighbor
      if (ey == Ny - 1) {
        uu = -uc;
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
                                PetscReal alpha, PetscReal dt) {
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
    PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));
    
    // Update DOWN edges
    for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
      // Skip boundary edges where we enforce BC
      if (ey == 0 || ey == Ny) continue;
      
      for (PetscInt ex = startx; ex < startx + nx; ++ex) {
        if (((ex + ey) & 1) != color)
          continue;

        const PetscScalar uc = aU[ey][ex][iuy];
        
        // Get neighbor values
        PetscScalar ul, ur, ud, uu;
        
        if (ex == 0) {
          ul = -uc;
        } else {
          ul = aU[ey][ex-1][iuy];
        }
        
        if (ex == Nx - 1) {
          ur = -uc;
        } else {
          ur = aU[ey][ex+1][iuy];
        }
        
        // Down neighbor (ey-1)
        if (ey == 1) {
          ud = 0.0;  // Boundary edge at ey=0 has u=0
        } else {
          ud = aU[ey-1][ex][iuy];
        }
        
        // Up neighbor (ey+1)
        if (ey == Ny - 1) {
          uu = 0.0;  // Boundary edge at ey=Ny has u=0
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
        // Skip boundary edges where we enforce BC
        if (ex == 0 || ex == Nx) continue;
        
        if (((ex + ey) & 1) != color)
          continue;

        const PetscScalar uc = aU[ey][ex][iux];
        
        // Get neighbor values
        PetscScalar ul, ur, ud, uu;
        
        // Left neighbor (ex-1)
        if (ex == 1) {
          ul = 0.0;
        } else {
          ul = aU[ey][ex-1][iux];
        }
        
        // Right neighbor (ex+1)
        if (ex == Nx - 1) {
          ur = 0.0;
        } else {
          ur = aU[ey][ex+1][iux];
        }
        
        // Down neighbor
        if (ey == 0) {
          ud = -uc;
        } else {
          ud = aU[ey-1][ex][iux];
        }
        
        // Up neighbor
        if (ey == Ny - 1) {
          uu = -uc;
        } else {
          uu = aU[ey+1][ex][iux];
        }

        const PetscScalar rhs = aUold[ey][ex][iux] + dt * aF[ey][ex][iux];
        const PetscScalar unew = (rhs + coef * (ix2 * (ul + ur) + iy2 * (ud + uu))) / diag;
        aU[ey][ex][iux] = unew;
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
  
  // DOWN edges (interior edges)
  for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
    // Skip boundary edges at ey=0 and ey=Ny (where we enforce BC)
    if (ey == 0 || ey == Ny) continue;
    
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      const PetscScalar x = cX[ex][icenter];
      const PetscScalar y = cY[ey][iprev];
      const PetscScalar uExact = _u_exact(x, y, t, alpha);
      const PetscScalar diff = aU[ey][ex][iuy] - uExact;
      localError2 += diff * diff;
    }
  }
  
  // LEFT edges (interior edges)
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx + nEx[0]; ++ex) {
      // Skip boundary edges at ex=0 and ex=Nx
      if (ex == 0 || ex == Nx) continue;
      
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
  PetscInt Nx = 64, Ny = 64;
  const PetscInt dof0 = 0, dof1 = 1, dof2 = 0; 
  const PetscInt stencilWidth = 1;
  const PetscReal alpha = 0.1;  // Thermal diffusivity
  PetscReal dt = 0.001;         // Time step
  PetscReal T_final = 0.1;      // Final time
  const PetscReal tol = 1e-8;   // Tolerance for GS iterations
  const PetscInt maxIts = 10000;
  int rank;
  PetscBool compute_error = PETSC_FALSE;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  // Get options from command line
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-nx", &Nx, NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-ny", &Ny, NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-dt", &dt, NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-T", &T_final, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-heat_check_error", &compute_error, NULL));

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
  PetscCall(DMGetLocalVector(dm, &uLocal));
  PetscCall(DMGetLocalVector(dm, &uOldLocal));
  PetscCall(DMGetLocalVector(dm, &fLocal));
  
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
                                 Nglob[0], Nglob[1], alpha, dt));
      PetscCall(ComputeResidualNorm(dm, u, uLocal, uOld, uOldLocal, f, fLocal, 
                                    Nglob[0], Nglob[1], alpha, dt, &resNorm));
      if (resNorm <= tol)
        break;
    }
    
    if (rank == 0 && (step % output_interval == 0 || step == 1)) {
      std::cout << "Step " << step << ", t=" << t << ", GS its=" << its 
                << ", ||r||=" << resNorm << std::endl;
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
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&uOld));
  PetscCall(VecDestroy(&f));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFinalize());
  return 0;
}
