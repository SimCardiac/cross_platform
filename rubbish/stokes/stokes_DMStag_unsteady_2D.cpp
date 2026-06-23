#include <petscdm.h>
#include <petscdmstag.h>
#include <petscksp.h>
#include <petscsys.h>
#include <cmath>
#include <iostream>

#include "../analytical/unsteady.h"

// Use time-dependent Stokes manufactured solution from library
using namespace UNSTEADY::TAYLOR_GREEN_2D;

// ============================================================================
// Setup right-hand side vector (force terms)
// ============================================================================
PetscErrorCode SetupRHS(const DM &dm, Vec &f, Vec &fLocal, PetscScalar t) {
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

  // Force at DOWN edges (v component)
  for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      aF[ey][ex][iuy] = fy_stokes(cX[ex][icenter], cY[ey][iprev], t);
    }
  }
  
  // Force at LEFT edges (u component)
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx + nEx[0]; ++ex) {
      aF[ey][ex][iux] = fx_stokes(cX[ex][iprev], cY[ey][icenter], t);
    }
  }
  
  // Pressure DOF (element center) - set to zero (not used in RHS)
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      aF[ey][ex][ip] = 0.0;
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
PetscErrorCode SetupInitialCondition(const DM &dm, Vec &sol, Vec &solLocal) {
  PetscFunctionBeginUser;
  PetscScalar ***aSol;
  PetscScalar **cX, **cY;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt iprev, icenter, ip, iux, iuy;
  
  PetscCall(DMGlobalToLocalBegin(dm, sol, INSERT_VALUES, solLocal));
  PetscCall(DMGlobalToLocalEnd(dm, sol, INSERT_VALUES, solLocal));
  PetscCall(DMStagVecGetArray(dm, solLocal, &aSol));
  PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));

  PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_LEFT, &iprev));
  
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy)); 
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));
  
  // Initial velocity v at DOWN edges
  for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      aSol[ey][ex][iuy] = v_exact(cX[ex][icenter], cY[ey][iprev], 0.0);
    }
  }
  
  // Initial velocity u at LEFT edges
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx + nEx[0]; ++ex) {
      aSol[ey][ex][iux] = u_exact(cX[ex][iprev], cY[ey][icenter], 0.0);
    }
  }
  
  // Initial pressure at element centers
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      aSol[ey][ex][ip] = p_exact(cX[ex][icenter], cY[ey][icenter], 0.0);
    }
  }
  
  PetscCall(DMStagVecRestoreArray(dm, solLocal, &aSol));
  PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(DMLocalToGlobal(dm, solLocal, INSERT_VALUES, sol));
  PetscFunctionReturn(0);
}

// ============================================================================
// Perform one red-black Gauss-Seidel sweep
// Solve: (I/dt - nu*Laplace)*u^{n+1} + grad(p^{n+1}) = u^n/dt + f
//        div(u^{n+1}) = 0
// ============================================================================
PetscErrorCode GaussSeidelSweep(const DM &dm, Vec &sol, Vec &solLocal, 
                                const Vec &solOld, const Vec &solOldLocal,
                                const Vec &f, const Vec &fLocal, 
                                PetscInt Nx, PetscInt Ny,
                                PetscReal nu, PetscReal dt) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  const PetscReal ix = 1.0 / hx;
  const PetscReal iy = 1.0 / hy;
  const PetscReal ix2 = 1.0 / (hx * hx);
  const PetscReal iy2 = 1.0 / (hy * hy);
  const PetscReal coef = nu * dt;
  const PetscReal diag_u = 1.0 + 2.0 * coef * (ix2 + iy2);
  
  for (int color = 0; color < 2; ++color) {
    PetscScalar ***aSol, ***aSolOld, ***aF;
    PetscInt startx, starty, nx, ny, nEx[2];
    PetscInt iux, iuy, ip;
    
    PetscCall(DMGlobalToLocalBegin(dm, sol, INSERT_VALUES, solLocal));
    PetscCall(DMGlobalToLocalEnd(dm, sol, INSERT_VALUES, solLocal));
    PetscCall(DMGlobalToLocalBegin(dm, solOld, INSERT_VALUES, solOldLocal));
    PetscCall(DMGlobalToLocalEnd(dm, solOld, INSERT_VALUES, solOldLocal));
    PetscCall(DMGlobalToLocalBegin(dm, f, INSERT_VALUES, fLocal));
    PetscCall(DMGlobalToLocalEnd(dm, f, INSERT_VALUES, fLocal));
    PetscCall(DMStagVecGetArray(dm, solLocal, &aSol));
    PetscCall(DMStagVecGetArray(dm, solOldLocal, &aSolOld));
    PetscCall(DMStagVecGetArray(dm, fLocal, &aF));
    PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));
    
    // Update v velocity (DOWN edges)
    for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
      if (ey == 0 || ey == Ny) continue;  // Boundary
      
      for (PetscInt ex = startx; ex < startx + nx; ++ex) {
        if (((ex + ey) & 1) != color) continue;

        const PetscScalar vc = aSol[ey][ex][iuy];
        
        // Velocity neighbors for Laplacian
        PetscScalar vl = (ex == 0) ? -vc : aSol[ey][ex-1][iuy];
        PetscScalar vr = (ex == Nx-1) ? -vc : aSol[ey][ex+1][iuy];
        PetscScalar vd = (ey == 1) ? 0.0 : aSol[ey-1][ex][iuy];
        PetscScalar vu = (ey == Ny-1) ? 0.0 : aSol[ey+1][ex][iuy];
        
        // Pressure gradient: ∂p/∂y at DOWN edge (between cells ey-1 and ey)
        PetscScalar dpdy = 0.0;
        if (ey > 0 && ey < Ny) {
          dpdy = iy * (aSol[ey][ex][ip] - aSol[ey-1][ex][ip]);
        }

        const PetscScalar rhs = aSolOld[ey][ex][iuy] / dt + aF[ey][ex][iuy] - dpdy;
        const PetscScalar vnew = (rhs + coef * (ix2 * (vl + vr) + iy2 * (vd + vu))) / diag_u;
        aSol[ey][ex][iuy] = vnew;
      }
    }
    
    // Update u velocity (LEFT edges)
    for (PetscInt ey = starty; ey < starty + ny; ++ey) {
      for (PetscInt ex = startx; ex < startx + nx + nEx[0]; ++ex) {
        if (ex == 0 || ex == Nx) continue;  // Boundary
        if (((ex + ey) & 1) != color) continue;

        const PetscScalar uc = aSol[ey][ex][iux];
        
        // Velocity neighbors for Laplacian
        PetscScalar ul = (ex == 1) ? 0.0 : aSol[ey][ex-1][iux];
        PetscScalar ur = (ex == Nx-1) ? 0.0 : aSol[ey][ex+1][iux];
        PetscScalar ud = (ey == 0) ? -uc : aSol[ey-1][ex][iux];
        PetscScalar uu = (ey == Ny-1) ? -uc : aSol[ey+1][ex][iux];
        
        // Pressure gradient: ∂p/∂x at LEFT edge (between cells ex-1 and ex)
        PetscScalar dpdx = 0.0;
        if (ex > 0 && ex < Nx) {
          dpdx = ix * (aSol[ey][ex][ip] - aSol[ey][ex-1][ip]);
        }

        const PetscScalar rhs = aSolOld[ey][ex][iux] / dt + aF[ey][ex][iux] - dpdx;
        const PetscScalar unew = (rhs + coef * (ix2 * (ul + ur) + iy2 * (ud + uu))) / diag_u;
        aSol[ey][ex][iux] = unew;
      }
    }
    
    // Update pressure (element centers) to enforce incompressibility
    // ∇·u = 0 → ∂u/∂x + ∂v/∂y = 0
    for (PetscInt ey = starty; ey < starty + ny; ++ey) {
      for (PetscInt ex = startx; ex < startx + nx; ++ex) {
        if (((ex + ey) & 1) != color) continue;
        
        // Compute divergence at cell center
        const PetscScalar dudx = ix * (aSol[ey][ex+1][iux] - aSol[ey][ex][iux]);
        const PetscScalar dvdy = iy * (aSol[ey+1][ex][iuy] - aSol[ey][ex][iuy]);
        const PetscScalar div = dudx + dvdy;
        
        // Pressure correction (relaxation parameter)
        const PetscReal omega = 0.8;
        aSol[ey][ex][ip] -= omega * div / (2.0 * (ix2 + iy2));
      }
    }

    PetscCall(DMStagVecRestoreArray(dm, solLocal, &aSol));
    PetscCall(DMStagVecRestoreArray(dm, solOldLocal, &aSolOld));
    PetscCall(DMStagVecRestoreArray(dm, fLocal, &aF));
    PetscCall(DMLocalToGlobal(dm, solLocal, INSERT_VALUES, sol));
  }
  PetscFunctionReturn(0);
}

// ============================================================================
// Compute L2 errors
// ============================================================================
PetscErrorCode ComputeL2Error(const DM &dm, const Vec &sol, Vec &solLocal,
                              PetscInt Nx, PetscInt Ny, PetscReal t, PetscReal nu,
                              PetscReal *u_error, PetscReal *v_error, PetscReal *p_error) {
  PetscFunctionBeginUser;
  PetscScalar ***aSol;
  PetscScalar **cX, **cY;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt iprev, icenter, iux, iuy, ip;
  
  PetscCall(DMGlobalToLocalBegin(dm, sol, INSERT_VALUES, solLocal));
  PetscCall(DMGlobalToLocalEnd(dm, sol, INSERT_VALUES, solLocal));
  PetscCall(DMStagVecGetArray(dm, solLocal, &aSol));
  PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));
  PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_LEFT, &iprev));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN, 0, &iuy));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_LEFT, 0, &iux));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));

  PetscReal localError_u = 0.0, localError_v = 0.0, localError_p = 0.0;
  
  // v component error (DOWN edges, interior only)
  for (PetscInt ey = starty; ey < starty + ny + nEx[1]; ++ey) {
    if (ey == 0 || ey == Ny) continue;
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      const PetscScalar x = cX[ex][icenter];
      const PetscScalar y = cY[ey][iprev];
      const PetscScalar diff = aSol[ey][ex][iuy] - v_exact(x, y, t);
      localError_v += diff * diff;
    }
  }
  
  // u component error (LEFT edges, interior only)
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx + nEx[0]; ++ex) {
      if (ex == 0 || ex == Nx) continue;
      const PetscScalar x = cX[ex][iprev];
      const PetscScalar y = cY[ey][icenter];
      const PetscScalar diff = aSol[ey][ex][iux] - u_exact(x, y, t);
      localError_u += diff * diff;
    }
  }
  
  // Pressure error (element centers)
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      const PetscScalar x = cX[ex][icenter];
      const PetscScalar y = cY[ey][icenter];
      const PetscScalar diff = aSol[ey][ex][ip] - p_exact(x, y, t);
      localError_p += diff * diff;
    }
  }
  
  PetscCall(DMStagVecRestoreArray(dm, solLocal, &aSol));
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
  PetscCall(PetscInitialize(&argc, &argv, NULL, "Unsteady Stokes equation solver using DMStag"));

  DM dm;
  Vec sol, solLocal, solOld, solOldLocal, f, fLocal;
  PetscInt Nx = 32, Ny = 32;
  const PetscInt dof0 = 0, dof1 = 1, dof2 = 1;  // edge velocities + element pressure
  const PetscInt stencilWidth = 1;
  const PetscReal nu = 1.0;     // Kinematic viscosity
  PetscReal dt = 0.001;         // Time step
  PetscReal T_final = 0.01;     // Final time
  const PetscReal tol = 1e-8;   // Tolerance for GS iterations
  const PetscInt maxIts = 10000;
  int rank;
  PetscBool compute_error = PETSC_FALSE;
  PetscBool convergence_test = PETSC_FALSE;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-nx", &Nx, NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-ny", &Ny, NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-dt", &dt, NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-T", &T_final, NULL));
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
  PetscCall(DMCreateGlobalVector(dm, &sol));
  PetscCall(DMCreateGlobalVector(dm, &solOld));
  PetscCall(DMCreateGlobalVector(dm, &f));
  PetscCall(DMGetLocalVector(dm, &solLocal));
  PetscCall(DMGetLocalVector(dm, &solOldLocal));
  PetscCall(DMGetLocalVector(dm, &fLocal));
  
  PetscInt Nglob[2];
  PetscCall(DMStagGetGlobalSizes(dm, &Nglob[0], &Nglob[1], NULL));
  
  if (rank == 0) {
    std::cout << "=== Unsteady Stokes Equation Solver ===" << std::endl;
    std::cout << "Grid: " << Nglob[0] << " x " << Nglob[1] << std::endl;
    std::cout << "Viscosity: " << nu << ", dt: " << dt << ", T_final: " << T_final << std::endl;
    std::cout << "Time steps: " << (int)(T_final / dt) << std::endl;
  }

  // Setup initial condition
  PetscCall(SetupInitialCondition(dm, sol, solLocal));
  
  // Time stepping loop
  PetscReal t = 0.0;
  PetscInt step = 0;
  const PetscInt output_interval = std::max(1, (int)(T_final / dt / 10));
  
  while (t < T_final) {
    step++;
    t += dt;
    
    // sol_old = sol^n
    PetscCall(VecCopy(sol, solOld));
    
    // Setup RHS (force terms at time t)
    PetscCall(SetupRHS(dm, f, fLocal, t));
    
    // Solve coupled system using Gauss-Seidel
    for (PetscInt its = 1; its <= maxIts; ++its) {
      PetscCall(GaussSeidelSweep(dm, sol, solLocal, solOld, solOldLocal, f, fLocal, 
                                 Nglob[0], Nglob[1], nu, dt));
    }
    
    if (rank == 0 && (step % output_interval == 0 || step == 1)) {
      std::cout << "Step " << step << ", t=" << t << std::endl;
    }
  }
  
  if (rank == 0) {
    std::cout << "\nTime integration complete. Final time: " << t << std::endl;
  }
  
  // Compute errors
  if (compute_error || convergence_test) {
    PetscReal u_err, v_err, p_err;
    PetscCall(ComputeL2Error(dm, sol, solLocal, Nglob[0], Nglob[1], t, nu, &u_err, &v_err, &p_err));
    
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
  PetscCall(DMRestoreLocalVector(dm, &solLocal));
  PetscCall(DMRestoreLocalVector(dm, &solOldLocal));
  PetscCall(DMRestoreLocalVector(dm, &fLocal));
  PetscCall(VecDestroy(&sol));
  PetscCall(VecDestroy(&solOld));
  PetscCall(VecDestroy(&f));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFinalize());
  return 0;
}
