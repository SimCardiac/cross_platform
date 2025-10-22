#include <petscdm.h>
#include <petscdmstag.h>
#include <petscksp.h>
#include <petscsys.h>

#include <cmath>
#include <iostream>

#include "DMStag_boundary_helpers.h"


namespace FUNC {
  using EXACT   = PetscScalar(*)(PetscScalar, PetscScalar, PetscScalar, PetscScalar);
  using INITIAL = PetscScalar(*)(PetscScalar, PetscScalar);
  using RHS     = PetscScalar(*)(PetscScalar, PetscScalar, PetscScalar);
}

static inline PetscScalar u_exact(PetscScalar x, PetscScalar y, PetscScalar t, PetscScalar alpha) {
  return std::exp(-2.0 * M_PI * M_PI * alpha * t) * std::sin(M_PI * x) * std::sin(M_PI * y);
}

static inline PetscScalar u_initial(PetscScalar x, PetscScalar y) {
  return std::sin(M_PI * x) * std::sin(M_PI * y);
}

static inline PetscScalar rhs_f(PetscScalar x, PetscScalar y, PetscScalar t) {
  return 0.0;
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
                                   PetscReal *residualNorm) {
  PetscFunctionBeginUser;
  PetscScalar ***aU, ***aUold, ***aF;
  PetscInt startx, starty, nx, ny, nEx[2];
  PetscInt icenter;
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
  /* enforce Dirichlet ghost values for u */
  PetscCall(DMStagVecGetArray(dm, uOldLocal, &aUold));
  PetscCall(DMStagVecGetArray(dm, fLocal, &aF));
  PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &icenter));
  
  PetscReal localSum = 0.0;
  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      const PetscScalar uc = aU[ey][ex][icenter];

      // Get neighbor values safely (reflecting Dirichlet at boundaries)
      const PetscScalar ul = DMStag_GetLeft(aU, ex, ey, icenter, Nx);
      const PetscScalar ur = DMStag_GetRight(aU, ex, ey, icenter, Nx);
      const PetscScalar ud = DMStag_GetDown(aU, ex, ey, icenter, Ny);
      const PetscScalar uu = DMStag_GetUp(aU, ex, ey, icenter, Ny);

      const PetscScalar lhs = diag * uc - coef * (ix2 * (ul + ur) + iy2 * (ud + uu));
      const PetscScalar rhs = aUold[ey][ex][icenter] + dt * aF[ey][ex][icenter];
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
    PetscInt icenter;
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
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &icenter));
    
    for (PetscInt ey = starty; ey < starty + ny; ++ey) {
      for (PetscInt ex = startx; ex < startx + nx; ++ex) {
        if (((ex + ey) & 1) != color)
          continue;

        const PetscScalar uc = aU[ey][ex][icenter];
        
        // Get neighbor values safely (reflecting Dirichlet at boundaries)
        const PetscScalar ul = DMStag_GetLeft(aU, ex, ey, icenter, Nx);
        const PetscScalar ur = DMStag_GetRight(aU, ex, ey, icenter, Nx);
        const PetscScalar ud = DMStag_GetDown(aU, ex, ey, icenter, Ny);
        const PetscScalar uu = DMStag_GetUp(aU, ex, ey, icenter, Ny);

        // RHS = u_old + dt * f
        const PetscScalar rhs = aUold[ey][ex][icenter] + dt * aF[ey][ex][icenter];
        // u^{n+1} = [rhs + dt*alpha*(1/hx^2*(ul+ur) + 1/hy^2*(ud+uu))] / diag
        const PetscScalar unew = (rhs + coef * (ix2 * (ul + ur) + iy2 * (ud + uu))) / diag;
        aU[ey][ex][icenter] = unew;
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
                              PetscReal *l2Error,FUNC::EXACT _u_exact) {
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
  PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));
  PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));

  for (PetscInt ey = starty; ey < starty + ny; ++ey) {
    for (PetscInt ex = startx; ex < startx + nx; ++ex) {
      const PetscScalar x = cX[ex][icenter];
      const PetscScalar y = cY[ey][icenter];
      aUe[ey][ex][ip] = _u_exact(x, y, t, alpha);
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
  *l2Error = nrm2 * std::sqrt(hx * hy);
  PetscCall(VecDestroy(&diff));
  PetscCall(VecDestroy(&uExact));
  PetscFunctionReturn(0);
}

int main(int argc, char **argv) {
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, "Heat equation solver using implicit Euler with DMStag"));

  DM dm;
  Vec u, uLocal, uOld, uOldLocal, f, fLocal;
  PetscInt Nx = 64, Ny = 64;
  const PetscInt dof0 = 0, dof1 = 0, dof2 = 1; 
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
