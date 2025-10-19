#include <petscdm.h>
#include <petscdmstag.h>
#include <petscksp.h>
#include <petscsys.h>

#include <cmath>
#include <iostream>

static inline PetscScalar u_exact(PetscScalar x, PetscScalar y) {
  return std::sin(M_PI * x) * std::sin(M_PI * y);
}

static inline PetscScalar rhs_f(PetscScalar x, PetscScalar y) {
  return 2.0 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
}

int main(int argc, char **argv) {
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, 
    "Matrix-free Gauss-Seidel for Poisson on DMStag (element-centered)"));

  DM dm;
  Vec u,uLocal, f, fLocal;            
  PetscInt Nx = 128, Ny = 128; 
  const PetscInt dof0 = 0, dof1 = 0, dof2 = 1;
  const PetscInt stencilWidth = 1;

  PetscReal tol = 1e-8;
  PetscInt maxIts = 40000;

  PetscCall(DMStagCreate2d(PETSC_COMM_WORLD,
                           DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                           Nx, Ny,
                           PETSC_DECIDE, PETSC_DECIDE,
                           dof0, dof1, dof2,
                           DMSTAG_STENCIL_BOX,
                           stencilWidth,
                           NULL, NULL,
                           &dm));
  PetscCall(DMSetFromOptions(dm));
  PetscCall(DMSetUp(dm));
  // Uniform coordinates [0,1] x [0,1]
  PetscCall(DMStagSetUniformCoordinatesProduct(dm, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0));

  // Create vectors
  PetscCall(DMCreateGlobalVector(dm, &u));
  PetscCall(DMCreateGlobalVector(dm, &f));
  PetscCall(DMGetLocalVector(dm, &uLocal));
  PetscCall(DMGetLocalVector(dm, &fLocal));
  PetscCall(VecSet(u, 0.0)); // initial guess satisfies Dirichlet 0 at boundary cells

  // Fill RHS f at element centers
  {
    PetscScalar ***aF;
    PetscScalar **cX, **cY;
    PetscInt startx, starty, nx, ny, nEx[2];
    PetscInt icenter, ip;
    PetscCall(DMGlobalToLocalBegin(dm, f, INSERT_VALUES, fLocal));
    PetscCall(DMGlobalToLocalEnd(dm, f, INSERT_VALUES, fLocal));
    PetscCall(DMStagVecGetArray(dm, fLocal, &aF));
    PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));
    PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
    PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));
    for (PetscInt ey = starty; ey < starty + ny; ++ey) {
      for (PetscInt ex = startx; ex < startx + nx; ++ex) {
        const PetscScalar x = cX[ex][icenter];
        const PetscScalar y = cY[ey][icenter];
        aF[ey][ex][ip] = rhs_f(x, y);
        // printf("RHS f at element (%d,%d) (x=%.8f,y=%.8f): %.6f\n", ex, ey, x, y, aF[ey][ex][icenter]);
      }
    }
    PetscCall(DMStagVecRestoreArray(dm, fLocal, &aF));
    PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
    PetscCall(DMLocalToGlobal(dm, fLocal, INSERT_VALUES, f));
  }

  // Discretization constants
  PetscInt Nglob[2];
  PetscCall(DMStagGetGlobalSizes(dm, &Nglob[0], &Nglob[1], NULL));
  const PetscReal hx = 1.0 / Nglob[0];
  const PetscReal hy = 1.0 / Nglob[1];
  const PetscReal ix2 = 1.0 / (hx * hx);
  const PetscReal iy2 = 1.0 / (hy * hy);
  const PetscReal diag = 2.0 * ix2 + 2.0 * iy2;

  // Iteration
  PetscInt its = 0;
  PetscReal resNorm = 0.0;
  PetscReal res0 = -1.0;
  const PetscReal atol = tol;
  int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // Helper: compute residual norm over interior cells
  auto compute_residual_norm = [&](PetscReal &outNorm) -> void {
    PetscScalar ***aU, ***aF;
    PetscInt startx, starty, nx, ny, nEx[2];
    PetscInt icenter;
    PetscCallAbort(PETSC_COMM_WORLD, DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
    PetscCallAbort(PETSC_COMM_WORLD, DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
    PetscCallAbort(PETSC_COMM_WORLD, DMGlobalToLocalBegin(dm, f, INSERT_VALUES, fLocal));
    PetscCallAbort(PETSC_COMM_WORLD, DMGlobalToLocalEnd(dm, f, INSERT_VALUES, fLocal));
    PetscCallAbort(PETSC_COMM_WORLD, DMStagVecGetArray(dm, uLocal, &aU));
    PetscCallAbort(PETSC_COMM_WORLD, DMStagVecGetArray(dm, fLocal, &aF));
    PetscCallAbort(PETSC_COMM_WORLD, DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));
    PetscCallAbort(PETSC_COMM_WORLD, DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &icenter));
    PetscReal localSum = 0.0;
    for (PetscInt ey = starty; ey < starty + ny; ++ey) {
      for (PetscInt ex = startx; ex < startx + nx; ++ex) {
        // For cell-centered grid, all cells are interior
        const PetscScalar uc = aU[ey][ex][icenter];
        
        // Get neighbor values with ghost cell BC
        PetscScalar ul, ur, ud, uu;
        if (ex == 0) ul = -uc; else ul = aU[ey][ex-1][icenter];
        if (ex == Nglob[0]-1) ur = -uc; else ur = aU[ey][ex+1][icenter];
        if (ey == 0) ud = -uc; else ud = aU[ey-1][ex][icenter];
        if (ey == Nglob[1]-1) uu = -uc; else uu = aU[ey+1][ex][icenter];
        
        const PetscScalar Au = diag * uc - ix2*(ul + ur) - iy2*(ud + uu);
        const PetscScalar ff = aF[ey][ex][icenter];
        const PetscScalar r = ff - Au; // residual for -Laplace(u)=f
        localSum += PetscRealPart(r*r);
      }
    }
    PetscCallAbort(PETSC_COMM_WORLD, DMStagVecRestoreArray(dm, uLocal, &aU));
    PetscCallAbort(PETSC_COMM_WORLD, DMStagVecRestoreArray(dm, fLocal, &aF));
    MPI_Allreduce(&localSum, &outNorm, 1, MPIU_REAL, MPIU_SUM, PETSC_COMM_WORLD);
    outNorm = std::sqrt(outNorm);
  };

  // Initial residual
  compute_residual_norm(resNorm);
  res0 = resNorm;
  if (rank == 0) std::cout << "GS start: ||r|| = " << resNorm << std::endl;

  // Red-black Gauss-Seidel iterations
  for (its = 1; its <= maxIts; ++its) {
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
      PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));
      PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &icenter));
      // printf("Processing color %d, nx=%d, ny=%d, Nglob[0]=%d element (%d, %d)\n", color, nx, ny, nEx[0], startx, starty);
      for (PetscInt ey = starty; ey < starty + ny; ++ey) {
        for (PetscInt ex = startx; ex < startx + nx; ++ex) {
          if ( ((ex + ey) & 1) != color ) continue;

          PetscScalar ul, ur, ud, uu;
          
          if (ex == 0) {
            ul = -aU[ey][ex][icenter]; // Ghost cell: u_{-1} = -u_0 for Dirichlet BC u=0
          } else {
            ul = aU[ey][ex-1][icenter];
          }
          
          if (ex == Nglob[0]-1) {
            ur = -aU[ey][ex][icenter]; // Ghost cell: u_{N} = -u_{N-1}
          } else {
            ur = aU[ey][ex+1][icenter];
          }
          
          if (ey == 0) {
            ud = -aU[ey][ex][icenter]; // Ghost cell: u_{-1} = -u_0
          } else {
            ud = aU[ey-1][ex][icenter];
          }
          
          if (ey == Nglob[1]-1) {
            uu = -aU[ey][ex][icenter]; // Ghost cell: u_{N} = -u_{N-1}
          } else {
            uu = aU[ey+1][ex][icenter];
          }
          
          const PetscScalar ff = aF[ey][ex][icenter];
          const PetscScalar unew = (ix2*(ul + ur) + iy2*(ud + uu) + ff) / diag;
          aU[ey][ex][icenter] = unew;
        }
      }

      PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
      PetscCall(DMStagVecRestoreArray(dm, fLocal, &aF));
      // scatter updated owned values back to global, then refresh ghosts for next color
      PetscCall(DMLocalToGlobal(dm, uLocal, INSERT_VALUES, u));
    }

    compute_residual_norm(resNorm);
    if (rank == 0 && (its % 5000 == 0 || resNorm <= atol)) {
      std::cout << "it=" << its << ", ||r||=" << resNorm << " (rel " << resNorm/(res0>0?res0:1.0) << ")" << std::endl;
    }
    if (resNorm <= atol) break;
  }

  if (rank == 0) {
    std::cout << "Done. it=" << its << ", ||r||=" << resNorm << std::endl;
  }

  // Optionally compute error vs exact solution
  PetscBool compute_error = PETSC_FALSE;
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-poisson_check_error", &compute_error, NULL));
  if (compute_error) {
    Vec uExact, diff;
    PetscCall(DMCreateGlobalVector(dm, &uExact));
    PetscCall(DMGetLocalVector(dm, &uLocal));
    PetscScalar ***aUe;
    PetscScalar **cX, **cY;
    PetscInt startx, starty, nx, ny, nEx[2];
    PetscInt icenter,ip;
    PetscCall(DMGlobalToLocalBegin(dm, uExact, INSERT_VALUES, uLocal));
    PetscCall(DMGlobalToLocalEnd(dm, uExact, INSERT_VALUES, uLocal));
    PetscCall(DMStagVecGetArray(dm, uLocal, &aUe));
    PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL, &nEx[0], &nEx[1], NULL));
    PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
    // PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &icenter));
    PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));

    for (PetscInt ey = starty; ey < starty + ny; ++ey) {
      for (PetscInt ex = startx; ex < startx + nx; ++ex) {
        const PetscScalar x = cX[ex][icenter];
        const PetscScalar y = cY[ey][icenter];
        // printf("Exact u at element (%d,%d) (x=%.8f,y=%.8f): %.6f\n", ex, ey, x, y, u_exact(x, y));
        aUe[ey][ex][ip] = u_exact(x,y);
      }
    }
    PetscCall(DMStagVecRestoreArray(dm, uLocal, &aUe));
    PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
    PetscCall(DMLocalToGlobal(dm, uLocal, INSERT_VALUES, uExact));

    PetscCall(VecDuplicate(u, &diff));
    PetscCall(VecCopy(u, diff));
    PetscCall(VecAXPY(diff, -1.0, uExact));
    PetscReal nrm2; PetscCall(VecNorm(diff, NORM_2, &nrm2));
    if (rank == 0) std::cout << "||u - u_exact||_2 = " << nrm2 * std::sqrt(hx * hy)<< std::endl;
    PetscCall(VecDestroy(&diff));
    PetscCall(VecDestroy(&uExact));
  }

  PetscCall(DMRestoreLocalVector(dm, &uLocal));
  PetscCall(DMRestoreLocalVector(dm, &fLocal));
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&f));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFinalize());
  return 0;
}


// it=10000, ||r||=1.40243e-08 (rel 2.22042e-11)
// Done. it=10001, ||r||=1.40243e-08
// ||u - u_exact||_2 = 1.44212