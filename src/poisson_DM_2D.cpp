#include <petscdm.h>
#include <petscdmda.h>
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
    "Matrix-free Gauss-Seidel for Poisson on DMDA (cell-centered)"));

  DM dm;
  Vec u,uLocal, f, fLocal;            
  PetscInt Nx = 32, Ny = 32; 
  const PetscInt dof = 1;
  const PetscInt stencilWidth = 1;

  PetscReal tol = 1e-8;
  PetscInt maxIts = 40000;

  PetscCall(DMDACreate2d(PETSC_COMM_WORLD,
                         DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                         DMDA_STENCIL_STAR,
                         Nx, Ny,
                         PETSC_DECIDE, PETSC_DECIDE,
                         dof, stencilWidth,
                         NULL, NULL,
                         &dm));
  PetscCall(DMSetFromOptions(dm));
  PetscCall(DMSetUp(dm));
  // Uniform coordinates [0,1] x [0,1]
  PetscCall(DMDASetUniformCoordinates(dm, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0));

  // Create vectors
  PetscCall(DMCreateGlobalVector(dm, &u));
  PetscCall(DMCreateGlobalVector(dm, &f));
  PetscCall(DMCreateLocalVector(dm, &uLocal));
  PetscCall(DMCreateLocalVector(dm, &fLocal));
  PetscCall(VecSet(u, 0.0)); // initial guess with Dirichlet BC u=0 on boundary

  // Fill RHS f at cell centers
  {
    DM dmCoord;
    Vec coords;
    PetscScalar **coordArr;
    PetscScalar **fArr;
    PetscInt xs, ys, xm, ym;
    
    PetscCall(DMGetCoordinateDM(dm, &dmCoord));
    PetscCall(DMGetCoordinates(dm, &coords));
    PetscCall(DMDAVecGetArray(dmCoord, coords, &coordArr));
    PetscCall(DMDAVecGetArray(dm, f, &fArr));
    PetscCall(DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL));
    
    for (PetscInt j = ys; j < ys + ym; j++) {
      for (PetscInt i = xs; i < xs + xm; i++) {
        // coordArr is indexed as [j][i*2+0] for x, [j][i*2+1] for y
        const PetscScalar x = coordArr[j][i*2];
        const PetscScalar y = coordArr[j][i*2+1];
        fArr[j][i] = rhs_f(x, y);
        printf("RHS f at (%d,%d) (x=%.3f,y=%.3f): %.6f\n", i, j, x, y, fArr[j][i]);
      }
    }
    
    PetscCall(DMDAVecRestoreArray(dmCoord, coords, &coordArr));
    PetscCall(DMDAVecRestoreArray(dm, f, &fArr));
  }

  // Get grid info
  PetscInt Mx, My;
  PetscCall(DMDAGetInfo(dm, NULL, &Mx, &My, NULL, NULL, NULL, NULL, 
                        NULL, NULL, NULL, NULL, NULL, NULL));
  const PetscReal hx = 1.0 / (Mx-1);
  const PetscReal hy = 1.0 / (My-1);
  const PetscReal ix2 = 1.0 / (hx * hx);
  const PetscReal iy2 = 1.0 / (hy * hy);
  const PetscReal diag = 2.0 * ix2 + 2.0 * iy2;

  // Iteration state
  PetscInt its = 0;
  PetscReal resNorm = 0.0, res0 = -1.0;
  const PetscReal atol = tol;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // Helper: compute residual norm over interior cells
  auto compute_residual_norm = [&](PetscReal &outNorm) -> void {
    PetscScalar **uArr, **fArr;
    PetscInt xs, ys, xm, ym;
    
    PetscCallAbort(PETSC_COMM_WORLD, DMGlobalToLocal(dm, u, INSERT_VALUES, uLocal));
    PetscCallAbort(PETSC_COMM_WORLD, DMGlobalToLocal(dm, f, INSERT_VALUES, fLocal));
    PetscCallAbort(PETSC_COMM_WORLD, DMDAVecGetArrayRead(dm, uLocal, &uArr));
    PetscCallAbort(PETSC_COMM_WORLD, DMDAVecGetArrayRead(dm, fLocal, &fArr));
    PetscCallAbort(PETSC_COMM_WORLD, DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL));
    
    PetscReal localSum = 0.0;
    // printf("Rank %d computing residual on local corners xs=%d ys=%d xm=%d ym=%d\n", 
    //        rank, xs, ys, xm, ym);
    for (PetscInt j = ys; j < ys + ym; j++) {
      for (PetscInt i = xs; i < xs + xm; i++) {
        // Skip boundary cells (Dirichlet u=0)
        const bool interior = (i > 0 && i < Mx-1 && j > 0 && j < My-1);
        if (!interior) continue;
        
        const PetscScalar uc = uArr[j][i];
        const PetscScalar ul = uArr[j][i-1];
        const PetscScalar ur = uArr[j][i+1];
        const PetscScalar ud = uArr[j-1][i];
        const PetscScalar uu = uArr[j+1][i];
        const PetscScalar Au = diag * uc - ix2*(ul + ur) - iy2*(ud + uu);
        const PetscScalar ff = fArr[j][i];
        const PetscScalar r = ff - Au; // residual for -Laplace(u)=f
        localSum += PetscRealPart(r * r);
      }
    }
    
    PetscCallAbort(PETSC_COMM_WORLD, DMDAVecRestoreArrayRead(dm, uLocal, &uArr));
    PetscCallAbort(PETSC_COMM_WORLD, DMDAVecRestoreArrayRead(dm, fLocal, &fArr));
    MPI_Allreduce(&localSum, &outNorm, 1, MPIU_REAL, MPIU_SUM, PETSC_COMM_WORLD);
    outNorm = std::sqrt(outNorm);
  };

  // Initial residual
  compute_residual_norm(resNorm);
  res0 = resNorm;
  if (rank == 0) {
    std::cout << "GS start: ||r|| = " << resNorm << std::endl;
  }

  // Red-black Gauss-Seidel iterations
  for (its = 1; its <= maxIts; its++) {
    for (int color = 0; color < 2; color++) {
      PetscScalar **uArr, **fArr;
      PetscInt xs, ys, xm, ym;
      
      PetscCall(DMGlobalToLocal(dm, u, INSERT_VALUES, uLocal));
      PetscCall(DMGlobalToLocal(dm, f, INSERT_VALUES, fLocal));
      PetscCall(DMDAVecGetArray(dm, uLocal, &uArr));
      PetscCall(DMDAVecGetArrayRead(dm, fLocal, &fArr));
      PetscCall(DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL));
      
      for (PetscInt j = ys; j < ys + ym; j++) {
        for (PetscInt i = xs; i < xs + xm; i++) {
          // Skip boundary (Dirichlet u=0)
          const bool interior = (i > 0 && i < Mx-1 && j > 0 && j < My-1);
          if (!interior) continue;
          
          // Red-black checkerboard pattern
          if (((i + j) & 1) != color) continue;
          
          const PetscScalar ul = uArr[j][i-1];
          const PetscScalar ur = uArr[j][i+1];
          const PetscScalar ud = uArr[j-1][i];
          const PetscScalar uu = uArr[j+1][i];
          const PetscScalar ff = fArr[j][i];
          const PetscScalar unew = (ix2*(ul + ur) + iy2*(ud + uu) + ff) / diag;
          uArr[j][i] = unew;
        }
      }
      
      PetscCall(DMDAVecRestoreArray(dm, uLocal, &uArr));
      PetscCall(DMDAVecRestoreArrayRead(dm, fLocal, &fArr));
      PetscCall(DMLocalToGlobal(dm, uLocal, INSERT_VALUES, u));
    }
    
    // Check convergence every 10 iterations
    if (its % 5000 == 0 || its == 1) {
      compute_residual_norm(resNorm);
      if (rank == 0) {
        std::cout << "it=" << its << ", ||r||=" << resNorm 
                  << " (rel " << resNorm/(res0>0?res0:1.0) << ")" << std::endl;
      }
      if (resNorm <= atol) break;
    }
  }

  // Final residual
  compute_residual_norm(resNorm);
  if (rank == 0) {
    std::cout << "Done. it=" << its << ", ||r||=" << resNorm << std::endl;
  }

  // Optionally compute error vs exact solution
  PetscBool compute_error = PETSC_FALSE;
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-poisson_check_error", &compute_error, NULL));
  if (compute_error) {
    Vec uExact, diff;
    DM dmCoord;
    Vec coords;
    PetscScalar **coordArr, **uExactArr;
    PetscInt xs, ys, xm, ym;
    
    PetscCall(DMCreateGlobalVector(dm, &uExact));
    PetscCall(DMGetCoordinateDM(dm, &dmCoord));
    PetscCall(DMGetCoordinates(dm, &coords));
    PetscCall(DMDAVecGetArrayRead(dmCoord, coords, &coordArr));
    PetscCall(DMDAVecGetArray(dm, uExact, &uExactArr));
    PetscCall(DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL));
    
    for (PetscInt j = ys; j < ys + ym; j++) {
      for (PetscInt i = xs; i < xs + xm; i++) {
        const PetscScalar x = coordArr[j][i*2];
        const PetscScalar y = coordArr[j][i*2+1];
        uExactArr[j][i] = u_exact(x, y);
      }
    }
    
    PetscCall(DMDAVecRestoreArrayRead(dmCoord, coords, &coordArr));
    PetscCall(DMDAVecRestoreArray(dm, uExact, &uExactArr));
    
    PetscCall(VecDuplicate(u, &diff));
    PetscCall(VecCopy(u, diff));
    PetscCall(VecAXPY(diff, -1.0, uExact));
    
    // Compute discrete L2 norm: sqrt(sum |diff|^2 * hx * hy)
    // which is sqrt(hx*hy) * ||diff||_l2
    PetscReal nrm2_discrete;
    PetscCall(VecNorm(diff, NORM_2, &nrm2_discrete));
    const PetscReal nrm2_L2 = std::sqrt(hx * hy) * nrm2_discrete;
    
    if (rank == 0) {
      std::cout << "||u - u_exact||_l2 (discrete) = " << nrm2_discrete << std::endl;
      std::cout << "||u - u_exact||_L2 (continuous) = " << nrm2_L2 << std::endl;
      std::cout << "Grid: " << Mx << " x " << My << ", h = (" << hx << ", " << hy << ")" << std::endl;
    }
    
    PetscCall(VecDestroy(&diff));
    PetscCall(VecDestroy(&uExact));
  }

  // Cleanup
  PetscCall(VecDestroy(&uLocal));
  PetscCall(VecDestroy(&fLocal));
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&f));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFinalize());
  return 0;
}

// Done. it=5000, ||r||=4.29277e-12
// ||u - u_exact||_l2 (discrete) = 0.0132724
// ||u - u_exact||_L2 (continuous) = 0.000428142
// Grid: 32 x 32, h = (0.0322581, 0.0322581)

// Done. it=15000, ||r||=3.15553e-11
// ||u - u_exact||_l2 (discrete) = 0.00652833
// ||u - u_exact||_L2 (continuous) = 0.000103624
// Grid: 64 x 64, h = (0.015873, 0.015873)

// Done. it=40001, ||r||=4.14498e-08
// ||u - u_exact||_l2 (discrete) = 0.00323816
// ||u - u_exact||_L2 (continuous) = 2.54973e-05
// Grid: 128 x 128, h = (0.00787402, 0.00787402)