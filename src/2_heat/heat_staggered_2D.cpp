#include <petscdm.h>
#include <petscdmstag.h>
#include <petscksp.h>
#include <cmath>
#include <iostream>

#include "common/boundary.h"
#include "common/mms.h"

// ============================================================================
// Staggered (mixed first-order) heat equation on DMStag (Item 6)
//
// First-order system:
//   gx - ∂u/∂x = 0       (left faces)
//   gy - ∂u/∂y = 0       (down faces)
//   ∂u/∂t - α(∂gx/∂x + ∂gy/∂y) = 0   (elements)
//
// Implicit Euler: [ I     0     -Gx      ] [gx^{n+1}]   [0    ]
//                 [ 0     I     -Gy      ] [gy^{n+1}] = [0    ]
//                 [-αΔtDx -αΔtDy  I     ] [u^{n+1}]    [u^n  ]
//
// Manufactured solution (SINPI_SCALAR_2D):
//   u = exp(-2π²αt) sin(πx) sin(πy)
// ============================================================================

static const PetscScalar alpha = 0.1;

// DOF indexing (same as staggered Poisson)
static inline PetscInt idx_u(PetscInt ex, PetscInt ey, PetscInt Nx, PetscInt Ny) {
  (void)Ny; return ey * Nx + ex;
}
static inline PetscInt idx_gx(PetscInt fx, PetscInt fy, PetscInt Nx, PetscInt Ny) {
  return Nx * Ny + fy * (Nx + 1) + fx;
}
static inline PetscInt idx_gy(PetscInt fx, PetscInt fy, PetscInt Nx, PetscInt Ny) {
  return Nx * Ny + (Nx + 1) * Ny + fx * (Ny + 1) + fy;
}
static inline PetscInt totalDOFs(PetscInt Nx, PetscInt Ny) {
  return Nx * Ny + (Nx + 1) * Ny + Nx * (Ny + 1);
}

// ============================================================================
// Assemble the monolithic system for one time step
//   [ I     0     -Gx      ] [gx]
//   [ 0     I     -Gy      ] [gy]
//   [-c*Dx  -c*Dy    I     ] [u ]    where c = alpha*dt
// ============================================================================
static PetscErrorCode AssembleSystem(Mat A, PetscInt Nx, PetscInt Ny, PetscReal dt,
                                      const BoundaryCondition &bc) {
  PetscFunctionBeginUser;
  const PetscReal h = 1.0 / Nx;
  const PetscReal c = alpha * dt;
  const PetscReal ih = 1.0 / h, ih2 = 2.0 / h;
  const PetscReal cih = c * ih;

  PetscInt rstart, rend;
  PetscCall(MatGetOwnershipRange(A, &rstart, &rend));

  // ---- gx equations: gx - Gx u = BC terms ----
  for (PetscInt fy = 0; fy < Ny; ++fy) {
    for (PetscInt fx = 0; fx <= Nx; ++fx) {
      PetscInt row = idx_gx(fx, fy, Nx, Ny);
      if (row < rstart || row >= rend) continue;
      if (fx == 0) {
        if (bc.left == BC_NEUMANN) {
          PetscCall(MatSetValue(A, row, row, 1.0, INSERT_VALUES));  // gx = -g
        } else {
          PetscCall(MatSetValue(A, row, row, 1.0, INSERT_VALUES));
          PetscCall(MatSetValue(A, row, idx_u(0, fy, Nx, Ny), -ih2, INSERT_VALUES));
        }
      } else if (fx == Nx) {
        if (bc.right == BC_NEUMANN) {
          PetscCall(MatSetValue(A, row, row, 1.0, INSERT_VALUES));  // gx = +g
        } else {
          PetscCall(MatSetValue(A, row, row, 1.0, INSERT_VALUES));
          PetscCall(MatSetValue(A, row, idx_u(Nx - 1, fy, Nx, Ny), ih2, INSERT_VALUES));
        }
      } else {
        PetscCall(MatSetValue(A, row, row, 1.0, INSERT_VALUES));
        PetscCall(MatSetValue(A, row, idx_u(fx, fy, Nx, Ny), -ih, INSERT_VALUES));
        PetscCall(MatSetValue(A, row, idx_u(fx - 1, fy, Nx, Ny), ih, INSERT_VALUES));
      }
    }
  }

  // ---- gy equations: gy - Gy u = BC terms ----
  for (PetscInt fx = 0; fx < Nx; ++fx) {
    for (PetscInt fy = 0; fy <= Ny; ++fy) {
      PetscInt row = idx_gy(fx, fy, Nx, Ny);
      if (row < rstart || row >= rend) continue;
      if (fy == 0) {
        if (bc.bottom == BC_NEUMANN) {
          PetscCall(MatSetValue(A, row, row, 1.0, INSERT_VALUES));  // gy = -g
        } else {
          PetscCall(MatSetValue(A, row, row, 1.0, INSERT_VALUES));
          PetscCall(MatSetValue(A, row, idx_u(fx, 0, Nx, Ny), -ih2, INSERT_VALUES));
        }
      } else if (fy == Ny) {
        if (bc.top == BC_NEUMANN) {
          PetscCall(MatSetValue(A, row, row, 1.0, INSERT_VALUES));  // gy = +g
        } else {
          PetscCall(MatSetValue(A, row, row, 1.0, INSERT_VALUES));
          PetscCall(MatSetValue(A, row, idx_u(fx, Ny - 1, Nx, Ny), ih2, INSERT_VALUES));
        }
      } else {
        PetscCall(MatSetValue(A, row, row, 1.0, INSERT_VALUES));
        PetscCall(MatSetValue(A, row, idx_u(fx, fy, Nx, Ny), -ih, INSERT_VALUES));
        PetscCall(MatSetValue(A, row, idx_u(fx, fy - 1, Nx, Ny), ih, INSERT_VALUES));
      }
    }
  }

  // ---- u equations: u - c*div(g) = u_old + dt*f ----
  for (PetscInt ey = 0; ey < Ny; ++ey) {
    for (PetscInt ex = 0; ex < Nx; ++ex) {
      PetscInt row = idx_u(ex, ey, Nx, Ny);
      if (row < rstart || row >= rend) continue;
      PetscCall(MatSetValue(A, row, row, 1.0, INSERT_VALUES));
      PetscCall(MatSetValue(A, row, idx_gx(ex + 1, ey, Nx, Ny), -cih, INSERT_VALUES));
      PetscCall(MatSetValue(A, row, idx_gx(ex, ey, Nx, Ny), cih, INSERT_VALUES));
      PetscCall(MatSetValue(A, row, idx_gy(ex, ey + 1, Nx, Ny), -cih, INSERT_VALUES));
      PetscCall(MatSetValue(A, row, idx_gy(ex, ey, Nx, Ny), cih, INSERT_VALUES));
    }
  }

  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  PetscFunctionReturn(0);
}

// ============================================================================
// Set initial condition: u at elements, gx/gy from MMS gradients
// ============================================================================
static PetscErrorCode SetInitialCondition(Vec x, PetscInt Nx, PetscInt Ny,
                                           const ManufacturedSolution &mms) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;
  PetscReal tf=1.0;
  if(mms.u_td){PetscScalar us=mms.u(0.25,0.25); if(us!=0) tf=mms.u_td(0.25,0.25,0,alpha)/us;}

  PetscInt rstart, rend;
  PetscCall(VecGetOwnershipRange(x, &rstart, &rend));

  for (PetscInt ey = 0; ey < Ny; ++ey) {
    for (PetscInt ex = 0; ex < Nx; ++ex) {
      PetscInt row = idx_u(ex, ey, Nx, Ny);
      if (row < rstart || row >= rend) continue;
      const PetscScalar xc = (ex + 0.5) * hx;
      const PetscScalar yc = (ey + 0.5) * hy;
      PetscCall(VecSetValue(x, row, mms.u_td(xc, yc, 0, alpha), INSERT_VALUES));
    }
  }

  for (PetscInt fy = 0; fy < Ny; ++fy) {
    for (PetscInt fx = 0; fx <= Nx; ++fx) {
      PetscInt row = idx_gx(fx, fy, Nx, Ny);
      if (row < rstart || row >= rend) continue;
      const PetscScalar xf = fx * hx;
      const PetscScalar yc = (fy + 0.5) * hy;
      PetscCall(VecSetValue(x, row, tf * mms.ux(xf, yc), INSERT_VALUES));
    }
  }

  for (PetscInt fx = 0; fx < Nx; ++fx) {
    for (PetscInt fy = 0; fy <= Ny; ++fy) {
      PetscInt row = idx_gy(fx, fy, Nx, Ny);
      if (row < rstart || row >= rend) continue;
      const PetscScalar xc = (fx + 0.5) * hx;
      const PetscScalar yf = fy * hy;
      PetscCall(VecSetValue(x, row, tf * mms.uy(xc, yf), INSERT_VALUES));
    }
  }

  PetscCall(VecAssemblyBegin(x));
  PetscCall(VecAssemblyEnd(x));
  PetscFunctionReturn(0);
}

// ============================================================================
// Build RHS: BC contributions + u_old + dt*f
// ============================================================================
static PetscErrorCode BuildRHS(Vec b, const Vec xOld, PetscInt Nx, PetscInt Ny,
                                PetscReal t, PetscReal dt,
                                const ManufacturedSolution &mms,
                                const BoundaryCondition &bc) {
  PetscFunctionBeginUser;
  const PetscReal h = 1.0 / Nx, c = alpha * dt;
  PetscReal t_eval = t + dt, tf = 1.0;
  if (mms.u_td) { PetscScalar us = mms.u(0.25, 0.25); if (us != 0) tf = mms.u_td(0.25, 0.25, t_eval, alpha) / us; }
  PetscCall(VecSet(b, 0.0));

  PetscInt rstart, rend;
  PetscCall(VecGetOwnershipRange(b, &rstart, &rend));

  // gx RHS: boundary contributions
  for (PetscInt fy = 0; fy < Ny; ++fy) {
    { PetscInt row = idx_gx(0, fy, Nx, Ny); if (row >= rstart && row < rend) {
        PetscScalar v = 0;
        if (bc.left == BC_NEUMANN) v = -tf * EvalBC(bc.g_left, fy * h + 0.5 * h);
        else                       v = -2.0 * tf * EvalBC(bc.g_left, fy * h + 0.5 * h) / h;
        PetscCall(VecSetValue(b, row, v, INSERT_VALUES));
    }}
    { PetscInt row = idx_gx(Nx, fy, Nx, Ny); if (row >= rstart && row < rend) {
        PetscScalar v = 0;
        if (bc.right == BC_NEUMANN) v = tf * EvalBC(bc.g_right, fy * h + 0.5 * h);
        else                        v = 2.0 * tf * EvalBC(bc.g_right, fy * h + 0.5 * h) / h;
        PetscCall(VecSetValue(b, row, v, INSERT_VALUES));
    }}
  }

  // gy RHS: boundary contributions
  for (PetscInt fx = 0; fx < Nx; ++fx) {
    { PetscInt row = idx_gy(fx, 0, Nx, Ny); if (row >= rstart && row < rend) {
        PetscScalar v = 0;
        if (bc.bottom == BC_NEUMANN) v = -tf * EvalBC(bc.g_bottom, fx * h + 0.5 * h);
        else                         v = -2.0 * tf * EvalBC(bc.g_bottom, fx * h + 0.5 * h) / h;
        PetscCall(VecSetValue(b, row, v, INSERT_VALUES));
    }}
    { PetscInt row = idx_gy(fx, Ny, Nx, Ny); if (row >= rstart && row < rend) {
        PetscScalar v = 0;
        if (bc.top == BC_NEUMANN) v = tf * EvalBC(bc.g_top, fx * h + 0.5 * h);
        else                      v = 2.0 * tf * EvalBC(bc.g_top, fx * h + 0.5 * h) / h;
        PetscCall(VecSetValue(b, row, v, INSERT_VALUES));
    }}
  }

  // u RHS: u_old + dt*f
  for (PetscInt ey = 0; ey < Ny; ++ey) {
    for (PetscInt ex = 0; ex < Nx; ++ex) {
      PetscInt row = idx_u(ex, ey, Nx, Ny);
      if (row < rstart || row >= rend) continue;
      PetscScalar uold; PetscCall(VecGetValues(xOld, 1, &row, &uold));
      PetscScalar val = uold + dt * mms.f_td((ex + 0.5) * h, (ey + 0.5) * h, t_eval, alpha);
      PetscCall(VecSetValue(b, row, val, INSERT_VALUES));
    }
  }

  PetscCall(VecAssemblyBegin(b));
  PetscCall(VecAssemblyEnd(b));
  PetscFunctionReturn(0);
}

// ============================================================================
// Compute L2 error for u (element DOFs) against time-dependent MMS
// ============================================================================
static PetscErrorCode ComputeError(Vec x, PetscInt Nx, PetscInt Ny,
                                     PetscReal t, const ManufacturedSolution &mms,
                                     PetscReal *error) {
  PetscFunctionBeginUser;
  const PetscReal hx = 1.0 / Nx;
  const PetscReal hy = 1.0 / Ny;

  PetscReal localSum = 0.0;
  PetscInt rstart, rend;
  PetscCall(VecGetOwnershipRange(x, &rstart, &rend));

  for (PetscInt ey = 0; ey < Ny; ++ey) {
    for (PetscInt ex = 0; ex < Nx; ++ex) {
      PetscInt row = idx_u(ex, ey, Nx, Ny);
      if (row < rstart || row >= rend) continue;
      PetscScalar val;
      PetscCall(VecGetValues(x, 1, &row, &val));
      const PetscScalar xc = (ex + 0.5) * hx;
      const PetscScalar yc = (ey + 0.5) * hy;
      const PetscScalar diff = val - mms.u_td(xc, yc, t, alpha);
      localSum += PetscRealPart(diff * diff);
    }
  }

  PetscReal globalSum;
  MPI_Allreduce(&localSum, &globalSum, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
  *error = std::sqrt(globalSum * hx * hy);
  PetscFunctionReturn(0);
}

// ============================================================================
// Main
// ============================================================================
int main(int argc, char **argv) {
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL,
                            "Staggered (mixed) heat equation on DMStag"));

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscInt Nx = 64, Ny = 64;
  PetscReal T_final = 0.05;
  PetscBool compute_error = PETSC_FALSE;
  PetscBool convergence_test = PETSC_FALSE;
  char mms_name[32]="sinpi"; PetscBool flg;
  PetscCall(PetscOptionsGetString(NULL,NULL,"-mms",mms_name,sizeof(mms_name),&flg));
  const ManufacturedSolution *mms=&MMS_SINPI;
  if(strcmp(mms_name,"poly2")==0) mms=&MMS_POLY2;
  else if(strcmp(mms_name,"cospi")==0) mms=&MMS_COSPI;

  char bc_str[8]=""; PetscCall(PetscOptionsGetString(NULL,NULL,"-bc_type",bc_str,sizeof(bc_str),&flg));
  BCType bl=BC_DIRICHLET, br=BC_DIRICHLET, bb=BC_DIRICHLET, bt=BC_DIRICHLET;
  if(strlen(bc_str)==4){
    bl=(bc_str[0]=='N')?BC_NEUMANN:BC_DIRICHLET; br=(bc_str[1]=='N')?BC_NEUMANN:BC_DIRICHLET;
    bb=(bc_str[2]=='N')?BC_NEUMANN:BC_DIRICHLET; bt=(bc_str[3]=='N')?BC_NEUMANN:BC_DIRICHLET;
  } else if(strcmp(mms_name,"cospi")==0){
    bl=br=bb=bt=BC_NEUMANN;
  }
  BoundaryCondition bc = MakeBCFromMMS(*mms, bl, br, bb, bt);

  PetscCall(PetscOptionsGetInt(NULL, NULL, "-nx", &Nx, NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-ny", &Ny, NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-T", &T_final, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-heat_check_error", &compute_error, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-convergence_test", &convergence_test, NULL));

  // Time step: dt ∝ h² for convergence test
  const PetscReal h = 1.0 / Nx;
  const PetscReal dt = h * h;
  const PetscInt Nsteps = static_cast<PetscInt>(std::ceil(T_final / dt));
  const PetscReal dtActual = T_final / Nsteps;

  PetscInt Ntotal = totalDOFs(Nx, Ny);

  // ---- Create monolithic matrix ----
  Mat A;
  PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
  PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, Ntotal, Ntotal));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatSetUp(A));
  PetscCall(MatSeqAIJSetPreallocation(A, 5, NULL));
  PetscCall(MatMPIAIJSetPreallocation(A, 5, NULL, 5, NULL));

  PetscCall(AssembleSystem(A, Nx, Ny, dtActual, bc));

  Vec x, xOld, b;
  PetscCall(MatCreateVecs(A, &x, &b));
  PetscCall(VecDuplicate(x, &xOld));

  // Initial condition
  PetscCall(SetInitialCondition(x, Nx, Ny, *mms));

  // ---- KSP solver ----
  KSP ksp;
  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
  PetscCall(KSPSetOperators(ksp, A, A));
  PetscCall(KSPSetTolerances(ksp, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT, 5000));
  { PC pc; KSPGetPC(ksp, &pc); PCSetType(pc, PCLU); }  // saddle-point needs direct solve
  PetscCall(KSPSetFromOptions(ksp));

  if (rank == 0) {
    std::cout << "Staggered heat: N=" << Nx << "x" << Ny
              << ", h=" << h << ", dt=" << dtActual
              << ", steps=" << Nsteps << ", DOFs=" << Ntotal
              << ", mms=" << mms->name << std::endl;
  }

  // ---- Time stepping ----
  PetscReal t = 0.0;
  for (PetscInt step = 1; step <= Nsteps; ++step) {
    t += dtActual;
    PetscCall(VecCopy(x, xOld));
    PetscCall(BuildRHS(b, xOld, Nx, Ny, t - dtActual, dtActual, *mms, bc));
    PetscCall(KSPSolve(ksp, b, x));

    if (rank == 0 && (step % 100 == 0 || step == Nsteps)) {
      PetscInt its;
      PetscCall(KSPGetIterationNumber(ksp, &its));
      std::cout << "  step " << step << "/" << Nsteps
                << ", KSP its=" << its << std::endl;
    }
  }

  // ---- Compute error ----
  PetscReal error = 0.0;
  if (compute_error) {
    PetscCall(ComputeError(x, Nx, Ny, t, *mms, &error));
    if (rank == 0) {
      std::cout << "||u(T) - u_exact(T)||_L2 = " << error << std::endl;
    }
  }

  if (convergence_test && rank == 0) {
    std::cout << "CONVERGENCE: " << Nx << " " << h << " " << error
              << " " << Nsteps << std::endl;
  }

  // ---- Cleanup ----
  PetscCall(KSPDestroy(&ksp));
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&xOld));
  PetscCall(VecDestroy(&b));
  PetscCall(MatDestroy(&A));
  PetscCall(PetscFinalize());
  return 0;
}
