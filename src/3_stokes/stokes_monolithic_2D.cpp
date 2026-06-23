#include <petsc.h>
#include <cmath>
#include <iostream>
#include "analytical/unsteady.h"
using namespace UNSTEADY::SIN2_STOKES_2D;

static const PetscScalar nu = 0.1;

static inline PetscInt iu(PetscInt fx, PetscInt fy, PetscInt Nx, PetscInt Ny) {
  return fy*(Nx+1) + fx;
}
static inline PetscInt iv(PetscInt fx, PetscInt fy, PetscInt Nx, PetscInt Ny) {
  return (Nx+1)*Ny + fy*Nx + fx;
}
static inline PetscInt ip(PetscInt ex, PetscInt ey, PetscInt Nx, PetscInt Ny) {
  return (Nx+1)*Ny + Nx*(Ny+1) + ey*Nx + ex;
}
static inline PetscInt Ndof(PetscInt Nx, PetscInt Ny) {
  return (Nx+1)*Ny + Nx*(Ny+1) + Nx*Ny;
}

int main(int argc, char **argv) {
  PetscCall(PetscInitialize(&argc, &argv, NULL, "Monolithic Stokes + pin"));
  int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscInt Nx = 32, Ny = 32; PetscReal Tf = 0.01;
  PetscBool chk = PETSC_FALSE, ct = PETSC_FALSE;
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-nx", &Nx, NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-ny", &Ny, NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-T", &Tf, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-stokes_check_error", &chk, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-convergence_test", &ct, NULL));

  const PetscReal h = 1.0 / Nx, dt = h * h;
  const PetscInt Ns = (PetscInt)ceil(Tf / dt);
  const PetscReal dta = Tf / Ns;
  const PetscInt Ntot = Ndof(Nx, Ny);

  // -- Assemble matrix --
  Mat A; MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, Ntot, Ntot);
  MatSetUp(A); MatSeqAIJSetPreallocation(A, 10, NULL);
  MatMPIAIJSetPreallocation(A, 10, NULL, 10, NULL);
  {
    const PetscReal ih = 1.0 / h, ih2 = 1.0 / (h * h);
    const PetscReal c = nu * dta, ch2 = c * ih2;
    PetscInt rs, re; MatGetOwnershipRange(A, &rs, &re);

    for (PetscInt fy = 0; fy < Ny; ++fy) for (PetscInt fx = 0; fx <= Nx; ++fx) {
      PetscInt row = iu(fx, fy, Nx, Ny); if (row < rs || row >= re) continue;
      if (fx == 0 || fx == Nx) { MatSetValue(A, row, row, 1.0, INSERT_VALUES); continue; }
      MatSetValue(A, row, row, 1.0 + ((fy == 0 || fy == Ny - 1) ? 5 : 4) * ch2, INSERT_VALUES);
      if (fx > 0)     MatSetValue(A, row, iu(fx - 1, fy, Nx, Ny), -ch2, INSERT_VALUES);
      if (fx < Nx)    MatSetValue(A, row, iu(fx + 1, fy, Nx, Ny), -ch2, INSERT_VALUES);
      if (fy > 0)     MatSetValue(A, row, iu(fx, fy - 1, Nx, Ny), -ch2, INSERT_VALUES);
      if (fy < Ny - 1) MatSetValue(A, row, iu(fx, fy + 1, Nx, Ny), -ch2, INSERT_VALUES);
      if (fx < Nx && fy < Ny) MatSetValue(A, row, ip(fx, fy, Nx, Ny),    dta * ih, INSERT_VALUES);
      if (fx > 0 && fy < Ny)  MatSetValue(A, row, ip(fx - 1, fy, Nx, Ny), -dta * ih, INSERT_VALUES);
    }
    for (PetscInt fy = 0; fy <= Ny; ++fy) for (PetscInt fx = 0; fx < Nx; ++fx) {
      PetscInt row = iv(fx, fy, Nx, Ny); if (row < rs || row >= re) continue;
      if (fy == 0 || fy == Ny) { MatSetValue(A, row, row, 1.0, INSERT_VALUES); continue; }
      MatSetValue(A, row, row, 1.0 + ((fx == 0 || fx == Nx - 1) ? 5 : 4) * ch2, INSERT_VALUES);
      if (fx > 0)     MatSetValue(A, row, iv(fx - 1, fy, Nx, Ny), -ch2, INSERT_VALUES);
      if (fx < Nx - 1) MatSetValue(A, row, iv(fx + 1, fy, Nx, Ny), -ch2, INSERT_VALUES);
      if (fy > 0)     MatSetValue(A, row, iv(fx, fy - 1, Nx, Ny), -ch2, INSERT_VALUES);
      if (fy < Ny)    MatSetValue(A, row, iv(fx, fy + 1, Nx, Ny), -ch2, INSERT_VALUES);
      if (fy < Ny && fx < Nx) MatSetValue(A, row, ip(fx, fy, Nx, Ny),    dta * ih, INSERT_VALUES);
      if (fy > 0 && fx < Nx)  MatSetValue(A, row, ip(fx, fy - 1, Nx, Ny), -dta * ih, INSERT_VALUES);
    }
    // p rows: (u_R-u_L)/h + (v_T-v_B)/h = 0, pin p(0,0)
    for (PetscInt ey = 0; ey < Ny; ++ey) for (PetscInt ex = 0; ex < Nx; ++ex) {
      PetscInt row = ip(ex, ey, Nx, Ny); if (row < rs || row >= re) continue;
      if (ex == 0 && ey == 0) {
        MatSetValue(A, row, row, 1.0, INSERT_VALUES); // pin
      } else {
        MatSetValue(A, row, row, 0.0, INSERT_VALUES);
        MatSetValue(A, row, iu(ex + 1, ey, Nx, Ny),  ih, INSERT_VALUES);
        MatSetValue(A, row, iu(ex, ey, Nx, Ny),     -ih, INSERT_VALUES);
        MatSetValue(A, row, iv(ex, ey + 1, Nx, Ny),  ih, INSERT_VALUES);
        MatSetValue(A, row, iv(ex, ey, Nx, Ny),     -ih, INSERT_VALUES);
      }
    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  }

  Vec x, b; MatCreateVecs(A, &x, &b);

  // Init
  { PetscInt rs, re; VecGetOwnershipRange(x, &rs, &re); VecSet(x, 0);
    for (PetscInt fy = 0; fy < Ny; ++fy) for (PetscInt fx = 0; fx <= Nx; ++fx) {
      PetscInt row = iu(fx, fy, Nx, Ny); if (row < rs || row >= re) continue;
      VecSetValue(x, row, u_exact(fx * h, (fy + 0.5) * h, 0, nu), INSERT_VALUES);
    }
    for (PetscInt fy = 0; fy <= Ny; ++fy) for (PetscInt fx = 0; fx < Nx; ++fx) {
      PetscInt row = iv(fx, fy, Nx, Ny); if (row < rs || row >= re) continue;
      VecSetValue(x, row, v_exact((fx + 0.5) * h, fy * h, 0, nu), INSERT_VALUES);
    }
    VecAssemblyBegin(x); VecAssemblyEnd(x);
  }

  KSP ksp; KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, A, A);
  KSPSetType(ksp, KSPPREONLY);
  { PC pc; KSPGetPC(ksp, &pc); PCSetType(pc, PCLU); }
  PetscOptionsInsertString(NULL, "-pc_factor_shift_type nonzero");
  KSPSetFromOptions(ksp);
  if (rank == 0) printf("Monolithic Stokes (pin): N=%dx%d h=%g dt=%g steps=%d DOFs=%d\n",
    Nx, Ny, h, dta, Ns, Ntot);

  PetscReal t = 0;
  for (PetscInt s = 1; s <= Ns; ++s) { t += dta;
    // RHS = x_old + Δt * f
    PetscCall(VecCopy(x, b));
    { PetscInt rs, re; VecGetOwnershipRange(b, &rs, &re);
      for (PetscInt fy = 0; fy < Ny; ++fy) for (PetscInt fx = 1; fx < Nx; ++fx) {
        PetscInt row = iu(fx, fy, Nx, Ny); if (row < rs || row >= re) continue;
        PetscScalar v; VecGetValues(b, 1, &row, &v);
        VecSetValue(b, row, v + dta * fx_stokes(fx * h, (fy + 0.5) * h, t, nu), INSERT_VALUES);
      }
      for (PetscInt fy = 1; fy < Ny; ++fy) for (PetscInt fx = 0; fx < Nx; ++fx) {
        PetscInt row = iv(fx, fy, Nx, Ny); if (row < rs || row >= re) continue;
        PetscScalar v; VecGetValues(b, 1, &row, &v);
        VecSetValue(b, row, v + dta * fy_stokes((fx + 0.5) * h, fy * h, t, nu), INSERT_VALUES);
      }
      // Dirichlet BC for u at x=0, x=1 and v at y=0, y=1
      for (PetscInt fy = 0; fy < Ny; ++fy) {
        PetscInt row = iu(0, fy, Nx, Ny);  if (row >= rs && row < re) VecSetValue(b, row, u_exact(0.0, (fy+0.5)*h, t, nu), INSERT_VALUES);
        row = iu(Nx, fy, Nx, Ny); if (row >= rs && row < re) VecSetValue(b, row, u_exact(1.0, (fy+0.5)*h, t, nu), INSERT_VALUES);
      }
      for (PetscInt fx = 0; fx < Nx; ++fx) {
        PetscInt row = iv(fx, 0, Nx, Ny);  if (row >= rs && row < re) VecSetValue(b, row, v_exact((fx+0.5)*h, 0.0, t, nu), INSERT_VALUES);
        row = iv(fx, Ny, Nx, Ny); if (row >= rs && row < re) VecSetValue(b, row, v_exact((fx+0.5)*h, 1.0, t, nu), INSERT_VALUES);
      }
      // Clear p and pin p(0,0)
      for (PetscInt ey = 0; ey < Ny; ++ey) for (PetscInt ex = 0; ex < Nx; ++ex) {
        PetscInt row = ip(ex, ey, Nx, Ny); if (row < rs || row >= re) continue;
        if (ex == 0 && ey == 0) VecSetValue(b, row, p_exact(0.5*h, 0.5*h, t, nu), INSERT_VALUES);
        else VecSetValue(b, row, 0.0, INSERT_VALUES);
      }
      VecAssemblyBegin(b); VecAssemblyEnd(b);
    }
    PetscCall(KSPSolve(ksp, b, x));
  }
  if (rank == 0) printf("  done %d steps\n", Ns);

  PetscReal eu = 0, ev = 0, ep = 0;
  if (chk) {
    double su = 0, sv = 0, sp = 0; PetscInt rs, re; VecGetOwnershipRange(x, &rs, &re);
    for (PetscInt fy = 0; fy < Ny; ++fy) for (PetscInt fx = 1; fx < Nx; ++fx) {
      PetscInt row = iu(fx, fy, Nx, Ny); if (row < rs || row >= re) continue;
      PetscScalar val; VecGetValues(x, 1, &row, &val);
      double d = val - u_exact(fx * h, (fy + 0.5) * h, t, nu); su += d * d;
    }
    for (PetscInt fy = 1; fy < Ny; ++fy) for (PetscInt fx = 0; fx < Nx; ++fx) {
      PetscInt row = iv(fx, fy, Nx, Ny); if (row < rs || row >= re) continue;
      PetscScalar val; VecGetValues(x, 1, &row, &val);
      double d = val - v_exact((fx + 0.5) * h, fy * h, t, nu); sv += d * d;
    }
    for (PetscInt ey = 0; ey < Ny; ++ey) for (PetscInt ex = 0; ex < Nx; ++ex) {
      PetscInt row = ip(ex, ey, Nx, Ny); if (row < rs || row >= re) continue;
      PetscScalar val; VecGetValues(x, 1, &row, &val);
      double d = val - p_exact((ex + 0.5) * h, (ey + 0.5) * h, t, nu); sp += d * d;
    }
    double g; MPI_Allreduce(&su, &g, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); eu = sqrt(g * h * h);
    MPI_Allreduce(&sv, &g, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); ev = sqrt(g * h * h);
    MPI_Allreduce(&sp, &g, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); ep = sqrt(g * h * h);
    if (rank == 0) printf("||u-u_ex||=%g ||v-v_ex||=%g ||p-p_ex||=%g\n", eu, ev, ep);
  }
  if (ct && rank == 0) printf("CONVERGENCE: %d %g %g %g %g %d\n", Nx, h, eu, ev, ep, Ns);

  KSPDestroy(&ksp); VecDestroy(&x); VecDestroy(&b); MatDestroy(&A); PetscFinalize(); return 0;
}
