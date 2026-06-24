#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <cmath>
#include <iostream>
#include "analytical/unsteady.h"
using namespace UNSTEADY::SINPI_SCALAR_2D;

// ============================================================================
// Standardized solver interface (steady: t=0, dt=0)
// ============================================================================

// --- Set initial/exact condition at time t ---
PetscErrorCode SetInitialCondition(DM dm, Vec u, PetscReal t) {
  PetscFunctionBeginUser;
  (void)t;
  PetscInt Nx,Ny; DMDAGetInfo(dm,NULL,&Nx,&Ny,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
  PetscReal hx=1.0/(Nx-1), hy=1.0/(Ny-1);
  VecSet(u,0);
  PetscInt rs,re; VecGetOwnershipRange(u,&rs,&re);
  for(PetscInt j=0;j<Ny;j++) for(PetscInt i=0;i<Nx;i++){
    PetscInt row=j*Nx+i; if(row<rs||row>=re) continue;
    VecSetValue(u,row,u_steady(i*hx,j*hy),INSERT_VALUES);
  }
  VecAssemblyBegin(u); VecAssemblyEnd(u);
  PetscFunctionReturn(0);
}

// --- Assemble matrix (BC in boundary rows, 5-point Laplacian interior) ---
PetscErrorCode AssembleSystem(DM dm, Mat A, PetscReal dt) {
  PetscFunctionBeginUser;
  (void)dt;
  PetscInt Nx, Ny; DMDAGetInfo(dm, NULL, &Nx, &Ny, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  PetscReal hx=1.0/(Nx-1), hy=1.0/(Ny-1), ix2=1.0/(hx*hx), iy2=1.0/(hy*hy), diag=2.0*(ix2+iy2);
  for(PetscInt j=0;j<Ny;j++) for(PetscInt i=0;i<Nx;i++){
    PetscInt row=j*Nx+i;
    if(i==0||i==Nx-1||j==0||j==Ny-1){
      MatSetValue(A,row,row,1.0,INSERT_VALUES);
    } else {
      MatSetValue(A,row,row,diag,INSERT_VALUES);
      MatSetValue(A,row,j*Nx+(i-1),-ix2,INSERT_VALUES);
      MatSetValue(A,row,j*Nx+(i+1),-ix2,INSERT_VALUES);
      MatSetValue(A,row,(j-1)*Nx+i,-iy2,INSERT_VALUES);
      MatSetValue(A,row,(j+1)*Nx+i,-iy2,INSERT_VALUES);
    }
  }
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  PetscFunctionReturn(0);
}

// --- Build RHS (t and uOld ignored for steady) ---
PetscErrorCode BuildRHS(DM dm, Vec b, const Vec uOld, PetscReal t, PetscReal dt) {
  PetscFunctionBeginUser;
  (void)uOld;(void)t;(void)dt;
  PetscInt Nx,Ny; DMDAGetInfo(dm,NULL,&Nx,&Ny,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
  PetscReal hx=1.0/(Nx-1);
  VecSet(b,0);
  PetscInt rs,re; VecGetOwnershipRange(b,&rs,&re);
  for(PetscInt j=0;j<Ny;j++) for(PetscInt i=0;i<Nx;i++){
    PetscInt row=j*Nx+i; if(row<rs||row>=re) continue;
    if(!(i==0||i==Nx-1||j==0||j==Ny-1))
      VecSetValue(b,row,f_poisson(i*hx,j*hx),INSERT_VALUES);
  }
  VecAssemblyBegin(b); VecAssemblyEnd(b);
  PetscFunctionReturn(0);
}

// --- Compute L2 error at time t ---
PetscErrorCode ComputeError(DM dm, const Vec u, PetscReal t, PetscReal *error) {
  PetscFunctionBeginUser;
  (void)t;
  PetscInt Nx,Ny; DMDAGetInfo(dm,NULL,&Nx,&Ny,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
  PetscReal hx=1.0/(Nx-1), hy=1.0/(Ny-1);
  Vec ue; VecDuplicate(u,&ue);
  PetscInt rs,re; VecGetOwnershipRange(u,&rs,&re);
  for(PetscInt j=0;j<Ny;j++) for(PetscInt i=0;i<Nx;i++){
    PetscInt row=j*Nx+i; if(row<rs||row>=re) continue;
    VecSetValue(ue,row,u_steady(i*hx,j*hy),INSERT_VALUES);
  }
  VecAssemblyBegin(ue); VecAssemblyEnd(ue);
  Vec d; VecDuplicate(u,&d); VecCopy(u,d); VecAXPY(d,-1.0,ue);
  PetscReal n2; VecNorm(d,NORM_2,&n2);
  *error = n2*sqrt(hx*hy);
  VecDestroy(&d); VecDestroy(&ue);
  PetscFunctionReturn(0);
}

int main(int argc,char**argv){
  PetscCall(PetscInitialize(&argc,&argv,NULL,"Vertex-centered Poisson"));
  int rank; MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  PetscInt Nx=64,Ny=64; PetscBool chk=PETSC_FALSE,ct=PETSC_FALSE;
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-nx",&Nx,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-ny",&Ny,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-poisson_check_error",&chk,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-convergence_test",&ct,NULL));
  DM dm; DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,Nx,Ny,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,&dm);
  DMSetUp(dm); DMDASetUniformCoordinates(dm,0,1,0,1,0,0);
  Mat A; DMCreateMatrix(dm,&A);
  MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
  Vec u,b; DMCreateGlobalVector(dm,&u); DMCreateGlobalVector(dm,&b);
  AssembleSystem(dm,A,0);
  BuildRHS(dm,b,NULL,0,0);
  KSP ksp; KSPCreate(PETSC_COMM_WORLD,&ksp); KSPSetOperators(ksp,A,A);
  KSPSetTolerances(ksp,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,10000);
  KSPSetFromOptions(ksp);
  if(rank==0) printf("Vertex Poisson: N=%dx%d h=%g\n",Nx,Ny,1.0/(Nx-1));
  KSPSolve(ksp,b,u);
  PetscReal err=0; if(chk){ ComputeError(dm,u,0,&err); if(rank==0) printf("||u-u_ex||=%g\n",err); }
  if(ct&&rank==0) printf("CONVERGENCE: %d %g %g 0\n",Nx,1.0/(Nx-1),err);
  KSPDestroy(&ksp); VecDestroy(&u); VecDestroy(&b); MatDestroy(&A); DMDestroy(&dm);
  PetscFinalize(); return 0;
}
