#include <petsc.h>
#include <cmath>
#include <iostream>
#include "analytical/unsteady.h"
using namespace UNSTEADY::SINPI_SCALAR_2D;

static inline PetscInt idx(PetscInt ex,PetscInt ey,PetscInt Nx,PetscInt Ny){ (void)Ny; return ey*Nx+ex; }

PetscErrorCode SetInitialCondition(DM dm, Vec u, PetscReal t) {
  PetscFunctionBeginUser; (void)dm;(void)t; VecSet(u,0); PetscFunctionReturn(0);
}

PetscErrorCode AssembleSystem(DM dm, Mat A, PetscReal dt) {
  PetscFunctionBeginUser; (void)dm;(void)dt;
  PetscInt Nx,Ny; DMStagGetGlobalSizes(dm,&Nx,&Ny,NULL);
  PetscReal hx=1.0/Nx, hy=1.0/Ny, ix2=1.0/(hx*hx), iy2=1.0/(hy*hy), diag=2.0*(ix2+iy2);
  for(PetscInt ey=0;ey<Ny;ey++) for(PetscInt ex=0;ex<Nx;ex++){
    PetscInt row=idx(ex,ey,Nx,Ny);
    PetscReal d=diag;
    if(ex==0) d+=ix2; if(ex==Nx-1) d+=ix2;  // ghost cell: u_ghost=-u_interior
    if(ey==0) d+=iy2; if(ey==Ny-1) d+=iy2;
    MatSetValue(A,row,row,d,INSERT_VALUES);
    if(ex>0)    MatSetValue(A,row,idx(ex-1,ey,Nx,Ny),-ix2,INSERT_VALUES);
    if(ex<Nx-1) MatSetValue(A,row,idx(ex+1,ey,Nx,Ny),-ix2,INSERT_VALUES);
    if(ey>0)    MatSetValue(A,row,idx(ex,ey-1,Nx,Ny),-iy2,INSERT_VALUES);
    if(ey<Ny-1) MatSetValue(A,row,idx(ex,ey+1,Nx,Ny),-iy2,INSERT_VALUES);
  }
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  PetscFunctionReturn(0);
}

PetscErrorCode BuildRHS(DM dm, Vec b, const Vec uOld, PetscReal t, PetscReal dt) {
  PetscFunctionBeginUser; (void)uOld;(void)t;(void)dt;
  PetscInt Nx,Ny; DMStagGetGlobalSizes(dm,&Nx,&Ny,NULL);
  PetscReal hx=1.0/Nx, hy=1.0/Ny; VecSet(b,0);
  for(PetscInt ey=0;ey<Ny;ey++) for(PetscInt ex=0;ex<Nx;ex++){
    PetscScalar x=(ex+0.5)*hx, y=(ey+0.5)*hy;
    VecSetValue(b,idx(ex,ey,Nx,Ny),f_poisson(x,y),INSERT_VALUES);
  }
  VecAssemblyBegin(b); VecAssemblyEnd(b);
  PetscFunctionReturn(0);
}

PetscErrorCode ComputeError(DM dm, const Vec u, PetscReal t, PetscReal *error) {
  PetscFunctionBeginUser; (void)dm;(void)t;
  PetscInt Nx,Ny; DMStagGetGlobalSizes(dm,&Nx,&Ny,NULL);
  PetscReal hx=1.0/Nx, hy=1.0/Ny; double s=0;
  for(PetscInt ey=0;ey<Ny;ey++) for(PetscInt ex=0;ex<Nx;ex++){
    PetscScalar v; PetscInt row=idx(ex,ey,Nx,Ny);
    VecGetValues(u,1,&row,&v);
    double d=v-u_steady((ex+0.5)*hx,(ey+0.5)*hy); s+=d*d;
  }
  double g; MPI_Allreduce(&s,&g,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  *error=sqrt(g*hx*hy); PetscFunctionReturn(0);
}

int main(int argc,char**argv){
  PetscCall(PetscInitialize(&argc,&argv,NULL,"Cell Poisson"));
  int rank; MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  PetscInt Nx=64,Ny=64; PetscBool chk=PETSC_FALSE,ct=PETSC_FALSE;
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-nx",&Nx,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-ny",&Ny,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-poisson_check_error",&chk,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-convergence_test",&ct,NULL));
  PetscInt Ntot=Nx*Ny;
  Mat A; MatCreate(PETSC_COMM_WORLD,&A); MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,Ntot,Ntot);
  MatSetUp(A); MatSeqAIJSetPreallocation(A,5,NULL); MatMPIAIJSetPreallocation(A,5,NULL,5,NULL);
  DM dm; DMStagCreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,Nx,Ny,PETSC_DECIDE,PETSC_DECIDE,1,0,0,DMSTAG_STENCIL_BOX,1,NULL,NULL,&dm);
  DMSetUp(dm); DMStagSetUniformCoordinatesProduct(dm,0,1,0,1,0,0);
  AssembleSystem(dm,A,0);
  Vec u,b; MatCreateVecs(A,&u,&b); BuildRHS(dm,b,NULL,0,0);
  KSP ksp; KSPCreate(PETSC_COMM_WORLD,&ksp); KSPSetOperators(ksp,A,A);
  KSPSetTolerances(ksp,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,10000); KSPSetFromOptions(ksp);
  if(rank==0) printf("Cell Poisson: N=%dx%d h=%g\n",Nx,Ny,1.0/Nx);
  KSPSolve(ksp,b,u);
  PetscReal err=0; if(chk){ ComputeError(dm,u,0,&err); if(rank==0) printf("||u-u_ex||=%g\n",err); }
  if(ct&&rank==0) printf("CONVERGENCE: %d %g %g 0\n",Nx,1.0/Nx,err);
  KSPDestroy(&ksp); VecDestroy(&u); VecDestroy(&b); MatDestroy(&A); DMDestroy(&dm);
  PetscFinalize(); return 0;
}
