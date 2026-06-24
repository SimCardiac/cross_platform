#include <petsc.h>
#include <cmath>
#include <iostream>
#include "analytical/unsteady.h"
using namespace UNSTEADY::SINPI_SCALAR_2D;

// Standard interface — staggered (mixed) Poisson
static inline PetscInt idx_u(PetscInt ex,PetscInt ey,PetscInt Nx,PetscInt Ny){ return ey*Nx+ex; }
static inline PetscInt idx_gx(PetscInt fx,PetscInt fy,PetscInt Nx,PetscInt Ny){ return Nx*Ny+fy*(Nx+1)+fx; }
static inline PetscInt idx_gy(PetscInt fx,PetscInt fy,PetscInt Nx,PetscInt Ny){ return Nx*Ny+(Nx+1)*Ny+fx*(Ny+1)+fy; }
static inline PetscInt Nd(PetscInt Nx,PetscInt Ny){ return Nx*Ny+(Nx+1)*Ny+Nx*(Ny+1); }

PetscErrorCode SetInitialCondition(DM dm, Vec u, PetscReal t) {
  PetscFunctionBeginUser; (void)dm;(void)t; VecSet(u,0); PetscFunctionReturn(0);
}

PetscErrorCode AssembleSystem(DM dm, Mat A, PetscReal dt) {
  PetscFunctionBeginUser; (void)dm;(void)dt;
  PetscInt Nx,Ny; DMStagGetGlobalSizes(dm,&Nx,&Ny,NULL);
  PetscReal h=1.0/Nx, ih=1.0/h, ih2=2.0/h, ix2=1.0/(h*h), iy2=ix2;
  PetscInt rs,re; MatGetOwnershipRange(A,&rs,&re);
  for(PetscInt fy=0;fy<Ny;fy++) for(PetscInt fx=0;fx<=Nx;fx++){
    PetscInt row=idx_gx(fx,fy,Nx,Ny); if(row<rs||row>=re) continue;
    MatSetValue(A,row,row,1.0,INSERT_VALUES);
    if(fx==0){ MatSetValue(A,row,idx_u(0,fy,Nx,Ny),-ih2,INSERT_VALUES); }
    else if(fx==Nx){ MatSetValue(A,row,idx_u(Nx-1,fy,Nx,Ny),ih2,INSERT_VALUES); }
    else{ MatSetValue(A,row,idx_u(fx,fy,Nx,Ny),-ih,INSERT_VALUES); MatSetValue(A,row,idx_u(fx-1,fy,Nx,Ny),ih,INSERT_VALUES); }
  }
  for(PetscInt fx=0;fx<Nx;fx++) for(PetscInt fy=0;fy<=Ny;fy++){
    PetscInt row=idx_gy(fx,fy,Nx,Ny); if(row<rs||row>=re) continue;
    MatSetValue(A,row,row,1.0,INSERT_VALUES);
    if(fy==0){ MatSetValue(A,row,idx_u(fx,0,Nx,Ny),-ih2,INSERT_VALUES); }
    else if(fy==Ny){ MatSetValue(A,row,idx_u(fx,Ny-1,Nx,Ny),ih2,INSERT_VALUES); }
    else{ MatSetValue(A,row,idx_u(fx,fy,Nx,Ny),-ih,INSERT_VALUES); MatSetValue(A,row,idx_u(fx,fy-1,Nx,Ny),ih,INSERT_VALUES); }
  }
  for(PetscInt ey=0;ey<Ny;ey++) for(PetscInt ex=0;ex<Nx;ex++){
    PetscInt row=idx_u(ex,ey,Nx,Ny); if(row<rs||row>=re) continue;
    MatSetValue(A,row,idx_gx(ex+1,ey,Nx,Ny),-ih,INSERT_VALUES);
    MatSetValue(A,row,idx_gx(ex,ey,Nx,Ny),ih,INSERT_VALUES);
    MatSetValue(A,row,idx_gy(ex,ey+1,Nx,Ny),-ih,INSERT_VALUES);
    MatSetValue(A,row,idx_gy(ex,ey,Nx,Ny),ih,INSERT_VALUES);
  }
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  PetscFunctionReturn(0);
}

PetscErrorCode BuildRHS(DM dm, Vec b, const Vec uOld, PetscReal t, PetscReal dt) {
  PetscFunctionBeginUser; (void)uOld;(void)t;(void)dt;
  PetscInt Nx,Ny; DMStagGetGlobalSizes(dm,&Nx,&Ny,NULL);
  PetscReal h=1.0/Nx; VecSet(b,0);
  PetscInt rs,re; VecGetOwnershipRange(b,&rs,&re);
  for(PetscInt ey=0;ey<Ny;ey++) for(PetscInt ex=0;ex<Nx;ex++){
    PetscInt row=idx_u(ex,ey,Nx,Ny); if(row<rs||row>=re) continue;
    VecSetValue(b,row,f_poisson((ex+0.5)*h,(ey+0.5)*h),INSERT_VALUES);
  }
  VecAssemblyBegin(b); VecAssemblyEnd(b);
  PetscFunctionReturn(0);
}

PetscErrorCode ComputeError(DM dm, const Vec x, PetscReal t, PetscReal *err) {
  PetscFunctionBeginUser; (void)dm;(void)t;
  PetscInt Nx,Ny; DMStagGetGlobalSizes(dm,&Nx,&Ny,NULL);
  PetscReal h=1.0/Nx; double s=0;
  PetscInt rs,re; VecGetOwnershipRange(x,&rs,&re);
  for(PetscInt ey=0;ey<Ny;ey++) for(PetscInt ex=0;ex<Nx;ex++){
    PetscInt row=idx_u(ex,ey,Nx,Ny); if(row<rs||row>=re) continue;
    PetscScalar v; VecGetValues(x,1,&row,&v);
    double d=v-u_steady((ex+0.5)*h,(ey+0.5)*h); s+=d*d;
  }
  double g; MPI_Allreduce(&s,&g,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  *err=sqrt(g*h*h); PetscFunctionReturn(0);
}

int main(int argc,char**argv){
  PetscCall(PetscInitialize(&argc,&argv,NULL,"Staggered Poisson"));
  int rank; MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  PetscInt Nx=64,Ny=64; PetscBool chk=PETSC_FALSE,ct=PETSC_FALSE;
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-nx",&Nx,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-ny",&Ny,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-poisson_check_error",&chk,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-convergence_test",&ct,NULL));
  DM dm; DMStagCreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,Nx,Ny,PETSC_DECIDE,PETSC_DECIDE,0,0,0,DMSTAG_STENCIL_BOX,1,NULL,NULL,&dm);
  DMSetUp(dm);
  PetscInt Ntot=Nd(Nx,Ny);
  Mat A; MatCreate(PETSC_COMM_WORLD,&A); MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,Ntot,Ntot);
  MatSetUp(A); MatSeqAIJSetPreallocation(A,5,NULL); MatMPIAIJSetPreallocation(A,5,NULL,5,NULL);
  AssembleSystem(dm,A,0);
  Vec x,b; MatCreateVecs(A,&x,&b); BuildRHS(dm,b,NULL,0,0);
  KSP ksp; KSPCreate(PETSC_COMM_WORLD,&ksp); KSPSetOperators(ksp,A,A); KSPSetType(ksp,KSPGMRES);
  KSPSetTolerances(ksp,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,10000);
  {PC pc;KSPGetPC(ksp,&pc);PCSetType(pc,PCNONE);} KSPSetFromOptions(ksp);
  if(rank==0) printf("Stagg Poisson: N=%dx%d h=%g DOFs=%d\n",Nx,Ny,1.0/Nx,Ntot);
  KSPSolve(ksp,b,x);
  PetscReal err=0; if(chk){ ComputeError(dm,x,0,&err); if(rank==0) printf("||u-u_ex||=%g\n",err); }
  if(ct&&rank==0) printf("CONVERGENCE: %d %g %g 0\n",Nx,1.0/Nx,err);
  KSPDestroy(&ksp); VecDestroy(&x); VecDestroy(&b); MatDestroy(&A); DMDestroy(&dm);
  PetscFinalize(); return 0;
}
