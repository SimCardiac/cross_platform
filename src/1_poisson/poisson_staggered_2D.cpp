#include <petsc.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include "analytical/unsteady.h"
#include "common/boundary.h"
#include "common/mms.h"

// Standard interface — staggered (mixed) Poisson
static inline PetscInt idx_u(PetscInt ex,PetscInt ey,PetscInt Nx,PetscInt Ny){ return ey*Nx+ex; }
static inline PetscInt idx_gx(PetscInt fx,PetscInt fy,PetscInt Nx,PetscInt Ny){ return Nx*Ny+fy*(Nx+1)+fx; }
static inline PetscInt idx_gy(PetscInt fx,PetscInt fy,PetscInt Nx,PetscInt Ny){ return Nx*Ny+(Nx+1)*Ny+fx*(Ny+1)+fy; }
static inline PetscInt Nd(PetscInt Nx,PetscInt Ny){ return Nx*Ny+(Nx+1)*Ny+Nx*(Ny+1); }

PetscErrorCode SetInitialCondition(DM dm, Vec u, PetscReal t) {
  PetscFunctionBeginUser; (void)dm;(void)t; VecSet(u,0); PetscFunctionReturn(0);
}

PetscErrorCode AssembleSystem(DM dm, Mat A, PetscReal dt, const ManufacturedSolution &mms, const BoundaryCondition &bc) {
  PetscFunctionBeginUser; (void)dm;(void)dt;(void)mms;
  PetscInt Nx,Ny; DMStagGetGlobalSizes(dm,&Nx,&Ny,NULL);
  PetscReal h=1.0/Nx, ih=1.0/h, ih2=2.0/h;
  PetscInt rs,re; MatGetOwnershipRange(A,&rs,&re);
  PetscBool allN=(bc.left==BC_NEUMANN && bc.right==BC_NEUMANN && bc.bottom==BC_NEUMANN && bc.top==BC_NEUMANN);

  // ---- gx equations: gx - Gx u = boundary terms ----
  for(PetscInt fy=0;fy<Ny;fy++) for(PetscInt fx=0;fx<=Nx;fx++){
    PetscInt row=idx_gx(fx,fy,Nx,Ny); if(row<rs||row>=re) continue;
    if(fx==0){
      if(bc.left==BC_NEUMANN) {
        MatSetValue(A,row,row,1.0,INSERT_VALUES);  // gx = -g (identity, no u connection)
      } else {
        MatSetValue(A,row,row,1.0,INSERT_VALUES);
        MatSetValue(A,row,idx_u(0,fy,Nx,Ny),-ih2,INSERT_VALUES);  // gx - 2u/h = -2g/h
      }
    } else if(fx==Nx){
      if(bc.right==BC_NEUMANN) {
        MatSetValue(A,row,row,1.0,INSERT_VALUES);  // gx = +g
      } else {
        MatSetValue(A,row,row,1.0,INSERT_VALUES);
        MatSetValue(A,row,idx_u(Nx-1,fy,Nx,Ny),ih2,INSERT_VALUES); // gx + 2u/h = 2g/h
      }
    } else {
      MatSetValue(A,row,row,1.0,INSERT_VALUES);
      MatSetValue(A,row,idx_u(fx,fy,Nx,Ny),-ih,INSERT_VALUES);
      MatSetValue(A,row,idx_u(fx-1,fy,Nx,Ny),ih,INSERT_VALUES);
    }
  }

  // ---- gy equations: gy - Gy u = boundary terms ----
  for(PetscInt fx=0;fx<Nx;fx++) for(PetscInt fy=0;fy<=Ny;fy++){
    PetscInt row=idx_gy(fx,fy,Nx,Ny); if(row<rs||row>=re) continue;
    if(fy==0){
      if(bc.bottom==BC_NEUMANN) {
        MatSetValue(A,row,row,1.0,INSERT_VALUES);  // gy = -g
      } else {
        MatSetValue(A,row,row,1.0,INSERT_VALUES);
        MatSetValue(A,row,idx_u(fx,0,Nx,Ny),-ih2,INSERT_VALUES);  // gy - 2u/h = -2g/h
      }
    } else if(fy==Ny){
      if(bc.top==BC_NEUMANN) {
        MatSetValue(A,row,row,1.0,INSERT_VALUES);  // gy = +g
      } else {
        MatSetValue(A,row,row,1.0,INSERT_VALUES);
        MatSetValue(A,row,idx_u(fx,Ny-1,Nx,Ny),ih2,INSERT_VALUES); // gy + 2u/h = 2g/h
      }
    } else {
      MatSetValue(A,row,row,1.0,INSERT_VALUES);
      MatSetValue(A,row,idx_u(fx,fy,Nx,Ny),-ih,INSERT_VALUES);
      MatSetValue(A,row,idx_u(fx,fy-1,Nx,Ny),ih,INSERT_VALUES);
    }
  }

  // ---- u equations: -div(g) = f (with pin for all-Neumann) ----
  for(PetscInt ey=0;ey<Ny;ey++) for(PetscInt ex=0;ex<Nx;ex++){
    PetscInt row=idx_u(ex,ey,Nx,Ny); if(row<rs||row>=re) continue;
    if(allN && ex==0 && ey==0) {
      MatSetValue(A,row,row,1.0,INSERT_VALUES);  // pin u(0,0) = exact
      continue;
    }
    MatSetValue(A,row,idx_gx(ex+1,ey,Nx,Ny),-ih,INSERT_VALUES);
    MatSetValue(A,row,idx_gx(ex,ey,Nx,Ny),ih,INSERT_VALUES);
    MatSetValue(A,row,idx_gy(ex,ey+1,Nx,Ny),-ih,INSERT_VALUES);
    MatSetValue(A,row,idx_gy(ex,ey,Nx,Ny),ih,INSERT_VALUES);
    MatSetValue(A,row,row,1e-8,INSERT_VALUES);  // regularization for saddle-point
  }
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  PetscFunctionReturn(0);
}

PetscErrorCode BuildRHS(DM dm, Vec b, const Vec uOld, PetscReal t, PetscReal dt, const ManufacturedSolution &mms, const BoundaryCondition &bc) {
  PetscFunctionBeginUser; (void)uOld;(void)t;(void)dt;
  PetscInt Nx,Ny; DMStagGetGlobalSizes(dm,&Nx,&Ny,NULL);
  PetscReal h=1.0/Nx; VecSet(b,0);
  PetscInt rs,re; VecGetOwnershipRange(b,&rs,&re);
  PetscBool allN=(bc.left==BC_NEUMANN && bc.right==BC_NEUMANN && bc.bottom==BC_NEUMANN && bc.top==BC_NEUMANN);

  // gx RHS: boundary contributions
  for(PetscInt fy=0;fy<Ny;fy++){
    { // fx=0 (left)
      PetscInt row=idx_gx(0,fy,Nx,Ny); if(row>=rs&&row<re){
        PetscScalar v=0;
        if(bc.left==BC_NEUMANN) v=-EvalBC(bc.g_left,fy*h+0.5*h);  // gx = -g_left
        else                    v=-2.0*EvalBC(bc.g_left,fy*h+0.5*h)/h; // gx -2u/h = -2g/h
        VecSetValue(b,row,v,INSERT_VALUES);
      }
    }
    { // fx=Nx (right)
      PetscInt row=idx_gx(Nx,fy,Nx,Ny); if(row>=rs&&row<re){
        PetscScalar v=0;
        if(bc.right==BC_NEUMANN) v=EvalBC(bc.g_right,fy*h+0.5*h);  // gx = +g_right
        else                     v=2.0*EvalBC(bc.g_right,fy*h+0.5*h)/h; // gx+2u/h = 2g/h
        VecSetValue(b,row,v,INSERT_VALUES);
      }
    }
  }

  // gy RHS: boundary contributions
  for(PetscInt fx=0;fx<Nx;fx++){
    { // fy=0 (bottom)
      PetscInt row=idx_gy(fx,0,Nx,Ny); if(row>=rs&&row<re){
        PetscScalar v=0;
        if(bc.bottom==BC_NEUMANN) v=-EvalBC(bc.g_bottom,fx*h+0.5*h);  // gy = -g_bottom
        else                      v=-2.0*EvalBC(bc.g_bottom,fx*h+0.5*h)/h; // gy-2u/h = -2g/h
        VecSetValue(b,row,v,INSERT_VALUES);
      }
    }
    { // fy=Ny (top)
      PetscInt row=idx_gy(fx,Ny,Nx,Ny); if(row>=rs&&row<re){
        PetscScalar v=0;
        if(bc.top==BC_NEUMANN) v=EvalBC(bc.g_top,fx*h+0.5*h);  // gy = +g_top
        else                   v=2.0*EvalBC(bc.g_top,fx*h+0.5*h)/h; // gy+2u/h = 2g/h
        VecSetValue(b,row,v,INSERT_VALUES);
      }
    }
  }

  // u RHS: source f (and pin for all-Neumann)
  for(PetscInt ey=0;ey<Ny;ey++) for(PetscInt ex=0;ex<Nx;ex++){
    PetscInt row=idx_u(ex,ey,Nx,Ny); if(row<rs||row>=re) continue;
    if(allN && ex==0 && ey==0)
      VecSetValue(b,row,mms.u(0.5*h,0.5*h),INSERT_VALUES);  // pin u
    else
      VecSetValue(b,row,mms.f((ex+0.5)*h,(ey+0.5)*h),INSERT_VALUES);
  }
  VecAssemblyBegin(b); VecAssemblyEnd(b);
  PetscFunctionReturn(0);
}

PetscErrorCode ComputeError(DM dm, const Vec x, PetscReal t, const ManufacturedSolution &mms, PetscReal *err) {
  PetscFunctionBeginUser; (void)dm;(void)t;(void)mms;
  PetscInt Nx,Ny; DMStagGetGlobalSizes(dm,&Nx,&Ny,NULL);
  PetscReal h=1.0/Nx; double s=0;
  PetscInt rs,re; VecGetOwnershipRange(x,&rs,&re);
  for(PetscInt ey=0;ey<Ny;ey++) for(PetscInt ex=0;ex<Nx;ex++){
    PetscInt row=idx_u(ex,ey,Nx,Ny); if(row<rs||row>=re) continue;
    PetscScalar v; VecGetValues(x,1,&row,&v);
    double d=v-mms.u((ex+0.5)*h,(ey+0.5)*h); s+=d*d;
  }
  double g; MPI_Allreduce(&s,&g,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  *err=sqrt(g*h*h); PetscFunctionReturn(0);
}

int main(int argc,char**argv){
  PetscCall(PetscInitialize(&argc,&argv,NULL,"Staggered Poisson"));
  int rank; MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  PetscInt Nx=64,Ny=64; PetscBool chk=PETSC_FALSE,ct=PETSC_FALSE;
  char mms_name[32]="sinpi"; PetscBool flg; PetscCall(PetscOptionsGetString(NULL,NULL,"-mms",mms_name,sizeof(mms_name),&flg));
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
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-nx",&Nx,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-ny",&Ny,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-poisson_check_error",&chk,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-convergence_test",&ct,NULL));
  DM dm; DMStagCreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,Nx,Ny,PETSC_DECIDE,PETSC_DECIDE,0,0,0,DMSTAG_STENCIL_BOX,1,NULL,NULL,&dm);
  DMSetUp(dm);
  PetscInt Ntot=Nd(Nx,Ny);
  Mat A; MatCreate(PETSC_COMM_WORLD,&A); MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,Ntot,Ntot);
  MatSetUp(A); MatSeqAIJSetPreallocation(A,5,NULL); MatMPIAIJSetPreallocation(A,5,NULL,5,NULL);
  AssembleSystem(dm,A,0, *mms, bc);
  Vec x,b; MatCreateVecs(A,&x,&b); BuildRHS(dm,b,NULL,0,0, *mms, bc);
  KSP ksp; KSPCreate(PETSC_COMM_WORLD,&ksp); KSPSetOperators(ksp,A,A);
  KSPSetTolerances(ksp,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,5000);
  {PC pc;KSPGetPC(ksp,&pc);PCSetType(pc,PCLU);}  // saddle-point needs direct solve
  KSPSetFromOptions(ksp);
  if(rank==0) printf("Stagg Poisson: N=%dx%d h=%g DOFs=%d\n",Nx,Ny,1.0/Nx,Ntot);
  KSPSolve(ksp,b,x);
  PetscReal err=0; if(chk){ ComputeError(dm,x,0,*mms,&err); if(rank==0) printf("||u-u_ex||=%g\n",err); }
  if(ct&&rank==0) printf("CONVERGENCE: %d %g %g 0\n",Nx,1.0/Nx,err);
  KSPDestroy(&ksp); VecDestroy(&x); VecDestroy(&b); MatDestroy(&A); DMDestroy(&dm);
  PetscFinalize(); return 0;
}
