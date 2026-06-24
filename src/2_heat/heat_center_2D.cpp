#include <petscdm.h>
#include <petscdmstag.h>
#include <petscksp.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include "common/boundary.h"
#include "common/mms.h"

static const PetscScalar alpha = 0.1;
static inline PetscInt idx(PetscInt ex,PetscInt ey,PetscInt Nx,PetscInt Ny){ (void)Ny; return ey*Nx+ex; }

// ============================================================================
// Standard interface: heat ∂u/∂t = αΔu + f, implicit Euler
// ============================================================================

PetscErrorCode SetInitialCondition(DM dm, Vec u, PetscReal t, const ManufacturedSolution &mms) {
  PetscFunctionBeginUser;
  PetscInt Nx,Ny; DMStagGetGlobalSizes(dm,&Nx,&Ny,NULL);
  PetscReal hx=1.0/Nx, hy=1.0/Ny;
  for(PetscInt ey=0;ey<Ny;ey++) for(PetscInt ex=0;ex<Nx;ex++)
    VecSetValue(u,idx(ex,ey,Nx,Ny),mms.u_td((ex+0.5)*hx,(ey+0.5)*hy,t,alpha),INSERT_VALUES);
  VecAssemblyBegin(u); VecAssemblyEnd(u);
  PetscFunctionReturn(0);
}

PetscErrorCode AssembleSystem(DM dm, Mat A, PetscReal dt, const ManufacturedSolution &mms,
                               const BoundaryCondition &bc) {
  PetscFunctionBeginUser; (void)dm;(void)mms;
  PetscInt Nx,Ny; DMStagGetGlobalSizes(dm,&Nx,&Ny,NULL);
  PetscReal hx=1.0/Nx, hy=1.0/Ny, c=alpha*dt, ix2=c/(hx*hx), iy2=c/(hy*hy), diag=1.0+2.0*(ix2+iy2);
  for(PetscInt ey=0;ey<Ny;ey++) for(PetscInt ex=0;ex<Nx;ex++){
    PetscInt row=idx(ex,ey,Nx,Ny);
    PetscReal d=diag;
    if(ex==0) { if(bc.left==BC_DIRICHLET) d+=ix2; else d-=ix2; }
    if(ex==Nx-1) { if(bc.right==BC_DIRICHLET) d+=ix2; else d-=ix2; }
    if(ey==0) { if(bc.bottom==BC_DIRICHLET) d+=iy2; else d-=iy2; }
    if(ey==Ny-1) { if(bc.top==BC_DIRICHLET) d+=iy2; else d-=iy2; }
    MatSetValue(A,row,row,d,INSERT_VALUES);
    if(ex>0)    MatSetValue(A,row,idx(ex-1,ey,Nx,Ny),-ix2,INSERT_VALUES);
    if(ex<Nx-1) MatSetValue(A,row,idx(ex+1,ey,Nx,Ny),-ix2,INSERT_VALUES);
    if(ey>0)    MatSetValue(A,row,idx(ex,ey-1,Nx,Ny),-iy2,INSERT_VALUES);
    if(ey<Ny-1) MatSetValue(A,row,idx(ex,ey+1,Nx,Ny),-iy2,INSERT_VALUES);
  }
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  PetscFunctionReturn(0);
}

PetscErrorCode BuildRHS(DM dm, Vec b, const Vec uOld, PetscReal t, PetscReal dt,
                         const ManufacturedSolution &mms, const BoundaryCondition &bc) {
  PetscFunctionBeginUser;
  PetscInt Nx,Ny; DMStagGetGlobalSizes(dm,&Nx,&Ny,NULL);
  PetscReal hx=1.0/Nx, hy=1.0/Ny, ix2=1.0/(hx*hx), iy2=1.0/(hy*hy), c=alpha*dt;
  PetscReal t_eval=t+dt, tf=1.0;
  if(mms.u_td){PetscScalar us=mms.u(0.25,0.25); if(us!=0) tf=mms.u_td(0.25,0.25,t_eval,alpha)/us;}
  VecSet(b,0);
  for(PetscInt ey=0;ey<Ny;ey++) for(PetscInt ex=0;ex<Nx;ex++){
    PetscScalar x=(ex+0.5)*hx, y=(ey+0.5)*hy;
    PetscInt row=idx(ex,ey,Nx,Ny);
    PetscScalar un; VecGetValues(uOld,1,&row,&un);
    PetscScalar val=un+dt*mms.f_td(x,y,t_eval,alpha);
    if(ex==0)   { if(bc.left==BC_DIRICHLET)   val+=2*c*tf*EvalBC(bc.g_left,y)*ix2;   else val+=c*tf*EvalBC(bc.g_left,y)/hx; }
    if(ex==Nx-1){ if(bc.right==BC_DIRICHLET)  val+=2*c*tf*EvalBC(bc.g_right,y)*ix2;  else val+=c*tf*EvalBC(bc.g_right,y)/hx; }
    if(ey==0)   { if(bc.bottom==BC_DIRICHLET) val+=2*c*tf*EvalBC(bc.g_bottom,x)*iy2; else val+=c*tf*EvalBC(bc.g_bottom,x)/hy; }
    if(ey==Ny-1){ if(bc.top==BC_DIRICHLET)    val+=2*c*tf*EvalBC(bc.g_top,x)*iy2;    else val+=c*tf*EvalBC(bc.g_top,x)/hy; }
    VecSetValue(b,row,val,INSERT_VALUES);
  }
  VecAssemblyBegin(b); VecAssemblyEnd(b);
  PetscFunctionReturn(0);
}

PetscErrorCode ComputeError(DM dm, const Vec u, PetscReal t, const ManufacturedSolution &mms, PetscReal *error) {
  PetscFunctionBeginUser;
  PetscInt Nx,Ny; DMStagGetGlobalSizes(dm,&Nx,&Ny,NULL);
  PetscReal hx=1.0/Nx, hy=1.0/Ny; double s=0;
  for(PetscInt ey=0;ey<Ny;ey++) for(PetscInt ex=0;ex<Nx;ex++){
    PetscScalar v; PetscInt row=idx(ex,ey,Nx,Ny); VecGetValues(u,1,&row,&v);
    double d=v-mms.u_td((ex+0.5)*hx,(ey+0.5)*hy,t,alpha); s+=d*d;
  }
  double g; MPI_Allreduce(&s,&g,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  *error=sqrt(g*hx*hy); PetscFunctionReturn(0);
}

int main(int argc,char**argv){
  PetscCall(PetscInitialize(&argc,&argv,NULL,"Cell Heat"));
  int rank; MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  PetscInt Nx=64,Ny=64; PetscReal Tfinal=0.05;
  PetscBool chk=PETSC_FALSE,ct=PETSC_FALSE;
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
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-nx",&Nx,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-ny",&Ny,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-T",&Tfinal,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-heat_check_error",&chk,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-convergence_test",&ct,NULL));

  PetscReal h=1.0/Nx, dt=h*h;
  PetscInt Nsteps=std::max(1,(int)std::ceil(Tfinal/dt));
  PetscReal dtActual=Tfinal/Nsteps;

  PetscInt Ntot=Nx*Ny;
  Mat A; MatCreate(PETSC_COMM_WORLD,&A); MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,Ntot,Ntot);
  MatSetUp(A); MatSeqAIJSetPreallocation(A,5,NULL); MatMPIAIJSetPreallocation(A,5,NULL,5,NULL);
  DM dm; DMStagCreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,Nx,Ny,PETSC_DECIDE,PETSC_DECIDE,1,0,0,DMSTAG_STENCIL_BOX,1,NULL,NULL,&dm);
  DMSetUp(dm); DMStagSetUniformCoordinatesProduct(dm,0,1,0,1,0,0);
  AssembleSystem(dm,A,dtActual,*mms,bc);
  Vec u,uOld,b; MatCreateVecs(A,&u,&b); VecDuplicate(u,&uOld);
  SetInitialCondition(dm,u,0,*mms);

  KSP ksp; KSPCreate(PETSC_COMM_WORLD,&ksp); KSPSetOperators(ksp,A,A);
  KSPSetTolerances(ksp,1e-12,PETSC_DEFAULT,PETSC_DEFAULT,2000); KSPSetFromOptions(ksp);

  if(rank==0) printf("Cell Heat: N=%dx%d h=%g dt=%g steps=%d T=%g\n",Nx,Ny,h,dtActual,Nsteps,Tfinal);

  for(PetscInt step=1;step<=Nsteps;step++){
    PetscReal tn=dtActual*(step-1);
    VecCopy(u,uOld);
    BuildRHS(dm,b,uOld,tn,dtActual,*mms,bc);
    KSPSolve(ksp,b,u);
  }

  PetscReal err=0; if(chk){ ComputeError(dm,u,Tfinal,*mms,&err); if(rank==0) printf("||u-u_ex||=%g\n",err); }
  if(ct&&rank==0) printf("CONVERGENCE: %d %g %g %d\n",Nx,h,err,Nsteps);
  KSPDestroy(&ksp); VecDestroy(&u); VecDestroy(&uOld); VecDestroy(&b); MatDestroy(&A); DMDestroy(&dm);
  PetscFinalize(); return 0;
}
