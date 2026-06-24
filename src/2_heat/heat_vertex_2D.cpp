#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include "analytical/unsteady.h"
#include "common/boundary.h"
#include "common/mms.h"

static const PetscScalar alpha = 0.1;  // thermal diffusivity

// ============================================================================
// Standard interface: heat equation ∂u/∂t = αΔu + f
// Implicit Euler: (I - αΔt Δ) u^{n+1} = u^n + Δt f
// ============================================================================

PetscErrorCode SetInitialCondition(DM dm, Vec u, PetscReal t, const ManufacturedSolution &mms) {
  PetscFunctionBeginUser;
  PetscInt Nx,Ny; DMDAGetInfo(dm,NULL,&Nx,&Ny,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
  PetscReal hx=1.0/(Nx-1), hy=1.0/(Ny-1);
  PetscInt xs,ys,xm,ym; DMDAGetCorners(dm,&xs,&ys,NULL,&xm,&ym,NULL);
  PetscScalar **arr; DMDAVecGetArray(dm,u,&arr);
  for(PetscInt j=ys;j<ys+ym;j++) for(PetscInt i=xs;i<xs+xm;i++)
    arr[j][i]=mms.u_td(i*hx,j*hy,t,alpha);  // time-dependent exact solution
  DMDAVecRestoreArray(dm,u,&arr);
  PetscFunctionReturn(0);
}

PetscErrorCode AssembleSystem(DM dm, Mat A, PetscReal dt, const ManufacturedSolution &mms,
                               const BoundaryCondition &bc) {
  PetscFunctionBeginUser;
  (void)mms;
  PetscInt Nx,Ny; DMDAGetInfo(dm,NULL,&Nx,&Ny,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
  PetscReal hx=1.0/(Nx-1), hy=1.0/(Ny-1);
  PetscReal c=alpha*dt, ix2=c/(hx*hx), iy2=c/(hy*hy), diag=1.0+2.0*(ix2+iy2);
  PetscInt rstart,rend; MatGetOwnershipRange(A,&rstart,&rend);
  for(PetscInt j=0;j<Ny;j++) for(PetscInt i=0;i<Nx;i++){
    PetscInt row=j*Nx+i; if(row<rstart||row>=rend) continue;
    PetscBool isB=(i==0||i==Nx-1||j==0||j==Ny-1);
    if(isB){
      BCType t=BC_DIRICHLET;
      if(i==0)t=bc.left; else if(i==Nx-1)t=bc.right; else if(j==0)t=bc.bottom; else t=bc.top;
      if(t==BC_DIRICHLET){
        MatSetValue(A,row,row,1.0,INSERT_VALUES);
      } else {
        // Neumann Helmholtz: same Laplacian stencil as Poisson, plus identity from I
        PetscReal d=1.0, left=0, right=0, down=0, up=0; // identity
        if(i==0)        { d+=2*ix2; right=-2*ix2; }
        else if(i==Nx-1){ d+=2*ix2; left=-2*ix2;  }
        else            { d+=2*ix2; left=-ix2; right=-ix2; }
        if(j==0)        { d+=2*iy2; up=-2*iy2;    }
        else if(j==Ny-1){ d+=2*iy2; down=-2*iy2;  }
        else            { d+=2*iy2; down=-iy2; up=-iy2; }
        MatSetValue(A,row,row,d,INSERT_VALUES);
        if(left)  MatSetValue(A,row,j*Nx+(i-1),left,INSERT_VALUES);
        if(right) MatSetValue(A,row,j*Nx+(i+1),right,INSERT_VALUES);
        if(down)  MatSetValue(A,row,(j-1)*Nx+i,down,INSERT_VALUES);
        if(up)    MatSetValue(A,row,(j+1)*Nx+i,up,INSERT_VALUES);
      }
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

PetscErrorCode BuildRHS(DM dm, Vec b, const Vec uOld, PetscReal t, PetscReal dt,
                         const ManufacturedSolution &mms, const BoundaryCondition &bc) {
  PetscFunctionBeginUser;
  PetscInt Nx,Ny; DMDAGetInfo(dm,NULL,&Nx,&Ny,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
  PetscReal hx=1.0/(Nx-1), hy=1.0/(Ny-1);
  // Time factor for unsteady BC: u_td/u_steady at a safe interior point
  PetscReal t_eval=t+dt, tf=1.0;
  if(mms.u_td) { PetscScalar us=mms.u(0.25,0.25); if(us!=0) tf=mms.u_td(0.25,0.25,t_eval,alpha)/us; }
  PetscInt rstart,rend; VecGetOwnershipRange(b,&rstart,&rend);
  for(PetscInt j=0;j<Ny;j++) for(PetscInt i=0;i<Nx;i++){
    PetscInt row=j*Nx+i; if(row<rstart||row>=rend) continue;
    PetscBool isB=(i==0||i==Nx-1||j==0||j==Ny-1);
    PetscScalar val;
    if(!isB){
      // Interior: u^n + dt*f (time-dependent source)
      PetscScalar un; VecGetValues(uOld,1,&row,&un);
      val = un + dt * mms.f_td(i*hx,j*hy,t_eval,alpha);
    } else {
      BCType t=BC_DIRICHLET; BCValueFunc g=NULL; PetscScalar coord=0;
      if(i==0)     { t=bc.left;   g=bc.g_left;   coord=j*hy; }
      else if(i==Nx-1){ t=bc.right; g=bc.g_right; coord=j*hy; }
      else if(j==0) { t=bc.bottom; g=bc.g_bottom; coord=i*hx; }
      else          { t=bc.top;    g=bc.g_top;    coord=i*hx; }
      if(t==BC_DIRICHLET){
        val = mms.u_td(i*hx,j*hy,t_eval,alpha);  // Dirichlet = exact at t^{n+1}
      } else {
        PetscScalar un; VecGetValues(uOld,1,&row,&un);
        val = un + dt * mms.f_td(i*hx,j*hy,t_eval,alpha);
        PetscReal c_=alpha*dt;
        if(i==0)      val += 2.0*c_*tf*EvalBC(bc.g_left,  j*hy)/hx;
        if(i==Nx-1)   val += 2.0*c_*tf*EvalBC(bc.g_right, j*hy)/hx;
        if(j==0)      val += 2.0*c_*tf*EvalBC(bc.g_bottom,i*hx)/hy;
        if(j==Ny-1)   val += 2.0*c_*tf*EvalBC(bc.g_top,   i*hx)/hy;
      }
    }
    VecSetValue(b,row,val,INSERT_VALUES);
  }
  VecAssemblyBegin(b); VecAssemblyEnd(b);
  PetscFunctionReturn(0);
}

PetscErrorCode ComputeError(DM dm, const Vec u, PetscReal t, const ManufacturedSolution &mms, PetscReal *error) {
  PetscFunctionBeginUser;
  PetscInt Nx,Ny; DMDAGetInfo(dm,NULL,&Nx,&Ny,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
  PetscReal hx=1.0/(Nx-1), hy=1.0/(Ny-1);
  Vec ue; VecDuplicate(u,&ue);
  PetscInt xs,ys,xm,ym; DMDAGetCorners(dm,&xs,&ys,NULL,&xm,&ym,NULL);
  PetscScalar **arr; DMDAVecGetArray(dm,ue,&arr);
  for(PetscInt j=ys;j<ys+ym;j++) for(PetscInt i=xs;i<xs+xm;i++)
    arr[j][i]=mms.u_td(i*hx,j*hy,t,alpha);
  DMDAVecRestoreArray(dm,ue,&arr);
  Vec d; VecDuplicate(u,&d); VecCopy(u,d); VecAXPY(d,-1.0,ue);
  PetscReal n2; VecNorm(d,NORM_2,&n2);
  *error = n2*sqrt(hx*hy);
  VecDestroy(&d); VecDestroy(&ue);
  PetscFunctionReturn(0);
}

int main(int argc,char**argv){
  PetscCall(PetscInitialize(&argc,&argv,NULL,"Vertex Heat"));
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

  PetscReal h=1.0/(Nx-1), dt=h*h; // Δt ∝ h²
  PetscInt Nsteps=std::max(1,(int)std::ceil(Tfinal/dt));
  PetscReal dtActual=Tfinal/Nsteps;

  DM dm; DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,Nx,Ny,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,&dm);
  DMSetUp(dm); DMDASetUniformCoordinates(dm,0,1,0,1,0,0);
  Mat A; DMCreateMatrix(dm,&A); MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
  Vec u,uOld,b; DMCreateGlobalVector(dm,&u); DMCreateGlobalVector(dm,&uOld); DMCreateGlobalVector(dm,&b);

  AssembleSystem(dm,A,dtActual,*mms,bc);
  SetInitialCondition(dm,u,0,*mms);

  KSP ksp; KSPCreate(PETSC_COMM_WORLD,&ksp); KSPSetOperators(ksp,A,A);
  KSPSetTolerances(ksp,1e-12,PETSC_DEFAULT,PETSC_DEFAULT,2000); KSPSetFromOptions(ksp);

  if(rank==0) printf("Vertex Heat: N=%dx%d h=%g dt=%g steps=%d T=%g\n",Nx,Ny,h,dtActual,Nsteps,Tfinal);

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
