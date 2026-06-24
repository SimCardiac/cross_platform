#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include "analytical/unsteady.h"
#include "common/boundary.h"
#include "common/mms.h"

// ============================================================================
// Standardized solver interface (steady: t=0, dt=0)
// ============================================================================

// --- Set initial/exact condition at time t ---
PetscErrorCode SetInitialCondition(DM dm, Vec u, PetscReal t, const ManufacturedSolution &mms) {
  PetscFunctionBeginUser;
  (void)t;
  PetscInt Nx,Ny; DMDAGetInfo(dm,NULL,&Nx,&Ny,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
  PetscReal hx=1.0/(Nx-1), hy=1.0/(Ny-1);
  VecSet(u,0);
  PetscInt rs,re; VecGetOwnershipRange(u,&rs,&re);
  for(PetscInt j=0;j<Ny;j++) for(PetscInt i=0;i<Nx;i++){
    PetscInt row=j*Nx+i; if(row<rs||row>=re) continue;
    VecSetValue(u,row,mms.u(i*hx,j*hy),INSERT_VALUES);
  }
  VecAssemblyBegin(u); VecAssemblyEnd(u);
  PetscFunctionReturn(0);
}

// --- Assemble matrix (BC in boundary rows, 5-point Laplacian interior) ---
PetscErrorCode AssembleSystem(DM dm, Mat A, PetscReal dt, const ManufacturedSolution &mms,
                               const BoundaryCondition &bc) {
  PetscFunctionBeginUser;
  (void)dt;
  PetscInt Nx, Ny; DMDAGetInfo(dm, NULL, &Nx, &Ny, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  PetscReal hx=1.0/(Nx-1), hy=1.0/(Ny-1), ix2=1.0/(hx*hx), iy2=1.0/(hy*hy), diag=2.0*(ix2+iy2);
  for(PetscInt j=0;j<Ny;j++) for(PetscInt i=0;i<Nx;i++){
    PetscInt row=j*Nx+i;
    PetscBool isBoundary=(i==0||i==Nx-1||j==0||j==Ny-1);
    if(isBoundary){
      BCType t=BC_DIRICHLET;
      if(i==0)t=bc.left; else if(i==Nx-1)t=bc.right; else if(j==0)t=bc.bottom; else t=bc.top;
      if(t==BC_DIRICHLET){
        MatSetValue(A,row,row,1.0,INSERT_VALUES);  // identity row
      } else {
        // Neumann: use boundary-aware stencil (no double-counting)
        PetscReal d=0, left=0, right=0, down=0, up=0;
        if(i==0)        { d+=2*ix2; right=-2*ix2; }  // Neumann left
        else if(i==Nx-1){ d+=2*ix2; left=-2*ix2;  }  // Neumann right
        else            { d+=2*ix2; left=-ix2; right=-ix2; }
        if(j==0)        { d+=2*iy2; up=-2*iy2;    }  // Neumann bottom
        else if(j==Ny-1){ d+=2*iy2; down=-2*iy2;  }  // Neumann top
        else            { d+=2*iy2; down=-iy2; up=-iy2; }
        // Pin one corner for all-Neumann (removes constant nullspace)
        PetscBool allN = (bc.left==BC_NEUMANN && bc.right==BC_NEUMANN && bc.bottom==BC_NEUMANN && bc.top==BC_NEUMANN);
        if(allN && i==0 && j==0) { d=1.0; left=right=down=up=0; }
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

// --- Build RHS ---
PetscErrorCode BuildRHS(DM dm, Vec b, const Vec uOld, PetscReal t, PetscReal dt, const ManufacturedSolution &mms,
                         const BoundaryCondition &bc) {
  PetscFunctionBeginUser;
  (void)uOld;(void)t;(void)dt;
  PetscInt Nx,Ny; DMDAGetInfo(dm,NULL,&Nx,&Ny,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
  PetscReal hx=1.0/(Nx-1), hy=1.0/(Ny-1);
  VecSet(b,0);
  PetscInt rs,re; VecGetOwnershipRange(b,&rs,&re);
  for(PetscInt j=0;j<Ny;j++) for(PetscInt i=0;i<Nx;i++){
    PetscInt row=j*Nx+i; if(row<rs||row>=re) continue;
    PetscBool isB=(i==0||i==Nx-1||j==0||j==Ny-1);
    if(!isB){
      VecSetValue(b,row,mms.f(i*hx,j*hx),INSERT_VALUES);
    } else {
      // Boundary: Dirichlet → g(coord); Neumann → source + flux
      BCType t=BC_DIRICHLET; BCValueFunc g=NULL; PetscScalar coord=0;
      if(i==0)     { t=bc.left;   g=bc.g_left;   coord=j*hy; }
      else if(i==Nx-1){ t=bc.right; g=bc.g_right; coord=j*hy; }
      else if(j==0) { t=bc.bottom; g=bc.g_bottom; coord=i*hx; }
      else          { t=bc.top;    g=bc.g_top;    coord=i*hx; }
      PetscScalar val=0;
      if(t==BC_DIRICHLET){
        val = EvalBC(g,coord);                     // u = g on boundary
      } else {
        // Neumann boundary (pin (0,0) for all-Neumann nullspace)
        PetscBool allN = (bc.left==BC_NEUMANN && bc.right==BC_NEUMANN && bc.bottom==BC_NEUMANN && bc.top==BC_NEUMANN);
        if(allN && i==0 && j==0)
          val = mms.u(0,0);
        else {
          val = mms.f(i*hx,j*hy);
          // Add flux for each Neumann edge: Au = f + 2g/h (ghost-cell correction)
          if(i==0)      val += 2.0*EvalBC(bc.g_left, j*hy)/hx;
          if(i==Nx-1)   val += 2.0*EvalBC(bc.g_right,j*hy)/hx;
          if(j==0)      val += 2.0*EvalBC(bc.g_bottom,i*hx)/hy;
          if(j==Ny-1)   val += 2.0*EvalBC(bc.g_top,  i*hx)/hy;
        }
      }
      VecSetValue(b,row,val,INSERT_VALUES);
    }
  }
  VecAssemblyBegin(b); VecAssemblyEnd(b);
  PetscFunctionReturn(0);
}

// --- Compute L2 error at time t ---
PetscErrorCode ComputeError(DM dm, const Vec u, PetscReal t, const ManufacturedSolution &mms, PetscReal *error) {
  PetscFunctionBeginUser;
  (void)t;
  PetscInt Nx,Ny; DMDAGetInfo(dm,NULL,&Nx,&Ny,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
  PetscReal hx=1.0/(Nx-1), hy=1.0/(Ny-1);
  Vec ue; VecDuplicate(u,&ue);
  PetscInt rs,re; VecGetOwnershipRange(u,&rs,&re);
  for(PetscInt j=0;j<Ny;j++) for(PetscInt i=0;i<Nx;i++){
    PetscInt row=j*Nx+i; if(row<rs||row>=re) continue;
    VecSetValue(ue,row,mms.u(i*hx,j*hy),INSERT_VALUES);
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
  char mms_name[32]="sinpi"; PetscBool flg; PetscCall(PetscOptionsGetString(NULL,NULL,"-mms",mms_name,sizeof(mms_name),&flg));
  const ManufacturedSolution *mms=&MMS_SINPI;
  if(strcmp(mms_name,"poly2")==0) mms=&MMS_POLY2;
  else if(strcmp(mms_name,"cospi")==0) mms=&MMS_COSPI;

  // BC selection: auto-detect from MMS, or override with -bc_type LRTB
  char bc_str[8]=""; PetscCall(PetscOptionsGetString(NULL,NULL,"-bc_type",bc_str,sizeof(bc_str),&flg));
  BCType bl=BC_DIRICHLET, br=BC_DIRICHLET, bb=BC_DIRICHLET, bt=BC_DIRICHLET;
  if(strlen(bc_str)==4){
    bl=(bc_str[0]=='N')?BC_NEUMANN:BC_DIRICHLET; br=(bc_str[1]=='N')?BC_NEUMANN:BC_DIRICHLET;
    bb=(bc_str[2]=='N')?BC_NEUMANN:BC_DIRICHLET; bt=(bc_str[3]=='N')?BC_NEUMANN:BC_DIRICHLET;
  } else if(strcmp(mms_name,"cospi")==0){
    // cospi: du/dn=0 → homogeneous Neumann by default
    bl=br=bb=bt=BC_NEUMANN;
  }
  BoundaryCondition bc = MakeBCFromMMS(*mms, bl, br, bb, bt);
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-nx",&Nx,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-ny",&Ny,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-poisson_check_error",&chk,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-convergence_test",&ct,NULL));
  DM dm; DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,Nx,Ny,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,&dm);
  DMSetUp(dm); DMDASetUniformCoordinates(dm,0,1,0,1,0,0);
  Mat A; DMCreateMatrix(dm,&A);
  MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);
  Vec u,b; DMCreateGlobalVector(dm,&u); DMCreateGlobalVector(dm,&b);
  AssembleSystem(dm,A,0, *mms, bc);
  BuildRHS(dm,b,NULL,0,0, *mms, bc);
  KSP ksp; KSPCreate(PETSC_COMM_WORLD,&ksp); KSPSetOperators(ksp,A,A);
  KSPSetTolerances(ksp,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,10000);
  KSPSetFromOptions(ksp);
  if(rank==0) printf("Vertex Poisson: N=%dx%d h=%g\n",Nx,Ny,1.0/(Nx-1));
  KSPSolve(ksp,b,u);
  PetscReal err=0; if(chk){ ComputeError(dm,u,0,*mms,&err); if(rank==0) printf("||u-u_ex||=%g\n",err); }
  if(ct&&rank==0) printf("CONVERGENCE: %d %g %g 0\n",Nx,1.0/(Nx-1),err);
  KSPDestroy(&ksp); VecDestroy(&u); VecDestroy(&b); MatDestroy(&A); DMDestroy(&dm);
  PetscFinalize(); return 0;
}
