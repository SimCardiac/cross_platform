#include <petscdm.h>
#include <petscdmstag.h>
#include <petscksp.h>
#include <petscsys.h>
#include <cmath>
#include <iostream>

// Item 7 — Stokes Projection Method (staggered grid)
// MMS: u=sin(πx)cos(πy)exp(-2π²νt), v=-cos(πx)sin(πy)exp(-2π²νt), p=0
static const PetscScalar nu = 0.1;

static inline PetscScalar u_ex(PetscScalar x, PetscScalar y, PetscScalar t) {
  return  std::sin(M_PI*x)*std::cos(M_PI*y)*std::exp(-2*M_PI*M_PI*nu*t);
}
static inline PetscScalar v_ex(PetscScalar x, PetscScalar y, PetscScalar t) {
  return -std::cos(M_PI*x)*std::sin(M_PI*y)*std::exp(-2*M_PI*M_PI*nu*t);
}

// ============================================================================
int main(int argc, char **argv) {
  PetscCall(PetscInitialize(&argc, &argv, NULL, "Stokes Projection"));
  int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscInt Nx=32, Ny=32; PetscReal Tf=0.01;
  PetscBool chk=PETSC_FALSE, ct=PETSC_FALSE;
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-nx",&Nx,NULL));
  PetscCall(PetscOptionsGetInt(NULL,NULL,"-ny",&Ny,NULL));
  PetscCall(PetscOptionsGetReal(NULL,NULL,"-T",&Tf,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-stokes_check_error",&chk,NULL));
  PetscCall(PetscOptionsGetBool(NULL,NULL,"-convergence_test",&ct,NULL));

  const PetscReal h=1.0/Nx, dt=h*h;
  PetscInt Ns=(PetscInt)ceil(Tf/dt); PetscReal dta=Tf/Ns;

  // DMStag: p=element, u=left-face, v=down-face
  DM dm; PetscInt su,sv,sp;
  PetscCall(DMStagCreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,
    Nx,Ny,PETSC_DECIDE,PETSC_DECIDE,1,1,1,DMSTAG_STENCIL_BOX,1,NULL,NULL,&dm));
  PetscCall(DMSetUp(dm));
  PetscCall(DMStagSetUniformCoordinatesProduct(dm,0,1,0,1,0,0));
  PetscCall(DMStagGetLocationSlot(dm,DMSTAG_LEFT,0,&su));
  PetscCall(DMStagGetLocationSlot(dm,DMSTAG_DOWN,0,&sv));
  PetscCall(DMStagGetLocationSlot(dm,DMSTAG_ELEMENT,0,&sp));

  Vec X,Xo; PetscCall(DMCreateGlobalVector(dm,&X)); PetscCall(DMCreateGlobalVector(dm,&Xo));
  Vec xl,xol; PetscCall(DMGetLocalVector(dm,&xl)); PetscCall(DMGetLocalVector(dm,&xol));
  Vec pl;     PetscCall(DMGetLocalVector(dm,&pl));

  PetscScalar ***a, ***ao, ***ap, **cX, **cY;
  PetscInt sx,sy,nx,ny,ne[2],ic,ip,id;
  PetscCall(DMStagGetCorners(dm,&sx,&sy,NULL,&nx,&ny,NULL,ne,ne+1,NULL));
  PetscCall(DMStagGetProductCoordinateArraysRead(dm,&cX,&cY,NULL));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm,DMSTAG_ELEMENT,&ic));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm,DMSTAG_LEFT,&ip));
  PetscCall(DMStagGetProductCoordinateLocationSlot(dm,DMSTAG_DOWN,&id));
  PetscInt xMax=sx+nx+ne[0]-1, yMax=sy+ny+ne[1]-1;

  // Init
  PetscCall(VecSet(X,0));
  PetscCall(DMGlobalToLocalBegin(dm,X,INSERT_VALUES,xl));
  PetscCall(DMGlobalToLocalEnd(dm,X,INSERT_VALUES,xl));
  PetscCall(DMStagVecGetArray(dm,xl,&a));
  for(PetscInt j=sy;j<=yMax;j++) for(PetscInt i=sx;i<=xMax;i++){
    a[j][i][su]=0; a[j][i][sv]=0; a[j][i][sp]=0;
  }
  // u at left faces: j 0..Ny-1, i 0..Nx
  for(PetscInt j=sy; j<sy+ny; j++) for(PetscInt i=sx; i<=sx+nx; i++)
    a[j][i][su]=u_ex(cX[i][ip],cY[j][ic],0);
  // v at down faces: j 0..Ny, i 0..Nx-1
  for(PetscInt j=sy; j<=sy+ny; j++) for(PetscInt i=sx; i<sx+nx; i++)
    a[j][i][sv]=v_ex(cX[i][ic],cY[j][id],0);
  PetscCall(DMStagVecRestoreArray(dm,xl,&a));
  PetscCall(DMLocalToGlobal(dm,xl,INSERT_VALUES,X));

  if(rank==0) printf("Stokes: N=%dx%d h=%g dt=%g steps=%d\n",Nx,Ny,h,dta,Ns);

  PetscReal t=0;
  const PetscReal c=nu*dta/(h*h), iha=1/h, diag=1+4*c;
  PetscInt GSu=50, GSv=50, GSp=5000;  // GS iterations per step

  for(PetscInt step=1; step<=Ns; step++){ t+=dta;
    PetscCall(VecCopy(X,Xo));

    // --- Helmholtz u (left faces, interior i=1..Nx-1, j=0..Ny-1) ---
    for(PetscInt gs=0; gs<GSu; gs++){
      for(int col=0;col<2;col++){
        PetscCall(DMGlobalToLocalBegin(dm,X,INSERT_VALUES,xl));
        PetscCall(DMGlobalToLocalEnd(dm,X,INSERT_VALUES,xl));
        PetscCall(DMGlobalToLocalBegin(dm,Xo,INSERT_VALUES,xol));
        PetscCall(DMGlobalToLocalEnd(dm,Xo,INSERT_VALUES,xol));
        PetscCall(DMStagVecGetArray(dm,xl,&a));
        PetscCall(DMStagVecGetArray(dm,xol,&ao));
        for(PetscInt j=sy; j<sy+ny; j++) for(PetscInt i=sx+1; i<sx+nx; i++){
          if(((i+j)&1)!=col) continue;
          PetscScalar ul=(i==sx+1)?0:a[j][i-1][su];
          PetscScalar ur=(i==sx+nx-1)?0:a[j][i+1][su];
          PetscScalar ud;
          if(j==sy)     ud= a[j][i][su];
          else          ud= a[j-1][i][su];
          PetscScalar uu;
          if(j==sy+ny-1) uu= a[j][i][su];
          else           uu= a[j+1][i][su];
          a[j][i][su]=(ao[j][i][su]+c*(ul+ur+ud+uu))/diag;
        }
        PetscCall(DMStagVecRestoreArray(dm,xol,&ao));
        PetscCall(DMStagVecRestoreArray(dm,xl,&a));
        PetscCall(DMLocalToGlobal(dm,xl,INSERT_VALUES,X));
      }
    }

    // --- Helmholtz v (down faces, interior j=1..Ny-1, i=0..Nx-1) ---
    for(PetscInt gs=0; gs<GSv; gs++){
      for(int col=0;col<2;col++){
        PetscCall(DMGlobalToLocalBegin(dm,X,INSERT_VALUES,xl));
        PetscCall(DMGlobalToLocalEnd(dm,X,INSERT_VALUES,xl));
        PetscCall(DMGlobalToLocalBegin(dm,Xo,INSERT_VALUES,xol));
        PetscCall(DMGlobalToLocalEnd(dm,Xo,INSERT_VALUES,xol));
        PetscCall(DMStagVecGetArray(dm,xl,&a));
        PetscCall(DMStagVecGetArray(dm,xol,&ao));
        for(PetscInt j=sy+1; j<sy+ny; j++) for(PetscInt i=sx; i<sx+nx; i++){
          if(((i+j)&1)!=col) continue;
          PetscScalar ul,ur;
          if(i==sx)     ul= a[j][i][sv];
          else          ul= a[j][i-1][sv];
          if(i==sx+nx-1) ur= a[j][i][sv];
          else           ur= a[j][i+1][sv];
          PetscScalar ud=(j==sy+1)?0:a[j-1][i][sv];
          PetscScalar uu=(j==sy+ny-1)?0:a[j+1][i][sv];
          a[j][i][sv]=(ao[j][i][sv]+c*(ul+ur+ud+uu))/diag;
        }
        PetscCall(DMStagVecRestoreArray(dm,xol,&ao));
        PetscCall(DMStagVecRestoreArray(dm,xl,&a));
        PetscCall(DMLocalToGlobal(dm,xl,INSERT_VALUES,X));
      }
    }

    // --- Pressure Poisson: Δp = (1/dt)∇·u* (GS, Neumann BC) ---
    // Zero p only (not u and v!)
    PetscCall(DMGlobalToLocalBegin(dm,X,INSERT_VALUES,xl));
    PetscCall(DMGlobalToLocalEnd(dm,X,INSERT_VALUES,xl));
    PetscCall(DMStagVecGetArray(dm,xl,&a));
    for(PetscInt j=sy; j<sy+ny; j++) for(PetscInt i=sx; i<sx+nx; i++) a[j][i][sp]=0;
    PetscCall(DMStagVecRestoreArray(dm,xl,&a));
    PetscCall(DMLocalToGlobal(dm,xl,INSERT_VALUES,X));

    PetscReal ih2=1.0/(h*h), pdiag=4*ih2;
    for(PetscInt gs=0; gs<GSp; gs++){
      for(int col=0;col<2;col++){
        PetscCall(DMGlobalToLocalBegin(dm,X,INSERT_VALUES,xl));
        PetscCall(DMGlobalToLocalEnd(dm,X,INSERT_VALUES,xl));
        PetscCall(DMStagVecGetArray(dm,xl,&a));
        for(PetscInt j=sy; j<sy+ny; j++) for(PetscInt i=sx; i<sx+nx; i++){
          if(((i+j)&1)!=col) continue;
          PetscScalar div=(a[j][i+1][su]-a[j][i][su])*iha
                         +(a[j+1][i][sv]-a[j][i][sv])*iha;
          PetscScalar pl=(i==sx)?a[j][i][sp]:a[j][i-1][sp];
          PetscScalar pr=(i==sx+nx-1)?a[j][i][sp]:a[j][i+1][sp];
          PetscScalar pd=(j==sy)?a[j][i][sp]:a[j-1][i][sp];
          PetscScalar pu=(j==sy+ny-1)?a[j][i][sp]:a[j+1][i][sp];
          a[j][i][sp]=(ih2*(pl+pr+pd+pu)-div/dta)/pdiag;
        }
        PetscCall(DMStagVecRestoreArray(dm,xl,&a));
        PetscCall(DMLocalToGlobal(dm,xl,INSERT_VALUES,X));
      }
    }

    // --- Velocity correction: u^{n+1}=u*-dt*∇p, v^{n+1}=v*-dt*∇p ---
    PetscCall(DMGlobalToLocalBegin(dm,X,INSERT_VALUES,xl));
    PetscCall(DMGlobalToLocalEnd(dm,X,INSERT_VALUES,xl));
    PetscCall(DMStagVecGetArray(dm,xl,&a));
    // u correction (interior left faces)
    for(PetscInt j=sy; j<sy+ny; j++) for(PetscInt i=sx+1; i<sx+nx; i++)
      a[j][i][su] -= dta*(a[j][i][sp]-a[j][i-1][sp])*iha;
    // v correction (interior down faces)
    for(PetscInt j=sy+1; j<sy+ny; j++) for(PetscInt i=sx; i<sx+nx; i++)
      a[j][i][sv] -= dta*(a[j][i][sp]-a[j-1][i][sp])*iha;
    PetscCall(DMStagVecRestoreArray(dm,xl,&a));
    PetscCall(DMLocalToGlobal(dm,xl,INSERT_VALUES,X));

    if(rank==0 && (step%100==0||step==Ns))
      printf("  step %d/%d\n",step,Ns);
  }

  // --- Errors ---
  PetscReal eu=0,ev=0;
  if(chk){
    PetscCall(DMGlobalToLocalBegin(dm,X,INSERT_VALUES,xl));
    PetscCall(DMGlobalToLocalEnd(dm,X,INSERT_VALUES,xl));
    PetscCall(DMStagVecGetArray(dm,xl,&a));
    double su2=0,sv2=0;
    for(PetscInt j=sy; j<sy+ny; j++) for(PetscInt i=sx+1; i<sx+nx; i++){
      double d=a[j][i][su]-u_ex(cX[i][ip],cY[j][ic],t); su2+=d*d;
    }
    for(PetscInt j=sy+1; j<sy+ny; j++) for(PetscInt i=sx; i<sx+nx; i++){
      double d=a[j][i][sv]-v_ex(cX[i][ic],cY[j][id],t); sv2+=d*d;
    }
    PetscCall(DMStagVecRestoreArray(dm,xl,&a));
    double gu,gv; MPI_Allreduce(&su2,&gu,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&sv2,&gv,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    eu=sqrt(gu*h*h); ev=sqrt(gv*h*h);
    if(rank==0){ printf("||u-u_ex||=%g\n||v-v_ex||=%g\n",eu,ev); }
  }
  if(ct && rank==0) printf("CONVERGENCE: %d %g %g %g %d\n",Nx,h,eu,ev,Ns);

  PetscCall(DMStagRestoreProductCoordinateArraysRead(dm,&cX,&cY,NULL));
  PetscCall(DMRestoreLocalVector(dm,&xl)); PetscCall(DMRestoreLocalVector(dm,&xol));
  PetscCall(DMRestoreLocalVector(dm,&pl));
  PetscCall(VecDestroy(&X)); PetscCall(VecDestroy(&Xo));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFinalize());
  return 0;
}
