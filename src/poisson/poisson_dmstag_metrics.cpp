#include "poisson/poisson_dmstag_metrics.hpp"
#include <mpi.h>
#include <cmath>
#include <algorithm>

namespace PoissonDMStagMetrics {

PetscErrorCode compute_l2_error(
    DM dm,
    Vec u,
    PetscScalar (*u_exact_func)(PetscScalar, PetscScalar),
    PetscReal* l2_error)
{
    PetscFunctionBeginUser;
    
    Vec uLocal, uExact, diff;
    PetscCall(DMGetLocalVector(dm, &uLocal));
    PetscCall(DMCreateGlobalVector(dm, &uExact));
    
    // Fill exact solution
    PetscScalar ***aUe;
    PetscScalar **cX, **cY;
    PetscInt startx, starty, nx, ny, nEx[2];
    PetscInt icenter, ip;
    
    PetscCall(DMGlobalToLocalBegin(dm, uExact, INSERT_VALUES, uLocal));
    PetscCall(DMGlobalToLocalEnd(dm, uExact, INSERT_VALUES, uLocal));
    PetscCall(DMStagVecGetArray(dm, uLocal, &aUe));
    PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                               &nEx[0], &nEx[1], NULL));
    PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
    PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));
    
    for (PetscInt ey = starty; ey < starty + ny; ++ey) {
        for (PetscInt ex = startx; ex < startx + nx; ++ex) {
            const PetscScalar x = cX[ex][icenter];
            const PetscScalar y = cY[ey][icenter];
            aUe[ey][ex][ip] = u_exact_func(x, y);
        }
    }
    
    PetscCall(DMStagVecRestoreArray(dm, uLocal, &aUe));
    PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
    PetscCall(DMLocalToGlobal(dm, uLocal, INSERT_VALUES, uExact));
    
    // Compute difference
    PetscCall(VecDuplicate(u, &diff));
    PetscCall(VecCopy(u, diff));
    PetscCall(VecAXPY(diff, -1.0, uExact));
    
    // Get L2 norm
    PetscReal nrm2;
    PetscCall(VecNorm(diff, NORM_2, &nrm2));
    
    // Get domain info for scaling
    PetscInt Nglob[2];
    PetscCall(DMStagGetGlobalSizes(dm, &Nglob[0], &Nglob[1], NULL));
    
    PetscReal bmin[3], bmax[3];
    PetscCall(DMGetBoundingBox(dm, bmin, bmax));
    const PetscReal hx = (bmax[0] - bmin[0]) / Nglob[0];
    const PetscReal hy = (bmax[1] - bmin[1]) / Nglob[1];
    
    *l2_error = nrm2 * std::sqrt(hx * hy);
    
    PetscCall(VecDestroy(&diff));
    PetscCall(VecDestroy(&uExact));
    PetscCall(DMRestoreLocalVector(dm, &uLocal));
    
    PetscFunctionReturn(0);
}

PetscErrorCode compute_linf_error(
    DM dm,
    Vec u,
    PetscScalar (*u_exact_func)(PetscScalar, PetscScalar),
    PetscReal* linf_error)
{
    PetscFunctionBeginUser;
    
    Vec uLocal;
    PetscCall(DMGetLocalVector(dm, &uLocal));
    PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
    PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
    
    PetscScalar ***aU;
    PetscScalar **cX, **cY;
    PetscInt startx, starty, nx, ny, nEx[2];
    PetscInt icenter, ip;
    
    PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
    PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                               &nEx[0], &nEx[1], NULL));
    PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
    PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));
    
    PetscReal local_max = 0.0;
    for (PetscInt ey = starty; ey < starty + ny; ++ey) {
        for (PetscInt ex = startx; ex < startx + nx; ++ex) {
            const PetscScalar x = cX[ex][icenter];
            const PetscScalar y = cY[ey][icenter];
            const PetscScalar u_computed = aU[ey][ex][ip];
            const PetscScalar u_exact = u_exact_func(x, y);
            const PetscReal error = std::abs(u_computed - u_exact);
            local_max = std::max(local_max, error);
        }
    }
    
    PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
    PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
    PetscCall(DMRestoreLocalVector(dm, &uLocal));
    
    // Global reduction
    PetscReal global_max;
    MPI_Allreduce(&local_max, &global_max, 1, MPIU_REAL, MPIU_MAX, PETSC_COMM_WORLD);
    *linf_error = global_max;
    
    PetscFunctionReturn(0);
}

PetscErrorCode compute_residual_norm(
    DM dm,
    Vec u,
    Vec f,
    PetscInt Nx,
    PetscInt Ny,
    PetscReal* residual_norm)
{
    PetscFunctionBeginUser;
    
    Vec uLocal, fLocal;
    PetscCall(DMGetLocalVector(dm, &uLocal));
    PetscCall(DMGetLocalVector(dm, &fLocal));
    
    PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
    PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
    PetscCall(DMGlobalToLocalBegin(dm, f, INSERT_VALUES, fLocal));
    PetscCall(DMGlobalToLocalEnd(dm, f, INSERT_VALUES, fLocal));
    
    PetscScalar ***aU, ***aF;
    PetscInt startx, starty, nx, ny, nEx[2];
    PetscInt icenter;
    
    PetscReal bmin[3], bmax[3];
    PetscCall(DMGetBoundingBox(dm, bmin, bmax));
    const PetscReal hx = (bmax[0] - bmin[0]) / Nx;
    const PetscReal hy = (bmax[1] - bmin[1]) / Ny;
    const PetscReal ix2 = 1.0 / (hx * hx);
    const PetscReal iy2 = 1.0 / (hy * hy);
    const PetscReal diag = 2.0 * ix2 + 2.0 * iy2;
    
    PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
    PetscCall(DMStagVecGetArray(dm, fLocal, &aF));
    PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                               &nEx[0], &nEx[1], NULL));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &icenter));
    
    PetscReal local_sum = 0.0;
    for (PetscInt ey = starty; ey < starty + ny; ++ey) {
        for (PetscInt ex = startx; ex < startx + nx; ++ex) {
            const PetscScalar uc = aU[ey][ex][icenter];
            
            // Get neighbor values with ghost cell BC (Dirichlet u=0)
            PetscScalar ul, ur, ud, uu;
            if (ex == 0) ul = -uc;
            else ul = aU[ey][ex - 1][icenter];
            
            if (ex == Nx - 1) ur = -uc;
            else ur = aU[ey][ex + 1][icenter];
            
            if (ey == 0) ud = -uc;
            else ud = aU[ey - 1][ex][icenter];
            
            if (ey == Ny - 1) uu = -uc;
            else uu = aU[ey + 1][ex][icenter];
            
            const PetscScalar Au = diag * uc - ix2 * (ul + ur) - iy2 * (ud + uu);
            const PetscScalar ff = aF[ey][ex][icenter];
            const PetscScalar r = ff - Au;
            local_sum += PetscRealPart(r * r);
        }
    }
    
    PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
    PetscCall(DMStagVecRestoreArray(dm, fLocal, &aF));
    PetscCall(DMRestoreLocalVector(dm, &uLocal));
    PetscCall(DMRestoreLocalVector(dm, &fLocal));
    
    PetscReal global_sum;
    MPI_Allreduce(&local_sum, &global_sum, 1, MPIU_REAL, MPIU_SUM, PETSC_COMM_WORLD);
    *residual_norm = std::sqrt(global_sum);
    
    PetscFunctionReturn(0);
}

} // namespace PoissonDMStagMetrics
