#include "poisson/poisson_dmstag_solver.hpp"
#include <iostream>
#include <cmath>
#include <string>

namespace PoissonDMStagSolver {

// Setup RHS vector
static PetscErrorCode SetupRHS(
    DM dm, Vec f, Vec fLocal,
    PetscScalar (*rhs_func)(PetscScalar, PetscScalar))
{
    PetscFunctionBeginUser;
    
    PetscScalar ***aF;
    PetscScalar **cX, **cY;
    PetscInt startx, starty, nx, ny, nEx[2];
    PetscInt icenter, ip;
    
    PetscCall(DMGlobalToLocalBegin(dm, f, INSERT_VALUES, fLocal));
    PetscCall(DMGlobalToLocalEnd(dm, f, INSERT_VALUES, fLocal));
    PetscCall(DMStagVecGetArray(dm, fLocal, &aF));
    PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                               &nEx[0], &nEx[1], NULL));
    PetscCall(DMStagGetProductCoordinateArraysRead(dm, &cX, &cY, NULL));
    PetscCall(DMStagGetProductCoordinateLocationSlot(dm, DMSTAG_ELEMENT, &icenter));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &ip));
    
    for (PetscInt ey = starty; ey < starty + ny; ++ey) {
        for (PetscInt ex = startx; ex < startx + nx; ++ex) {
            const PetscScalar x = cX[ex][icenter];
            const PetscScalar y = cY[ey][icenter];
            aF[ey][ex][ip] = rhs_func(x, y);
        }
    }
    
    PetscCall(DMStagVecRestoreArray(dm, fLocal, &aF));
    PetscCall(DMStagRestoreProductCoordinateArraysRead(dm, &cX, &cY, NULL));
    PetscCall(DMLocalToGlobal(dm, fLocal, INSERT_VALUES, f));
    
    PetscFunctionReturn(0);
}

// Gauss-Seidel sweep (red-black)
static PetscErrorCode GaussSeidelSweep(
    DM dm, Vec u, Vec uLocal, Vec f, Vec fLocal,
    PetscInt Nx, PetscInt Ny,
    PetscReal hx, PetscReal hy)
{
    PetscFunctionBeginUser;
    
    const PetscReal ix2 = 1.0 / (hx * hx);
    const PetscReal iy2 = 1.0 / (hy * hy);
    const PetscReal diag = 2.0 * ix2 + 2.0 * iy2;
    
    for (int color = 0; color < 2; ++color) {
        PetscScalar ***aU, ***aF;
        PetscInt startx, starty, nx, ny, nEx[2];
        PetscInt icenter;
        
        PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
        PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
        PetscCall(DMGlobalToLocalBegin(dm, f, INSERT_VALUES, fLocal));
        PetscCall(DMGlobalToLocalEnd(dm, f, INSERT_VALUES, fLocal));
        PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
        PetscCall(DMStagVecGetArray(dm, fLocal, &aF));
        PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                                   &nEx[0], &nEx[1], NULL));
        PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &icenter));
        
        for (PetscInt ey = starty; ey < starty + ny; ++ey) {
            for (PetscInt ex = startx; ex < startx + nx; ++ex) {
                if (((ex + ey) & 1) != color) continue;
                
                PetscScalar ul, ur, ud, uu;
                
                // Ghost cell BC (Dirichlet u=0)
                if (ex == 0) ul = -aU[ey][ex][icenter];
                else ul = aU[ey][ex - 1][icenter];
                
                if (ex == Nx - 1) ur = -aU[ey][ex][icenter];
                else ur = aU[ey][ex + 1][icenter];
                
                if (ey == 0) ud = -aU[ey][ex][icenter];
                else ud = aU[ey - 1][ex][icenter];
                
                if (ey == Ny - 1) uu = -aU[ey][ex][icenter];
                else uu = aU[ey + 1][ex][icenter];
                
                const PetscScalar ff = aF[ey][ex][icenter];
                const PetscScalar unew = (ix2 * (ul + ur) + iy2 * (ud + uu) + ff) / diag;
                aU[ey][ex][icenter] = unew;
            }
        }
        
        PetscCall(DMStagVecRestoreArray(dm, uLocal, &aU));
        PetscCall(DMStagVecRestoreArray(dm, fLocal, &aF));
        PetscCall(DMLocalToGlobal(dm, uLocal, INSERT_VALUES, u));
    }
    
    PetscFunctionReturn(0);
}

// Compute residual norm
static PetscErrorCode ComputeResidualNorm(
    DM dm, Vec u, Vec uLocal, Vec f, Vec fLocal,
    PetscInt Nx, PetscInt Ny,
    PetscReal hx, PetscReal hy,
    PetscReal* residual_norm)
{
    PetscFunctionBeginUser;
    
    const PetscReal ix2 = 1.0 / (hx * hx);
    const PetscReal iy2 = 1.0 / (hy * hy);
    const PetscReal diag = 2.0 * ix2 + 2.0 * iy2;
    
    PetscScalar ***aU, ***aF;
    PetscInt startx, starty, nx, ny, nEx[2];
    PetscInt icenter;
    
    PetscCall(DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uLocal));
    PetscCall(DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uLocal));
    PetscCall(DMGlobalToLocalBegin(dm, f, INSERT_VALUES, fLocal));
    PetscCall(DMGlobalToLocalEnd(dm, f, INSERT_VALUES, fLocal));
    PetscCall(DMStagVecGetArray(dm, uLocal, &aU));
    PetscCall(DMStagVecGetArray(dm, fLocal, &aF));
    PetscCall(DMStagGetCorners(dm, &startx, &starty, NULL, &nx, &ny, NULL,
                               &nEx[0], &nEx[1], NULL));
    PetscCall(DMStagGetLocationSlot(dm, DMSTAG_ELEMENT, 0, &icenter));
    
    PetscReal local_sum = 0.0;
    for (PetscInt ey = starty; ey < starty + ny; ++ey) {
        for (PetscInt ex = startx; ex < startx + nx; ++ex) {
            const PetscScalar uc = aU[ey][ex][icenter];
            
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
    
    PetscReal global_sum;
    MPI_Allreduce(&local_sum, &global_sum, 1, MPIU_REAL, MPIU_SUM, PETSC_COMM_WORLD);
    *residual_norm = std::sqrt(global_sum);
    
    PetscFunctionReturn(0);
}

PetscErrorCode run_poisson_case(
    DM* dm_out,
    Vec* u_out,
    Vec* f_out,
    PetscInt nx, PetscInt ny,
    PetscReal x0, PetscReal x1,
    PetscReal y0, PetscReal y1,
    PetscScalar (*f_rhs)(PetscScalar, PetscScalar),
    PetscScalar (*u_exact)(PetscScalar, PetscScalar),
    const char* bc_west,
    const char* bc_east,
    const char* bc_south,
    const char* bc_north,
    PetscReal tol,
    PetscInt max_its)
{
    PetscFunctionBeginUser;
    
    // Create DMStag
    DM dm;
    const PetscInt dof0 = 0, dof1 = 0, dof2 = 1;
    const PetscInt stencilWidth = 1;
    
    PetscCall(DMStagCreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                             nx, ny, PETSC_DECIDE, PETSC_DECIDE,
                             dof0, dof1, dof2, DMSTAG_STENCIL_BOX, stencilWidth,
                             NULL, NULL, &dm));
    PetscCall(DMSetUp(dm));
    PetscCall(DMStagSetUniformCoordinatesProduct(dm, x0, x1, y0, y1, 0.0, 0.0));
    
    // Create vectors
    Vec u, f, uLocal, fLocal;
    PetscCall(DMCreateGlobalVector(dm, &u));
    PetscCall(DMCreateGlobalVector(dm, &f));
    PetscCall(DMGetLocalVector(dm, &uLocal));
    PetscCall(DMGetLocalVector(dm, &fLocal));
    PetscCall(VecSet(u, 0.0));
    
    // Setup RHS
    PetscCall(SetupRHS(dm, f, fLocal, f_rhs));
    
    // Get grid spacing
    const PetscReal hx = (x1 - x0) / nx;
    const PetscReal hy = (y1 - y0) / ny;
    
    // Solve with Gauss-Seidel
    PetscReal res_norm;
    PetscInt its;
    
    PetscCall(ComputeResidualNorm(dm, u, uLocal, f, fLocal, nx, ny, hx, hy, &res_norm));
    
    for (its = 1; its <= max_its; ++its) {
        PetscCall(GaussSeidelSweep(dm, u, uLocal, f, fLocal, nx, ny, hx, hy));
        PetscCall(ComputeResidualNorm(dm, u, uLocal, f, fLocal, nx, ny, hx, hy, &res_norm));
        
        if (res_norm <= tol) break;
    }
    
    // Cleanup local vectors
    PetscCall(DMRestoreLocalVector(dm, &uLocal));
    PetscCall(DMRestoreLocalVector(dm, &fLocal));
    
    // Return outputs
    *dm_out = dm;
    *u_out = u;
    *f_out = f;
    
    PetscFunctionReturn(0);
}

} // namespace PoissonDMStagSolver
