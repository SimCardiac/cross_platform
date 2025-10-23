#include "poisson_solver.hpp"
#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <cstring>

// Assemble 5-point stencil Laplacian matrix on cell-centered grid
static PetscErrorCode assemble_laplacian(DM dm, Mat A, 
                                          int nx, int ny,
                                          double hx, double hy,
                                          const char* bc_west,
                                          const char* bc_east,
                                          const char* bc_south,
                                          const char* bc_north) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(dm, &info); CHKERRQ(ierr);
    
    MatStencil row, col[5];
    PetscScalar v[5];
    PetscInt ncols;
    
    double hx2_inv = 1.0 / (hx * hx);
    double hy2_inv = 1.0 / (hy * hy);
    double diag = 2.0 * hx2_inv + 2.0 * hy2_inv;
    
    for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
        for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
            row.i = i; row.j = j;
            ncols = 0;
            
            // Check if on Dirichlet boundary
            bool on_west = (i == 0) && (strcmp(bc_west, "Dirichlet") == 0);
            bool on_east = (i == nx - 1) && (strcmp(bc_east, "Dirichlet") == 0);
            bool on_south = (j == 0) && (strcmp(bc_south, "Dirichlet") == 0);
            bool on_north = (j == ny - 1) && (strcmp(bc_north, "Dirichlet") == 0);
            bool is_dirichlet = on_west || on_east || on_south || on_north;
            
            if (is_dirichlet) {
                // Identity row for Dirichlet BC
                col[ncols].i = i; col[ncols].j = j; v[ncols++] = 1.0;
            } else {
                // Center
                col[ncols].i = i; col[ncols].j = j; v[ncols++] = diag;
                
                // West neighbor
                if (i > 0) {
                    col[ncols].i = i-1; col[ncols].j = j; v[ncols++] = -hx2_inv;
                }
                
                // East neighbor
                if (i < nx - 1) {
                    col[ncols].i = i+1; col[ncols].j = j; v[ncols++] = -hx2_inv;
                }
                
                // South neighbor
                if (j > 0) {
                    col[ncols].i = i; col[ncols].j = j-1; v[ncols++] = -hy2_inv;
                }
                
                // North neighbor
                if (j < ny - 1) {
                    col[ncols].i = i; col[ncols].j = j+1; v[ncols++] = -hy2_inv;
                }
            }
            
            ierr = MatSetValuesStencil(A, 1, &row, ncols, col, v, INSERT_VALUES); CHKERRQ(ierr);
        }
    }
    
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode run_poisson_case(
    DM* dm_out, Vec* u_out, Mat* A_out, Vec* f_out,
    int nx, int ny,
    double x0, double x1, double y0, double y1,
    double (*rhs_func)(double, double),
    double (*bc_func)(double, double),
    const char* bc_west, const char* bc_east,
    const char* bc_south, const char* bc_north,
    double ksp_rtol, int ksp_max_it) {
    
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // Create DMDA (cell-centered, 2D)
    DM dm;
    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                        DMDA_STENCIL_STAR, nx, ny,
                        PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &dm); CHKERRQ(ierr);
    ierr = DMSetFromOptions(dm); CHKERRQ(ierr);
    ierr = DMSetUp(dm); CHKERRQ(ierr);
    
    // Set coordinates
    ierr = DMDASetUniformCoordinates(dm, x0, x1, y0, y1, 0.0, 0.0); CHKERRQ(ierr);
    
    double hx = (x1 - x0) / nx;
    double hy = (y1 - y0) / ny;
    
    // Create vectors
    Vec u, f;
    ierr = DMCreateGlobalVector(dm, &u); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(dm, &f); CHKERRQ(ierr);
    ierr = VecSet(u, 0.0); CHKERRQ(ierr);
    
    // Fill RHS and apply boundary conditions
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(dm, &info); CHKERRQ(ierr);
    
    PetscScalar **f_arr, **u_arr;
    ierr = DMDAVecGetArray(dm, f, &f_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(dm, u, &u_arr); CHKERRQ(ierr);
    
    for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
        for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
            double x = x0 + (i + 0.5) * hx;
            double y = y0 + (j + 0.5) * hy;
            
            // Check if on boundary
            bool on_west = (i == 0);
            bool on_east = (i == nx - 1);
            bool on_south = (j == 0);
            bool on_north = (j == ny - 1);
            bool on_boundary = on_west || on_east || on_south || on_north;
            
            if (on_boundary) {
                // For Dirichlet boundaries, set u directly to boundary value
                bool is_dirichlet = false;
                if (on_west && strcmp(bc_west, "Dirichlet") == 0) is_dirichlet = true;
                if (on_east && strcmp(bc_east, "Dirichlet") == 0) is_dirichlet = true;
                if (on_south && strcmp(bc_south, "Dirichlet") == 0) is_dirichlet = true;
                if (on_north && strcmp(bc_north, "Dirichlet") == 0) is_dirichlet = true;
                
                if (is_dirichlet) {
                    u_arr[j][i] = bc_func(x, y);
                    f_arr[j][i] = bc_func(x, y);  // RHS = boundary value for identity row
                } else {
                    f_arr[j][i] = rhs_func(x, y);
                }
            } else {
                f_arr[j][i] = rhs_func(x, y);
            }
        }
    }
    
    ierr = DMDAVecRestoreArray(dm, f, &f_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(dm, u, &u_arr); CHKERRQ(ierr);
    
    // Create matrix
    Mat A;
    ierr = DMCreateMatrix(dm, &A); CHKERRQ(ierr);
    ierr = assemble_laplacian(dm, A, nx, ny, hx, hy, bc_west, bc_east, bc_south, bc_north); CHKERRQ(ierr);
    
    // Setup KSP solver
    KSP ksp;
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPCG); CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, ksp_rtol, PETSC_DEFAULT, PETSC_DEFAULT, ksp_max_it); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
    
    // Handle nullspace for pure Neumann (T010)
    bool all_neumann = (strcmp(bc_west, "Neumann") == 0 &&
                        strcmp(bc_east, "Neumann") == 0 &&
                        strcmp(bc_south, "Neumann") == 0 &&
                        strcmp(bc_north, "Neumann") == 0);
    
    if (all_neumann) {
        MatNullSpace nullsp;
        ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nullsp); CHKERRQ(ierr);
        ierr = MatSetNullSpace(A, nullsp); CHKERRQ(ierr);
        ierr = MatNullSpaceRemove(nullsp, f); CHKERRQ(ierr);
        ierr = MatNullSpaceDestroy(&nullsp); CHKERRQ(ierr);
    }
    
    // Solve
    ierr = KSPSolve(ksp, f, u); CHKERRQ(ierr);
    
    // Output
    *dm_out = dm;
    *u_out = u;
    *f_out = f;
    if (A_out) {
        *A_out = A;
    } else {
        ierr = MatDestroy(&A); CHKERRQ(ierr);
    }
    
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}
