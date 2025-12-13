#ifndef POISSON_SOLVER_HPP
#define POISSON_SOLVER_HPP

#include <petsc.h>

// Wrapper for running a single Poisson test case
// Creates DM, assembles system, solves, returns solution and operators
PetscErrorCode run_poisson_case(
    DM* dm_out,           // Output: DM object
    Vec* u_out,           // Output: solution vector
    Mat* A_out,           // Output: system matrix (optional, can be NULL)
    Vec* f_out,           // Output: RHS vector
    int nx, int ny,       // Grid resolution
    double x0, double x1, double y0, double y1,  // Domain bounds
    double (*rhs_func)(double, double),           // RHS function
    double (*bc_func)(double, double),            // Boundary condition function
    const char* bc_west,  // "Dirichlet" or "Neumann"
    const char* bc_east,
    const char* bc_south,
    const char* bc_north,
    double ksp_rtol = 1e-10,  // KSP relative tolerance
    int ksp_max_it = 10000    // KSP max iterations
);

#endif // POISSON_SOLVER_HPP
