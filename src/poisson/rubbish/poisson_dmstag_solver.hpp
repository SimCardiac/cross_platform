#ifndef POISSON_DMSTAG_SOLVER_HPP
#define POISSON_DMSTAG_SOLVER_HPP

#include <petscdm.h>
#include <petscdmstag.h>
#include <petscvec.h>
#include <petsc.h>

namespace PoissonDMStagSolver {

// Run a complete Poisson test case with DMStag
// Returns: 0 on success, PETSc error code otherwise
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
    PetscReal tol = 1e-8,
    PetscInt max_its = 40000);

} // namespace PoissonDMStagSolver

#endif // POISSON_DMSTAG_SOLVER_HPP
