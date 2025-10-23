#ifndef POISSON_DMSTAG_METRICS_HPP
#define POISSON_DMSTAG_METRICS_HPP

#include <petscdm.h>
#include <petscdmstag.h>
#include <petscvec.h>
#include <petsc.h>
#include <cmath>

namespace PoissonDMStagMetrics {

// Compute L2 error: ||u - u_exact||_L2
PetscErrorCode compute_l2_error(
    DM dm,
    Vec u,
    PetscScalar (*u_exact_func)(PetscScalar, PetscScalar),
    PetscReal* l2_error);

// Compute L-infinity error: ||u - u_exact||_∞
PetscErrorCode compute_linf_error(
    DM dm,
    Vec u,
    PetscScalar (*u_exact_func)(PetscScalar, PetscScalar),
    PetscReal* linf_error);

// Compute residual norm: ||f - Au||_2
PetscErrorCode compute_residual_norm(
    DM dm,
    Vec u,
    Vec f,
    PetscInt Nx,
    PetscInt Ny,
    PetscReal* residual_norm);

// Get grid spacing
inline void get_grid_spacing(PetscInt Nx, PetscInt Ny,
                            PetscReal x0, PetscReal x1,
                            PetscReal y0, PetscReal y1,
                            PetscReal* hx, PetscReal* hy) {
    *hx = (x1 - x0) / Nx;
    *hy = (y1 - y0) / Ny;
}

} // namespace PoissonDMStagMetrics

#endif // POISSON_DMSTAG_METRICS_HPP
