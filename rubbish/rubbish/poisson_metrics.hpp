#ifndef POISSON_METRICS_HPP
#define POISSON_METRICS_HPP

#include <petsc.h>

namespace PoissonMetrics {

// Compute L2 norm of error between numerical solution and exact function
PetscErrorCode compute_l2_error(DM dm, Vec u_numerical, 
                                 double (*u_exact)(double, double),
                                 double* l2_error);

// Compute L-infinity norm of error
PetscErrorCode compute_linf_error(DM dm, Vec u_numerical,
                                   double (*u_exact)(double, double),
                                   double* linf_error);

// Compute residual norm ||Au - f||
PetscErrorCode compute_residual_norm(Mat A, Vec u, Vec f, double* residual_norm);

// Compute grid spacing for convergence rate
void compute_grid_spacing(DM dm, double* hx, double* hy);

} // namespace PoissonMetrics

#endif // POISSON_METRICS_HPP
