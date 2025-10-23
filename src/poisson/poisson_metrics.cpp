#include "poisson_metrics.hpp"
#include <cmath>
#include <petscdm.h>
#include <petscdmda.h>

namespace PoissonMetrics {

PetscErrorCode compute_l2_error(DM dm, Vec u_numerical,
                                 double (*u_exact)(double, double),
                                 double* l2_error) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // Get DM info
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(dm, &info); CHKERRQ(ierr);
    
    // Get domain bounds (PETSc uses array syntax)
    PetscReal bmin[3], bmax[3];
    ierr = DMGetBoundingBox(dm, bmin, bmax); CHKERRQ(ierr);
    PetscReal xmin = bmin[0], xmax = bmax[0];
    PetscReal ymin = bmin[1], ymax = bmax[1];
    
    PetscReal hx = (xmax - xmin) / info.mx;
    PetscReal hy = (ymax - ymin) / info.my;
    
    // Access solution array
    PetscScalar **u_arr;
    ierr = DMDAVecGetArrayRead(dm, u_numerical, &u_arr); CHKERRQ(ierr);
    
    PetscReal local_sum = 0.0;
    for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
        for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
            PetscReal x = xmin + (i + 0.5) * hx;
            PetscReal y = ymin + (j + 0.5) * hy;
            PetscReal exact = u_exact(x, y);
            PetscReal diff = u_arr[j][i] - exact;
            local_sum += diff * diff;
        }
    }
    
    ierr = DMDAVecRestoreArrayRead(dm, u_numerical, &u_arr); CHKERRQ(ierr);
    
    // Global reduction
    PetscReal global_sum;
    ierr = MPIU_Allreduce(&local_sum, &global_sum, 1, MPIU_REAL, MPIU_SUM, PetscObjectComm((PetscObject)dm)); CHKERRQ(ierr);
    
    *l2_error = std::sqrt(global_sum * hx * hy);
    
    PetscFunctionReturn(0);
}

PetscErrorCode compute_linf_error(DM dm, Vec u_numerical,
                                   double (*u_exact)(double, double),
                                   double* linf_error) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(dm, &info); CHKERRQ(ierr);
    
    PetscReal bmin[3], bmax[3];
    ierr = DMGetBoundingBox(dm, bmin, bmax); CHKERRQ(ierr);
    PetscReal xmin = bmin[0], xmax = bmax[0];
    PetscReal ymin = bmin[1], ymax = bmax[1];
    
    PetscReal hx = (xmax - xmin) / info.mx;
    PetscReal hy = (ymax - ymin) / info.my;
    
    PetscScalar **u_arr;
    ierr = DMDAVecGetArrayRead(dm, u_numerical, &u_arr); CHKERRQ(ierr);
    
    PetscReal local_max = 0.0;
    for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
        for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
            PetscReal x = xmin + (i + 0.5) * hx;
            PetscReal y = ymin + (j + 0.5) * hy;
            PetscReal exact = u_exact(x, y);
            PetscReal diff = std::abs(u_arr[j][i] - exact);
            if (diff > local_max) local_max = diff;
        }
    }
    
    ierr = DMDAVecRestoreArrayRead(dm, u_numerical, &u_arr); CHKERRQ(ierr);
    
    PetscReal global_max;
    ierr = MPIU_Allreduce(&local_max, &global_max, 1, MPIU_REAL, MPIU_MAX, PetscObjectComm((PetscObject)dm)); CHKERRQ(ierr);
    
    *linf_error = global_max;
    
    PetscFunctionReturn(0);
}

PetscErrorCode compute_residual_norm(Mat A, Vec u, Vec f, double* residual_norm) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    Vec residual;
    ierr = VecDuplicate(f, &residual); CHKERRQ(ierr);
    
    // residual = A*u - f
    ierr = MatMult(A, u, residual); CHKERRQ(ierr);
    ierr = VecAXPY(residual, -1.0, f); CHKERRQ(ierr);
    
    PetscReal norm;
    ierr = VecNorm(residual, NORM_2, &norm); CHKERRQ(ierr);
    
    ierr = VecDestroy(&residual); CHKERRQ(ierr);
    
    *residual_norm = norm;
    
    PetscFunctionReturn(0);
}

void compute_grid_spacing(DM dm, double* hx, double* hy) {
    DMDALocalInfo info;
    DMDAGetLocalInfo(dm, &info);
    
    PetscReal bmin[3], bmax[3];
    DMGetBoundingBox(dm, bmin, bmax);
    
    *hx = (bmax[0] - bmin[0]) / info.mx;
    *hy = (bmax[1] - bmin[1]) / info.my;
}

} // namespace PoissonMetrics
