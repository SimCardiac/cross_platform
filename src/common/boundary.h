#pragma once
#include <petscsys.h>
#include "mms.h"

// ============================================================================
// Uniform boundary condition specification
// ============================================================================

enum BCType { BC_DIRICHLET = 0, BC_NEUMANN = 1 };

typedef PetscScalar (*BCValueFunc)(PetscScalar coord);

struct BoundaryCondition {
    BCType     left, right, bottom, top;
    BCValueFunc g_left, g_right, g_bottom, g_top;
};

// Evaluate BC value (returns 0 if NULL)
static inline PetscScalar EvalBC(BCValueFunc g, PetscScalar coord) {
    return g ? g(coord) : 0.0;
}

// Convenience: all homogeneous Dirichlet
static inline BoundaryCondition HomogeneousDirichletBC() {
    BoundaryCondition bc;
    bc.left = bc.right = bc.bottom = bc.top = BC_DIRICHLET;
    bc.g_left = bc.g_right = bc.g_bottom = bc.g_top = NULL;
    return bc;
}

// ============================================================================
// Auto-generate BC value functions from a ManufacturedSolution
// Uses a static pointer — call SetBCMMS() before using the BC funcs.
// ============================================================================
static const ManufacturedSolution *bc_mms_ptr = NULL;
static inline void SetBCMMS(const ManufacturedSolution &m) { bc_mms_ptr = &m; }

// Dirichlet value wrappers
static inline PetscScalar _bc_dir_l(PetscScalar y) { return bc_mms_ptr->u(0, y); }
static inline PetscScalar _bc_dir_r(PetscScalar y) { return bc_mms_ptr->u(1, y); }
static inline PetscScalar _bc_dir_b(PetscScalar x) { return bc_mms_ptr->u(x, 0); }
static inline PetscScalar _bc_dir_t(PetscScalar x) { return bc_mms_ptr->u(x, 1); }

// Neumann value wrappers (∂u/∂n, with outward normal sign)
static inline PetscScalar _bc_neu_l(PetscScalar y) { return -bc_mms_ptr->ux(0, y); } // n=(-1,0)
static inline PetscScalar _bc_neu_r(PetscScalar y) { return  bc_mms_ptr->ux(1, y); } // n=(+1,0)
static inline PetscScalar _bc_neu_b(PetscScalar x) { return -bc_mms_ptr->uy(x, 0); } // n=(0,-1)
static inline PetscScalar _bc_neu_t(PetscScalar x) { return  bc_mms_ptr->uy(x, 1); } // n=(0,+1)

// Pick the right value function for a side + BC type
static inline BCValueFunc _bc_pick(BCType t, BCValueFunc d, BCValueFunc n) {
    return (t == BC_DIRICHLET) ? d : n;
}

// Create a BoundaryCondition from an MMS + per-side types
// The BC values are automatically derived from the MMS
static inline BoundaryCondition MakeBCFromMMS(const ManufacturedSolution &mms,
                                               BCType l, BCType r, BCType b, BCType t) {
    SetBCMMS(mms);
    BoundaryCondition bc;
    bc.left = l;   bc.right = r;   bc.bottom = b;   bc.top = t;
    bc.g_left   = _bc_pick(l, _bc_dir_l, _bc_neu_l);
    bc.g_right  = _bc_pick(r, _bc_dir_r, _bc_neu_r);
    bc.g_bottom = _bc_pick(b, _bc_dir_b, _bc_neu_b);
    bc.g_top    = _bc_pick(t, _bc_dir_t, _bc_neu_t);
    return bc;
}
