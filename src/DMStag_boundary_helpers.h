#pragma once

#include <petscdm.h>
#include <petscdmstag.h>

// Enforce Dirichlet ghosts for DMStag element-centered DOFs.
// For cell-centered DOFs located on element centers, and also DOFs located on
// element faces (LEFT/DOWN), this helper fills ghost DOFs using reflection:
// u_ghost = -u_interior which enforces zero value at the physical boundary.
//
// The function expects the DMStag local array pointer (PetscScalar ***arr)
// and will modify in-place the ghost bands. It supports both element DOFs
// and the adjacent face DOFs (LEFT, DOWN) layout produced by DMStag product
// coordinates.

// Safe neighbor getters that return the neighbor value or the reflected ghost
// value (i.e. -center) when the neighbor lies outside the physical domain.
// These operate on a DMStag local array returned by DMStagVecGetArray and
// expect the element slot index for the DOF to be used (for cell-centered DOFs
// this is typically the element slot returned by DMStagGetLocationSlot).

inline PetscScalar DMStag_GetLeft(PetscScalar ***arr, PetscInt ex, PetscInt ey, PetscInt slot, PetscInt Nx) {
	if (ex == 0) return -arr[ey][ex][slot];
	return arr[ey][ex-1][slot];
}

inline PetscScalar DMStag_GetRight(PetscScalar ***arr, PetscInt ex, PetscInt ey, PetscInt slot, PetscInt Nx) {
	if (ex == Nx-1) return -arr[ey][ex][slot];
	return arr[ey][ex+1][slot];
}

inline PetscScalar DMStag_GetDown(PetscScalar ***arr, PetscInt ex, PetscInt ey, PetscInt slot, PetscInt Ny) {
	if (ey == 0) return -arr[ey][ex][slot];
	return arr[ey-1][ex][slot];
}

inline PetscScalar DMStag_GetUp(PetscScalar ***arr, PetscInt ex, PetscInt ey, PetscInt slot, PetscInt Ny) {
	if (ey == Ny-1) return -arr[ey][ex][slot];
	return arr[ey+1][ex][slot];
}

// Backwards-compatible placeholder (no-op) kept so existing code can include header
// but the original EnforceDirichletGhosts is removed.
static inline PetscErrorCode EnforceDirichletGhosts(DM, PetscScalar ***, PetscInt, PetscInt) { return 0; }
