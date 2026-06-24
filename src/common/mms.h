#pragma once
#include <petscsys.h>
#include <cmath>

// ============================================================================
// Manufactured Solution (MMS) — exact solution + source + boundary derivatives
// ============================================================================

struct ManufacturedSolution {
    const char *name;
    PetscScalar (*u)(PetscScalar x, PetscScalar y);      // exact solution (steady)
    PetscScalar (*f)(PetscScalar x, PetscScalar y);      // source: f = -Δu (Poisson)
    PetscScalar (*ux)(PetscScalar x, PetscScalar y);     // du/dx (for Neumann BC)
    PetscScalar (*uy)(PetscScalar x, PetscScalar y);     // du/dy (for Neumann BC)
    // Time-dependent variants (nullptr → use steady version)
    PetscScalar (*u_td)(PetscScalar x, PetscScalar y, PetscScalar t, PetscScalar alpha);
    PetscScalar (*f_td)(PetscScalar x, PetscScalar y, PetscScalar t, PetscScalar alpha);
};

// ============================================================================
// MMS 1: sinpi — u = sin(πx)sin(πy), f = 2π² sin(πx)sin(πy)
// Poisson: u=0 on boundary. Heat: u(x,y,t)=e^{-2π²αt} sin(πx)sin(πy), f=0.
// ============================================================================
inline PetscScalar sinpi_u(PetscScalar x, PetscScalar y) {
    return std::sin(M_PI * x) * std::sin(M_PI * y);
}
inline PetscScalar sinpi_f(PetscScalar x, PetscScalar y) {
    return 2.0 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
}
inline PetscScalar sinpi_ux(PetscScalar x, PetscScalar y) {
    return M_PI * std::cos(M_PI * x) * std::sin(M_PI * y);
}
inline PetscScalar sinpi_uy(PetscScalar x, PetscScalar y) {
    return M_PI * std::sin(M_PI * x) * std::cos(M_PI * y);
}
inline PetscScalar sinpi_u_heat(PetscScalar x, PetscScalar y, PetscScalar t, PetscScalar alpha) {
    return std::exp(-2.0 * M_PI * M_PI * alpha * t) * std::sin(M_PI * x) * std::sin(M_PI * y);
}
inline PetscScalar sinpi_f_heat(PetscScalar, PetscScalar, PetscScalar, PetscScalar) {
    return 0.0;  // f_heat = u_t - αΔu = 0
}
static const ManufacturedSolution MMS_SINPI = {"sinpi", sinpi_u, sinpi_f, sinpi_ux, sinpi_uy, sinpi_u_heat, sinpi_f_heat};

// ============================================================================
// MMS 2: poly2 — u = x² + y², f = -4
// Heat: steady, u_t=0, f_heat = -4α
// ============================================================================
inline PetscScalar poly2_u(PetscScalar x, PetscScalar y) {
    return x * x + y * y;
}
inline PetscScalar poly2_f(PetscScalar, PetscScalar) {
    return -4.0;
}
inline PetscScalar poly2_ux(PetscScalar x, PetscScalar y) {
    (void)y; return 2.0 * x;
}
inline PetscScalar poly2_uy(PetscScalar x, PetscScalar y) {
    (void)x; return 2.0 * y;
}
inline PetscScalar poly2_u_heat(PetscScalar x, PetscScalar y, PetscScalar, PetscScalar) {
    return x * x + y * y;  // steady
}
inline PetscScalar poly2_f_heat(PetscScalar, PetscScalar, PetscScalar, PetscScalar alpha) {
    return -4.0 * alpha;  // f_heat = u_t - αΔu = 0 - α*(-4) = 4α? No: -αΔu = -α*(4) = -4α
}
static const ManufacturedSolution MMS_POLY2 = {"poly2", poly2_u, poly2_f, poly2_ux, poly2_uy, poly2_u_heat, poly2_f_heat};

// ============================================================================
// MMS 3: cospi — u = cos(πx)cos(πy), f = 2π² cos(πx)cos(πy)
// ∂u/∂n = 0 on all boundaries (homogeneous Neumann compatible)
// Heat: u(x,y,t)=e^{-2π²αt} cos(πx)cos(πy), f=0.
// ============================================================================
inline PetscScalar cospi_u(PetscScalar x, PetscScalar y) {
    return std::cos(M_PI * x) * std::cos(M_PI * y);
}
inline PetscScalar cospi_f(PetscScalar x, PetscScalar y) {
    return 2.0 * M_PI * M_PI * std::cos(M_PI * x) * std::cos(M_PI * y);
}
inline PetscScalar cospi_ux(PetscScalar x, PetscScalar y) {
    return -M_PI * std::sin(M_PI * x) * std::cos(M_PI * y);
}
inline PetscScalar cospi_uy(PetscScalar x, PetscScalar y) {
    return -M_PI * std::cos(M_PI * x) * std::sin(M_PI * y);
}
inline PetscScalar cospi_u_heat(PetscScalar x, PetscScalar y, PetscScalar t, PetscScalar alpha) {
    return std::exp(-2.0 * M_PI * M_PI * alpha * t) * std::cos(M_PI * x) * std::cos(M_PI * y);
}
inline PetscScalar cospi_f_heat(PetscScalar, PetscScalar, PetscScalar, PetscScalar) {
    return 0.0;
}
static const ManufacturedSolution MMS_COSPI = {"cospi", cospi_u, cospi_f, cospi_ux, cospi_uy, cospi_u_heat, cospi_f_heat};

// ============================================================================
// All available MMS pairs
// ============================================================================
static const ManufacturedSolution *ALL_MMS[] = {&MMS_SINPI, &MMS_POLY2, &MMS_COSPI, NULL};
