#ifndef POISSON_DMSTAG_MMS_HPP
#define POISSON_DMSTAG_MMS_HPP

#include <petsc.h>
#include <string>
#include <stdexcept>
#include <cmath>

namespace PoissonDMStagMMS {

// Manufactured solution pair
struct MMSPair {
    std::string id;
    std::string description;
    PetscScalar (*u_exact)(PetscScalar x, PetscScalar y);
    PetscScalar (*f_rhs)(PetscScalar x, PetscScalar y);
    PetscScalar (*du_dx)(PetscScalar x, PetscScalar y);  // for Neumann BC
    PetscScalar (*du_dy)(PetscScalar x, PetscScalar y);  // for Neumann BC
};

// Polynomial: u = x^2 + y^2, -Δu = -4
inline PetscScalar poly2_u(PetscScalar x, PetscScalar y) {
    return x * x + y * y;
}

inline PetscScalar poly2_f(PetscScalar x, PetscScalar y) {
    return -4.0;
}

inline PetscScalar poly2_du_dx(PetscScalar x, PetscScalar y) {
    return 2.0 * x;
}

inline PetscScalar poly2_du_dy(PetscScalar x, PetscScalar y) {
    return 2.0 * y;
}

// Trigonometric: u = sin(πx)sin(πy), -Δu = 2π²sin(πx)sin(πy)
inline PetscScalar sinpi_u(PetscScalar x, PetscScalar y) {
    return std::sin(M_PI * x) * std::sin(M_PI * y);
}

inline PetscScalar sinpi_f(PetscScalar x, PetscScalar y) {
    return 2.0 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
}

inline PetscScalar sinpi_du_dx(PetscScalar x, PetscScalar y) {
    return M_PI * std::cos(M_PI * x) * std::sin(M_PI * y);
}

inline PetscScalar sinpi_du_dy(PetscScalar x, PetscScalar y) {
    return M_PI * std::sin(M_PI * x) * std::cos(M_PI * y);
}

// Trigonometric: u = cos(πx)cos(πy), -Δu = 2π²cos(πx)cos(πy)
inline PetscScalar cospi_u(PetscScalar x, PetscScalar y) {
    return std::cos(M_PI * x) * std::cos(M_PI * y);
}

inline PetscScalar cospi_f(PetscScalar x, PetscScalar y) {
    return 2.0 * M_PI * M_PI * std::cos(M_PI * x) * std::cos(M_PI * y);
}

inline PetscScalar cospi_du_dx(PetscScalar x, PetscScalar y) {
    return -M_PI * std::sin(M_PI * x) * std::cos(M_PI * y);
}

inline PetscScalar cospi_du_dy(PetscScalar x, PetscScalar y) {
    return -M_PI * std::cos(M_PI * x) * std::sin(M_PI * y);
}

// Get MMS pair by ID
inline const MMSPair* get_pair(const std::string& id) {
    static const MMSPair pairs[] = {
        {"poly2", "Polynomial u=x²+y²", poly2_u, poly2_f, poly2_du_dx, poly2_du_dy},
        {"sinpi", "Trigonometric u=sin(πx)sin(πy)", sinpi_u, sinpi_f, sinpi_du_dx, sinpi_du_dy},
        {"cospi", "Trigonometric u=cos(πx)cos(πy)", cospi_u, cospi_f, cospi_du_dx, cospi_du_dy}
    };
    
    for (const auto& pair : pairs) {
        if (pair.id == id) return &pair;
    }
    
    throw std::runtime_error("Unknown MMS pair ID: " + id);
}

} // namespace PoissonDMStagMMS

#endif // POISSON_DMSTAG_MMS_HPP
