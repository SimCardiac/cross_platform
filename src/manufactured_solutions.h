#pragma once
#include <petscsys.h>
#include <cmath>

// ============================================================================
// Manufactured Solutions Library
// Contains exact solutions for various PDEs with different dimensions
// ============================================================================

// ============================================================================
// POISSON EQUATION: -∇²u = f
// ============================================================================
namespace POISSON {

// Function type aliases for 2D
namespace FUNC_2D {
  using EXACT = PetscScalar(*)(PetscScalar, PetscScalar);
  using RHS   = PetscScalar(*)(PetscScalar, PetscScalar);
}

// Function type aliases for 3D
namespace FUNC_3D {
  using EXACT = PetscScalar(*)(PetscScalar, PetscScalar, PetscScalar);
  using RHS   = PetscScalar(*)(PetscScalar, PetscScalar, PetscScalar);
}

// Solution 1: Trigonometric solution (homogeneous Dirichlet BC)
namespace TRIG_2D {
  static inline PetscScalar u_exact(PetscScalar x, PetscScalar y) {
    return std::sin(M_PI * x) * std::sin(M_PI * y);
  }
  
  static inline PetscScalar rhs_f(PetscScalar x, PetscScalar y) {
    return 2.0 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
  }
}

namespace TRIG_3D {
  static inline PetscScalar u_exact(PetscScalar x, PetscScalar y, PetscScalar z) {
    return std::sin(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * z);
  }
  
  static inline PetscScalar rhs_f(PetscScalar x, PetscScalar y, PetscScalar z) {
    return 3.0 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * z);
  }
}

// Solution 2: Polynomial solution
namespace POLY_2D {
  static inline PetscScalar u_exact(PetscScalar x, PetscScalar y) {
    return x * (1.0 - x) * y * (1.0 - y);
  }
  
  static inline PetscScalar rhs_f(PetscScalar x, PetscScalar y) {
    return 2.0 * (x * (1.0 - x) + y * (1.0 - y));
  }
}

namespace POLY_3D {
  static inline PetscScalar u_exact(PetscScalar x, PetscScalar y, PetscScalar z) {
    return x * (1.0 - x) * y * (1.0 - y) * z * (1.0 - z);
  }
  
  static inline PetscScalar rhs_f(PetscScalar x, PetscScalar y, PetscScalar z) {
    return 2.0 * (y * (1.0 - y) * z * (1.0 - z) +
                  x * (1.0 - x) * z * (1.0 - z) +
                  x * (1.0 - x) * y * (1.0 - y));
  }
}

} // namespace POISSON


// ============================================================================
// HEAT EQUATION: ∂u/∂t - α∇²u = f
// ============================================================================
namespace HEAT {

// Function type aliases for 2D
namespace FUNC_2D {
  using EXACT   = PetscScalar(*)(PetscScalar, PetscScalar, PetscScalar, PetscScalar);
  using INITIAL = PetscScalar(*)(PetscScalar, PetscScalar);
  using RHS     = PetscScalar(*)(PetscScalar, PetscScalar, PetscScalar);
}

// Function type aliases for 3D
namespace FUNC_3D {
  using EXACT   = PetscScalar(*)(PetscScalar, PetscScalar, PetscScalar, PetscScalar, PetscScalar);
  using INITIAL = PetscScalar(*)(PetscScalar, PetscScalar, PetscScalar);
  using RHS     = PetscScalar(*)(PetscScalar, PetscScalar, PetscScalar, PetscScalar);
}

// Solution 1: Exponential decay (homogeneous BC, zero source)
namespace DECAY_2D {
  static inline PetscScalar u_exact(PetscScalar x, PetscScalar y, PetscScalar t, PetscScalar alpha) {
    return std::exp(-2.0 * M_PI * M_PI * alpha * t) * std::sin(M_PI * x) * std::sin(M_PI * y);
  }
  
  static inline PetscScalar u_initial(PetscScalar x, PetscScalar y) {
    return std::sin(M_PI * x) * std::sin(M_PI * y);
  }
  
  static inline PetscScalar rhs_f(PetscScalar x, PetscScalar y, PetscScalar t) {
    return 0.0;
  }
}

namespace DECAY_3D {
  static inline PetscScalar u_exact(PetscScalar x, PetscScalar y, PetscScalar z, PetscScalar t, PetscScalar alpha) {
    return std::exp(-3.0 * M_PI * M_PI * alpha * t) * 
           std::sin(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * z);
  }
  
  static inline PetscScalar u_initial(PetscScalar x, PetscScalar y, PetscScalar z) {
    return std::sin(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * z);
  }
  
  static inline PetscScalar rhs_f(PetscScalar x, PetscScalar y, PetscScalar z, PetscScalar t) {
    return 0.0;
  }
}

// Solution 2: Manufactured solution with source term
namespace MMS_2D {
  static inline PetscScalar u_exact(PetscScalar x, PetscScalar y, PetscScalar t, PetscScalar alpha) {
    return std::sin(M_PI * x) * std::sin(M_PI * y) * std::cos(M_PI * t);
  }
  
  static inline PetscScalar u_initial(PetscScalar x, PetscScalar y) {
    return std::sin(M_PI * x) * std::sin(M_PI * y);
  }
  
  static inline PetscScalar rhs_f(PetscScalar x, PetscScalar y, PetscScalar t) {
    const PetscScalar alpha = 0.1; // Default value, should match solver
    return -M_PI * std::sin(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * t) +
           2.0 * M_PI * M_PI * alpha * std::sin(M_PI * x) * std::sin(M_PI * y) * std::cos(M_PI * t);
  }
}

namespace MMS_3D {
  static inline PetscScalar u_exact(PetscScalar x, PetscScalar y, PetscScalar z, PetscScalar t, PetscScalar alpha) {
    return std::sin(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * z) * std::cos(M_PI * t);
  }
  
  static inline PetscScalar u_initial(PetscScalar x, PetscScalar y, PetscScalar z) {
    return std::sin(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * z);
  }
  
  static inline PetscScalar rhs_f(PetscScalar x, PetscScalar y, PetscScalar z, PetscScalar t) {
    const PetscScalar alpha = 0.1; // Default value, should match solver
    return -M_PI * std::sin(M_PI * x) * std::sin(M_PI * y) * std::sin(M_PI * z) * std::sin(M_PI * t) +
           3.0 * M_PI * M_PI * alpha * std::sin(M_PI * x) * std::sin(M_PI * y) * 
           std::sin(M_PI * z) * std::cos(M_PI * t);
  }
}

// Solution 3: Non-zero Dirichlet BC (polynomial solution)
namespace NONZERO_DIRICHLET_2D {
  static inline PetscScalar u_exact(PetscScalar x, PetscScalar y, PetscScalar t, PetscScalar alpha) {
    // u(x,y,t) = (1 + 0.5*sin(pi*t)) * (x^2 + y^2)
    // This gives u = 0 at (0,0) but non-zero on other boundaries
    // For unit square: u = 1 + 0.5*sin(pi*t) on boundaries where x=1 or y=1
    return (1.0 + 0.5 * std::sin(M_PI * t)) * (x * x + y * y);
  }
  
  static inline PetscScalar u_initial(PetscScalar x, PetscScalar y) {
    return x * x + y * y;  // u(x,y,0) = x^2 + y^2
  }
  
  static inline PetscScalar rhs_f(PetscScalar x, PetscScalar y, PetscScalar t) {
    const PetscScalar alpha = 0.1; // Default value, should match solver
    // ∂u/∂t = 0.5*π*cos(πt) * (x^2 + y^2)
    // ∇²u = 4 * (1 + 0.5*sin(πt))
    // f = ∂u/∂t - α*∇²u
    return 0.5 * M_PI * std::cos(M_PI * t) * (x * x + y * y) - 
           4.0 * alpha * (1.0 + 0.5 * std::sin(M_PI * t));
  }
}

} // namespace HEAT


// ============================================================================
// STOKES EQUATION (time-dependent): ∂u/∂t - ν∇²u + ∇p = f, ∇·u = 0
// ============================================================================
namespace STOKES {

// Function type aliases for 2D
namespace FUNC_2D {
  using VELOCITY_EXACT = PetscScalar(*)(PetscScalar, PetscScalar, PetscScalar);
  using PRESSURE_EXACT = PetscScalar(*)(PetscScalar, PetscScalar, PetscScalar);
  using VELOCITY_INITIAL = PetscScalar(*)(PetscScalar, PetscScalar);
  using FORCE = PetscScalar(*)(PetscScalar, PetscScalar, PetscScalar);
}

// Function type aliases for 3D
namespace FUNC_3D {
  using VELOCITY_EXACT = PetscScalar(*)(PetscScalar, PetscScalar, PetscScalar, PetscScalar);
  using PRESSURE_EXACT = PetscScalar(*)(PetscScalar, PetscScalar, PetscScalar, PetscScalar);
  using VELOCITY_INITIAL = PetscScalar(*)(PetscScalar, PetscScalar, PetscScalar);
  using FORCE = PetscScalar(*)(PetscScalar, PetscScalar, PetscScalar, PetscScalar);
}

// Solution 1: Taylor-Green vortex (2D, time-dependent)
namespace TAYLOR_GREEN_2D {
  static inline PetscScalar u_exact(PetscScalar x, PetscScalar y, PetscScalar t) {
    const PetscScalar nu = 1.0; // Kinematic viscosity
    return -std::cos(M_PI * x) * std::sin(M_PI * y) * std::exp(-2.0 * M_PI * M_PI * nu * t);
  }
  
  static inline PetscScalar v_exact(PetscScalar x, PetscScalar y, PetscScalar t) {
    const PetscScalar nu = 1.0;
    return std::sin(M_PI * x) * std::cos(M_PI * y) * std::exp(-2.0 * M_PI * M_PI * nu * t);
  }
  
  static inline PetscScalar p_exact(PetscScalar x, PetscScalar y, PetscScalar t) {
    const PetscScalar nu = 1.0;
    return -0.25 * (std::cos(2.0 * M_PI * x) + std::cos(2.0 * M_PI * y)) * 
           std::exp(-4.0 * M_PI * M_PI * nu * t);
  }
  
  static inline PetscScalar u_initial(PetscScalar x, PetscScalar y) {
    return -std::cos(M_PI * x) * std::sin(M_PI * y);
  }
  
  static inline PetscScalar v_initial(PetscScalar x, PetscScalar y) {
    return std::sin(M_PI * x) * std::cos(M_PI * y);
  }
  
  static inline PetscScalar fx(PetscScalar x, PetscScalar y, PetscScalar t) {
    return 0.0; // Zero forcing for natural decay
  }
  
  static inline PetscScalar fy(PetscScalar x, PetscScalar y, PetscScalar t) {
    return 0.0;
  }
}

// Solution 2: Steady-state Taylor-Green (for steady Stokes)
namespace TAYLOR_GREEN_STEADY_2D {
  static inline PetscScalar u_exact(PetscScalar x, PetscScalar y, PetscScalar t = 0.0) {
    return -std::cos(M_PI * x) * std::sin(M_PI * y);
  }
  
  static inline PetscScalar v_exact(PetscScalar x, PetscScalar y, PetscScalar t = 0.0) {
    return std::sin(M_PI * x) * std::cos(M_PI * y);
  }
  
  static inline PetscScalar p_exact(PetscScalar x, PetscScalar y, PetscScalar t = 0.0) {
    return -0.25 * (std::cos(2.0 * M_PI * x) + std::cos(2.0 * M_PI * y));
  }
  
  static inline PetscScalar fx(PetscScalar x, PetscScalar y, PetscScalar t = 0.0) {
    const PetscScalar nu = 1.0;
    return 2.0 * M_PI * M_PI * nu * u_exact(x, y, t) + 0.5 * M_PI * std::sin(2.0 * M_PI * x);
  }
  
  static inline PetscScalar fy(PetscScalar x, PetscScalar y, PetscScalar t = 0.0) {
    const PetscScalar nu = 1.0;
    return 2.0 * M_PI * M_PI * nu * v_exact(x, y, t) + 0.5 * M_PI * std::sin(2.0 * M_PI * y);
  }
}

// Solution 3: 3D Taylor-Green vortex
namespace TAYLOR_GREEN_3D {
  static inline PetscScalar u_exact(PetscScalar x, PetscScalar y, PetscScalar z, PetscScalar t) {
    const PetscScalar nu = 1.0;
    return std::sin(M_PI * x) * std::cos(M_PI * y) * std::cos(M_PI * z) * 
           std::exp(-3.0 * M_PI * M_PI * nu * t);
  }
  
  static inline PetscScalar v_exact(PetscScalar x, PetscScalar y, PetscScalar z, PetscScalar t) {
    const PetscScalar nu = 1.0;
    return -std::cos(M_PI * x) * std::sin(M_PI * y) * std::cos(M_PI * z) * 
           std::exp(-3.0 * M_PI * M_PI * nu * t);
  }
  
  static inline PetscScalar w_exact(PetscScalar x, PetscScalar y, PetscScalar z, PetscScalar t) {
    return 0.0;
  }
  
  static inline PetscScalar p_exact(PetscScalar x, PetscScalar y, PetscScalar z, PetscScalar t) {
    const PetscScalar nu = 1.0;
    return (std::cos(2.0 * M_PI * x) + std::cos(2.0 * M_PI * y)) * 
           (std::cos(2.0 * M_PI * z) + 2.0) / 16.0 * 
           std::exp(-6.0 * M_PI * M_PI * nu * t);
  }
  
  static inline PetscScalar u_initial(PetscScalar x, PetscScalar y, PetscScalar z) {
    return std::sin(M_PI * x) * std::cos(M_PI * y) * std::cos(M_PI * z);
  }
  
  static inline PetscScalar v_initial(PetscScalar x, PetscScalar y, PetscScalar z) {
    return -std::cos(M_PI * x) * std::sin(M_PI * y) * std::cos(M_PI * z);
  }
  
  static inline PetscScalar w_initial(PetscScalar x, PetscScalar y, PetscScalar z) {
    return 0.0;
  }
  
  static inline PetscScalar fx(PetscScalar x, PetscScalar y, PetscScalar z, PetscScalar t) {
    return 0.0;
  }
  
  static inline PetscScalar fy(PetscScalar x, PetscScalar y, PetscScalar z, PetscScalar t) {
    return 0.0;
  }
  
  static inline PetscScalar fz(PetscScalar x, PetscScalar y, PetscScalar z, PetscScalar t) {
    return 0.0;
  }
}

// Solution 4: Kovasznay flow (2D steady-state with non-zero Reynolds number)
namespace KOVASZNAY_2D {
  static inline PetscScalar lambda(PetscScalar Re) {
    return Re / 2.0 - std::sqrt(Re * Re / 4.0 + 4.0 * M_PI * M_PI);
  }
  
  static inline PetscScalar u_exact(PetscScalar x, PetscScalar y, PetscScalar t = 0.0) {
    const PetscScalar Re = 40.0;
    const PetscScalar lam = lambda(Re);
    return 1.0 - std::exp(lam * x) * std::cos(2.0 * M_PI * y);
  }
  
  static inline PetscScalar v_exact(PetscScalar x, PetscScalar y, PetscScalar t = 0.0) {
    const PetscScalar Re = 40.0;
    const PetscScalar lam = lambda(Re);
    return lam / (2.0 * M_PI) * std::exp(lam * x) * std::sin(2.0 * M_PI * y);
  }
  
  static inline PetscScalar p_exact(PetscScalar x, PetscScalar y, PetscScalar t = 0.0) {
    const PetscScalar Re = 40.0;
    const PetscScalar lam = lambda(Re);
    return 0.5 * (1.0 - std::exp(2.0 * lam * x));
  }
  
  // Requires computing forcing term based on the nonlinear terms
  static inline PetscScalar fx(PetscScalar x, PetscScalar y, PetscScalar t = 0.0) {
    // For Kovasznay flow, f = 0 in the original formulation
    return 0.0;
  }
  
  static inline PetscScalar fy(PetscScalar x, PetscScalar y, PetscScalar t = 0.0) {
    return 0.0;
  }
}

} // namespace STOKES


// ============================================================================
// HELPER TEMPLATE FOR EASY ACCESS
// ============================================================================
template<typename Solution>
struct MMS {
  using Type = Solution;
};
