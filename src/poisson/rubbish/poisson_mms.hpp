#ifndef POISSON_MMS_HPP
#define POISSON_MMS_HPP

#include <string>
#include <cmath>

namespace PoissonMMS {

// Manufactured solution pair interface
struct MMSPair {
    std::string id;
    double (*u_exact)(double x, double y);
    double (*f_rhs)(double x, double y);
    double (*du_dx)(double x, double y);  // For Neumann BCs
    double (*du_dy)(double x, double y);  // For Neumann BCs
};

// Available manufactured solutions
namespace Solutions {
    // Pair 1: u(x,y) = x^2 + y^2, f(x,y) = -4
    double poly2_u(double x, double y);
    double poly2_f(double x, double y);
    double poly2_du_dx(double x, double y);
    double poly2_du_dy(double x, double y);
    
    // Pair 2: u(x,y) = sin(πx) sin(πy), f(x,y) = 2π^2 sin(πx) sin(πy)
    double sinpi_u(double x, double y);
    double sinpi_f(double x, double y);
    double sinpi_du_dx(double x, double y);
    double sinpi_du_dy(double x, double y);
    
    // Pair 3: u(x,y) = cos(πx) cos(πy), f(x,y) = 2π^2 cos(πx) cos(πy)
    double cospi_u(double x, double y);
    double cospi_f(double x, double y);
    double cospi_du_dx(double x, double y);
    double cospi_du_dy(double x, double y);
}

// Get MMS pair by ID
const MMSPair* get_pair(const std::string& id);

} // namespace PoissonMMS

#endif // POISSON_MMS_HPP
