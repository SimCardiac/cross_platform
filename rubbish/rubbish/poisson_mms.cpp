#include "poisson_mms.hpp"
#include <cmath>
#include <stdexcept>

namespace PoissonMMS {

constexpr double PI = 3.14159265358979323846;

namespace Solutions {
    // Pair 1: poly2
    double poly2_u(double x, double y) { return x*x + y*y; }
    double poly2_f(double x, double y) { return -4.0; }
    double poly2_du_dx(double x, double y) { return 2.0 * x; }
    double poly2_du_dy(double x, double y) { return 2.0 * y; }
    
    // Pair 2: sinpi
    double sinpi_u(double x, double y) { return std::sin(PI * x) * std::sin(PI * y); }
    double sinpi_f(double x, double y) { return 2.0 * PI * PI * std::sin(PI * x) * std::sin(PI * y); }
    double sinpi_du_dx(double x, double y) { return PI * std::cos(PI * x) * std::sin(PI * y); }
    double sinpi_du_dy(double x, double y) { return PI * std::sin(PI * x) * std::cos(PI * y); }
    
    // Pair 3: cospi
    double cospi_u(double x, double y) { return std::cos(PI * x) * std::cos(PI * y); }
    double cospi_f(double x, double y) { return 2.0 * PI * PI * std::cos(PI * x) * std::cos(PI * y); }
    double cospi_du_dx(double x, double y) { return -PI * std::sin(PI * x) * std::cos(PI * y); }
    double cospi_du_dy(double x, double y) { return -PI * std::cos(PI * x) * std::sin(PI * y); }
}

// Static pair definitions
static const MMSPair poly2_pair = {
    "poly2",
    Solutions::poly2_u,
    Solutions::poly2_f,
    Solutions::poly2_du_dx,
    Solutions::poly2_du_dy
};

static const MMSPair sinpi_pair = {
    "sinpi",
    Solutions::sinpi_u,
    Solutions::sinpi_f,
    Solutions::sinpi_du_dx,
    Solutions::sinpi_du_dy
};

static const MMSPair cospi_pair = {
    "cospi",
    Solutions::cospi_u,
    Solutions::cospi_f,
    Solutions::cospi_du_dx,
    Solutions::cospi_du_dy
};

const MMSPair* get_pair(const std::string& id) {
    if (id == "poly2") return &poly2_pair;
    if (id == "sinpi") return &sinpi_pair;
    if (id == "cospi") return &cospi_pair;
    throw std::runtime_error("Unknown MMS pair ID: " + id);
}

} // namespace PoissonMMS
