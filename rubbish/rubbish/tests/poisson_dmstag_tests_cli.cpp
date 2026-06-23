#include <petsc.h>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include "poisson/poisson_dmstag_mms.hpp"
#include "poisson/poisson_dmstag_metrics.hpp"
#include "poisson/poisson_dmstag_solver.hpp"

struct TestConfig {
    int nx, ny;
    double x0, x1, y0, y1;
    std::vector<std::string> pairs;
    std::string bc_west, bc_east, bc_south, bc_north;
    double tol_l2, tol_linf;
    int refine_levels;
    std::string output_path;
};

void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << " [options]\n"
              << "Options:\n"
              << "  --grid nx ny          Grid size (required)\n"
              << "  --domain x0 x1 y0 y1  Domain bounds (default: 0 1 0 1)\n"
              << "  --pair ID             MMS pair: poly2|sinpi|cospi (repeatable)\n"
              << "  --bc W E S N          Boundary conditions (default: Dirichlet x4)\n"
              << "  --tol l2 linf         Tolerances (default: 1e-10 1e-8)\n"
              << "  --refine N            Refinement levels for convergence test\n"
              << "  --output path         JSON summary output path\n";
}

bool parse_args(int argc, char** argv, TestConfig& config) {
    // Defaults
    config.nx = config.ny = -1;
    config.x0 = 0.0; config.x1 = 1.0;
    config.y0 = 0.0; config.y1 = 1.0;
    config.bc_west = config.bc_east = config.bc_south = config.bc_north = "Dirichlet";
    config.tol_l2 = 1e-10;
    config.tol_linf = 1e-8;
    config.refine_levels = 0;
    
    for (int i = 1; i < argc; ) {
        if (strcmp(argv[i], "--grid") == 0 && i + 2 < argc) {
            config.nx = std::atoi(argv[i+1]);
            config.ny = std::atoi(argv[i+2]);
            i += 3;
        } else if (strcmp(argv[i], "--domain") == 0 && i + 4 < argc) {
            config.x0 = std::atof(argv[i+1]);
            config.x1 = std::atof(argv[i+2]);
            config.y0 = std::atof(argv[i+3]);
            config.y1 = std::atof(argv[i+4]);
            i += 5;
        } else if (strcmp(argv[i], "--pair") == 0 && i + 1 < argc) {
            config.pairs.push_back(argv[i+1]);
            i += 2;
        } else if (strcmp(argv[i], "--bc") == 0 && i + 4 < argc) {
            config.bc_west = argv[i+1];
            config.bc_east = argv[i+2];
            config.bc_south = argv[i+3];
            config.bc_north = argv[i+4];
            i += 5;
        } else if (strcmp(argv[i], "--tol") == 0 && i + 2 < argc) {
            config.tol_l2 = std::atof(argv[i+1]);
            config.tol_linf = std::atof(argv[i+2]);
            i += 3;
        } else if (strcmp(argv[i], "--refine") == 0 && i + 1 < argc) {
            config.refine_levels = std::atoi(argv[i+1]);
            i += 2;
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            config.output_path = argv[i+1];
            i += 2;
        } else {
            std::cerr << "Unknown or incomplete option: " << argv[i] << "\n";
            return false;
        }
    }
    
    if (config.nx <= 0 || config.ny <= 0) {
        std::cerr << "Error: --grid nx ny is required\n";
        return false;
    }
    
    if (config.pairs.empty()) {
        std::cerr << "Error: At least one --pair must be specified\n";
        return false;
    }
    
    return true;
}

int main(int argc, char** argv) {
    PetscErrorCode ierr;
    ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
    
    TestConfig config;
    if (!parse_args(argc, argv, config)) {
        print_usage(argv[0]);
        PetscFinalize();
        return 1;
    }
    
    int total_tests = 0;
    int passed_tests = 0;
    
    std::cout << "=== Poisson DMStag Test Runner ===\n";
    std::cout << "Grid: " << config.nx << " x " << config.ny << "\n";
    std::cout << "Domain: [" << config.x0 << "," << config.x1 << "] x ["
              << config.y0 << "," << config.y1 << "]\n";
    std::cout << "Pairs: ";
    for (const auto& p : config.pairs) std::cout << p << " ";
    std::cout << "\n\n";
    
    for (const auto& pair_id : config.pairs) {
        total_tests++;
        std::cout << "--- Test: " << pair_id << " ---\n";
        
        try {
            const PoissonDMStagMMS::MMSPair* pair = PoissonDMStagMMS::get_pair(pair_id);
            
            DM dm;
            Vec u, f;
            
            // Solve the Poisson problem with DMStag
            ierr = PoissonDMStagSolver::run_poisson_case(
                &dm, &u, &f,
                config.nx, config.ny,
                config.x0, config.x1, config.y0, config.y1,
                pair->f_rhs, pair->u_exact,
                config.bc_west.c_str(), config.bc_east.c_str(),
                config.bc_south.c_str(), config.bc_north.c_str());
            CHKERRQ(ierr);
            
            // Compute errors
            double l2_error, linf_error, residual_norm;
            ierr = PoissonDMStagMetrics::compute_l2_error(dm, u, pair->u_exact, &l2_error); 
            CHKERRQ(ierr);
            ierr = PoissonDMStagMetrics::compute_linf_error(dm, u, pair->u_exact, &linf_error); 
            CHKERRQ(ierr);
            ierr = PoissonDMStagMetrics::compute_residual_norm(dm, u, f, config.nx, config.ny, &residual_norm); 
            CHKERRQ(ierr);
            
            std::cout << "  L2 error:   " << l2_error << "\n";
            std::cout << "  Linf error: " << linf_error << "\n";
            std::cout << "  Residual:   " << residual_norm << "\n";
            
            bool pass = (l2_error < config.tol_l2 && linf_error < config.tol_linf);
            if (pass) {
                passed_tests++;
                std::cout << "  Status: PASS\n";
            } else {
                std::cout << "  Status: FAIL\n";
            }
            
            // Cleanup
            ierr = VecDestroy(&u); CHKERRQ(ierr);
            ierr = VecDestroy(&f); CHKERRQ(ierr);
            ierr = DMDestroy(&dm); CHKERRQ(ierr);
            
        } catch (const std::exception& e) {
            std::cerr << "  ERROR: " << e.what() << "\n";
        }
    }
    
    std::cout << "\n=== Summary ===\n";
    std::cout << "Total: " << total_tests << " | Passed: " << passed_tests 
              << " | Failed: " << (total_tests - passed_tests) << "\n";
    
    ierr = PetscFinalize(); CHKERRQ(ierr);
    return (passed_tests == total_tests) ? 0 : 1;
}
