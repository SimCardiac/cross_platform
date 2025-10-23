# Source Code Organization

This directory contains the source code for the cross-platform PETSc solver project.

## Directory Structure

```
src/
├── poisson/                  # Poisson equation solver components
│   ├── poisson_mms.{hpp,cpp}        # Manufactured solution pairs
│   ├── poisson_metrics.{hpp,cpp}    # Error metrics (L2, Linf, residual)
│   ├── poisson_solver.{hpp,cpp}     # Solver wrapper (KSP-based)
│   └── tests/                       # Poisson-specific tests
│       └── poisson_tests_cli.cpp    # CLI test runner
│
├── io/                       # I/O utilities
│   ├── vti.{hpp,cpp}                # Parallel VTI writer
│   └── a.h
│
├── dmstag/                   # DMStag utilities
├── space_filling_curve/      # Space-filling curve implementations
├── examples/                 # Example programs
│   └── vti_example.cpp
│
├── wyh/                      # Specialized solvers
│   ├── lid_driven_cavity_flow/
│   └── oscillating_disk/
│
├── poisson_DMStag_2D.cpp     # Original DMStag Poisson solver
├── poisson_DM_2D.cpp         # Original DMDA Poisson solver
├── heat_DMStag_2D.cpp        # Heat equation solver
├── steady_stokes_DMStag_2D.cpp  # Stokes solver
├── DMStag_boundary_helpers.{cpp,h}
├── DMStag_2D_demo.cpp
└── main.cpp

## Module Organization

### Poisson Module (`src/poisson/`)

Organized testing infrastructure for Poisson equation validation:

- **MMS Module**: Provides manufactured solution pairs (poly2, sinpi, cospi)
- **Metrics Module**: Computes L2/Linf errors and residual norms
- **Solver Module**: Wraps PETSc KSP for solving Poisson problems
- **Tests**: CLI-based test runner with CTest integration

### Build Targets

| Target | Source | Description |
|--------|--------|-------------|
| `poisson_mms` | `poisson/poisson_mms.cpp` | Static library for MMS pairs |
| `poisson_metrics` | `poisson/poisson_metrics.cpp` | Static library for error metrics |
| `poisson_solver` | `poisson/poisson_solver.cpp` | Static library for solver wrapper |
| `poisson_tests_cli` | `poisson/tests/poisson_tests_cli.cpp` | Executable test runner |
| `ex_poisson_stagger` | `poisson_DMStag_2D.cpp` | Original DMStag solver |
| `ex_poisson_center` | `poisson_DM_2D.cpp` | Original DMDA solver |

## Design Principles

1. **Modularity**: Each physics/test component in its own subdirectory
2. **Separation**: Original solvers and test infrastructure separated
3. **Reusability**: Libraries can be used by multiple test executables
4. **Testability**: Test code isolated in `tests/` subdirectories

## Future Organization

As the project grows, consider:
- `src/heat/` for heat equation components
- `src/stokes/` for Stokes flow components
- `src/common/` for shared utilities
- Moving original solvers to `examples/` or `benchmarks/`
