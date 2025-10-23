# Research: Poisson Equation Unit Tests

## Decisions

- Testing harness: Use CTest integrated with CMake for simple orchestration; individual tests may be CLI-driven or use a small C++ test target (final choice can be CTest-only unless gtest adds clear value)
- Manufactured solution pairs (u_exact, f):
  1) u(x,y) = x^2 + y^2 → f(x,y) = -4
  2) u(x,y) = sin(πx) sin(πy) → f(x,y) = 2π^2 sin(πx) sin(πy)
  3) u(x,y) = cos(πx) cos(πy) → f(x,y) = 2π^2 cos(πx) cos(πy)
- CLI-only configuration: All parameters (grid size, pair selection, BC types, tolerances) supplied via command line
- Build flow: `source ~/.bashrc` → out-of-source build via CMake + Make
- Tolerances and targets: residual L2 ≤ 1e-10; identity checks ≤ 1e-8 relative (≤ 1e-12 absolute if needed); refinement order within 10% of expected

## Rationale

- CTest is natively supported by CMake, lowering integration friction and making it easy to run on CI
- The selected manufactured solutions are standard, smooth, and exercise both polynomial and trigonometric behaviors
- CLI-only inputs enable reproducibility without config files, matching user constraints
- Sourcing ~/.bashrc ensures PETSc/MPI environment variables are available prior to configuration

## Alternatives Considered

- GoogleTest: richer assertions and test structure; rejected initially to keep integration minimal; can be added later if structure grows
- JSON config files: rejected per requirement for command-line only; could be added as optional override in the future
- Domain definitions via physical lengths (Lx, Ly) and non-uniform grids: deferred for later complexity; start with unit square and uniform grids
