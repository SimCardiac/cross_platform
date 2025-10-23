# Quickstart: Poisson DMStag Unit Tests

## Overview

The DMStag (Staggered Grid) version of Poisson equation tests uses element-centered (cell-centered) discretization with Gauss-Seidel iterative solver. Compared to DMDA version, DMStag provides better numerical accuracy and boundary condition handling.

## Build

1) Ensure PETSc and MPI are available in your environment
2) Source your shell init so environment variables are present:
   ```bash
   source ~/.bashrc
   ```
3) Configure and build (out-of-source):
   ```bash
   mkdir -p build && cd build
   cmake ..
   make -j
   ```

## Run

### Direct CLI execution

- Single pair on 32x32 grid:
  ```bash
  ./poisson_dmstag_tests_cli --grid 32 32 --pair sinpi
  ```

- Multiple pairs, custom tolerances:
  ```bash
  ./poisson_dmstag_tests_cli --grid 16 16 --pair sinpi --pair poly2 --tol 1e-10 1e-8
  ```

- Mixed boundary conditions:
  ```bash
  ./poisson_dmstag_tests_cli --grid 32 32 --pair cospi --bc Dirichlet Neumann Dirichlet Neumann
  ```

- With refinement for convergence testing:
  ```bash
  ./poisson_dmstag_tests_cli --grid 16 16 --pair sinpi --refine 2
  ```

### Test via CTest

From build directory:

- Run all Poisson tests (including DMStag):
  ```bash
  ctest -R poisson
  ```

- Run only DMStag tests:
  ```bash
  ctest -L dmstag
  ```

- Run with verbose output:
  ```bash
  ctest -L dmstag -V
  ```

- Run tests in parallel:
  ```bash
  ctest -L dmstag -j4
  ```

- Run specific test:
  ```bash
  ctest -R poisson_dmstag_smoke_poly2
  ```

## Expected Output

The CLI runner prints per-pair metrics:
```
=== Poisson DMStag Test Runner ===
Grid: 32 x 32
Domain: [0,1] x [0,1]
Pairs: sinpi 

--- Test: sinpi ---
  L2 error:   2.51e-05
  Linf error: 3.14e-05
  Residual:   6.11e-08
  Status: PASS

=== Summary ===
Total: 1 | Passed: 1 | Failed: 0
```

CTest output:
```
Test project /path/to/build
    Start 1: poisson_dmstag_smoke_poly2
1/3 Test #1: poisson_dmstag_smoke_poly2 .......   Passed    0.45 sec
    Start 2: poisson_dmstag_smoke_sinpi
2/3 Test #2: poisson_dmstag_smoke_sinpi .......   Passed    0.58 sec
    Start 3: poisson_dmstag_multiple_pairs
3/3 Test #3: poisson_dmstag_multiple_pairs ....   Passed    1.12 sec

100% tests passed, 0 tests failed out of 3
```

## DMStag vs DMDA

### DMStag Advantages
- **Better boundary condition handling**: Cell-centered discretization avoids special treatment at boundaries
- **Higher numerical accuracy**: Typically achieves better error orders for smooth solutions
- **Cleaner implementation**: Gauss-Seidel iteration without complex boundary corrections

### DMDA Advantages
- **KSP solver support**: Can use PETSc's advanced solvers (GMRES, CG, etc.)
- **More flexible BC options**: Supports more types of boundary conditions
- **Mature ecosystem**: More examples and documentation

### Recommendation
- For simple elliptic problems, **recommend DMStag**
- Use DMDA when advanced solvers or complex BCs needed
- For performance-critical applications, test both and choose faster

## Manufactured Solution (MMS) Pairs

### poly2 - Polynomial Solution
- Exact solution: `u(x,y) = x² + y²`
- RHS: `f(x,y) = -4`
- Feature: Constant RHS, suitable for testing basic operator correctness

### sinpi - Trigonometric Solution
- Exact solution: `u(x,y) = sin(πx)sin(πy)`
- RHS: `f(x,y) = 2π²sin(πx)sin(πy)`
- Feature: Zero boundary conditions, classic test case

### cospi - Cosine Solution
- Exact solution: `u(x,y) = cos(πx)cos(πy)`
- RHS: `f(x,y) = 2π²cos(πx)cos(πy)`
- Feature: Non-zero boundary values, tests BC implementation

## Notes

- All parameters must be provided via command line; no config files
- Exit code: 0 if all tests pass, non-zero if any fail
- Tests are deterministic; re-running produces identical results within tolerance
- DMStag uses Gauss-Seidel iteration, convergence rate sensitive to grid size
- For large grids (>128x128), consider increasing `max_its` or reducing `tol`

## Troubleshooting

### Test fails: L2 error too large
- Check grid is fine enough (at least 16x16)
- Verify tolerance settings are reasonable (default 1e-10 for L2)
- Confirm Gauss-Seidel converges (residual should be < 1e-8)

### Test timeout
- Reduce grid size or increase timeout limit
- DMStag's Gauss-Seidel can be slow on large grids
- Consider using MPI parallelism: `mpirun -np 4 ./poisson_dmstag_tests_cli ...`

### Residual doesn't converge
- Check RHS function is correct
- Verify boundary conditions match MMS pair
- For pure Neumann BCs, ensure compatibility condition (∫f = 0)
