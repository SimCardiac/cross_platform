# Quickstart: Poisson Unit Tests

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
  ./poisson_tests_cli --grid 32 32 --pair sinpi
  ```

- Multiple pairs, custom tolerances:
  ```bash
  ./poisson_tests_cli --grid 16 16 --pair sinpi --pair poly2 --tol 1e-10 1e-8
  ```

- Mixed boundary conditions:
  ```bash
  ./poisson_tests_cli --grid 32 32 --pair cospi --bc Dirichlet Neumann Dirichlet Neumann
  ```

- With refinement for convergence testing:
  ```bash
  ./poisson_tests_cli --grid 16 16 --pair sinpi --refine 2
  ```

### Test via CTest

From build directory:

- Run all Poisson tests:
  ```bash
  ctest -R poisson
  ```

- Run with verbose output:
  ```bash
  ctest -R poisson -V
  ```

- Run tests in parallel:
  ```bash
  ctest -R poisson -j4
  ```

- Run specific test:
  ```bash
  ctest -R poisson_smoke_poly2
  ```

## Expected Output

The CLI runner prints per-pair metrics:
```
=== Poisson Test Runner ===
Grid: 32 x 32
Domain: [0,1] x [0,1]
Pairs: sinpi 

--- Test: sinpi ---
  L2 error:   1.23e-05
  Linf error: 4.56e-06
  Residual:   2.34e-11
  Status: PASS

=== Summary ===
Total: 1 | Passed: 1 | Failed: 0
```

CTest output:
```
Test project /path/to/build
    Start 1: poisson_smoke_poly2
1/4 Test #1: poisson_smoke_poly2 ..............   Passed    0.52 sec
    Start 2: poisson_smoke_sinpi
2/4 Test #2: poisson_smoke_sinpi ..............   Passed    0.61 sec
...
100% tests passed, 0 tests failed out of 4
```

## Notes

- All parameters must be provided via command line; no config files
- Exit code: 0 if all tests pass, non-zero if any fail
- Tests are deterministic; re-running produces identical results within tolerance
