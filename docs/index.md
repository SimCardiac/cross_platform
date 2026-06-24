# Cross-Platform PDE Solver Framework

PETSc-based structured grid PDE solvers with manufactured solution (MMS) verification — supporting **Poisson**, **Heat**, and **Stokes** equations across 3 grid types.

[![PETSc](https://img.shields.io/badge/PETSc-3.24-blue)](https://petsc.org/)
[![C++](https://img.shields.io/badge/C++-17-blue)](https://isocpp.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

---

## Quick Start

```bash
# Build all solvers
cd build && cmake .. && make -j8

# Run a solver with convergence check
./poisson_vertex_2D -nx 64 -ny 64 -mms sinpi -poisson_check_error
./heat_vertex_2D -nx 64 -ny 64 -mms cospi -bc_type NNNN -heat_check_error
./stokes_projection_2D -nx 32 -ny 32 -stokes_check_error
```

## CLI Reference

All solvers accept these common flags:

| Flag | Description | Default |
|------|-------------|---------|
| `-nx N` | Grid points in x | 64 |
| `-ny N` | Grid points in y | 64 |
| `-mms NAME` | Manufactured solution: `sinpi`, `poly2`, `cospi` | `sinpi` |
| `-bc_type LRTB` | BC types: D=Dirichlet, N=Neumann (e.g. `DNDN`) | auto-detected |
| `-{solver}_check_error` | Print L2 error at end | off |
| `-convergence_test` | Output machine-readable convergence data | off |

## Solvers

| # | Solver | Grid | Equation | BC Types |
|---|--------|------|----------|----------|
| 1 | `poisson_vertex_2D` | Vertex (DMDA) | $-\Delta u = f$ | 18/18 ✅ |
| 2 | `poisson_DMStag_2D` | Cell (DMStag) | $-\Delta u = f$ | 18/18 ✅ |
| 3 | `poisson_staggered_2D` | Staggered | Mixed 1st-order | 18/18 ✅ |
| 4 | `heat_vertex_2D` | Vertex (DMDA) | $u_t = \alpha\Delta u + f$ | 18/18 ✅ |
| 5 | `heat_center_2D` | Cell (DMStag) | $u_t = \alpha\Delta u + f$ | 18/18 ✅ |
| 6 | `heat_staggered_2D` | Staggered | Mixed 1st-order | 18/18 ✅ |
| 7 | `stokes_projection_2D` | Staggered | Stokes projection | Built-in |
| 8 | `stokes_monolithic_2D` | Staggered | Stokes monolithic | Built-in |

## Manufactured Solutions

| Name | $u(x,y)$ | $f = -\Delta u$ | Default BC |
|------|----------|----------------|------------|
| `sinpi` | $\sin(\pi x)\sin(\pi y)$ | $2\pi^2\sin(\pi x)\sin(\pi y)$ | DDDD |
| `poly2` | $x^2 + y^2$ | $-4$ | DDDD |
| `cospi` | $\cos(\pi x)\cos(\pi y)$ | $2\pi^2\cos(\pi x)\cos(\pi y)$ | NNNN |

All solvers verified with **second-order convergence** across all BC configurations.

## Documentation

- [Project Overview](project_overview.md)
- [Installation Guide](installation.md)
- [Boundary Conditions](poisson_neumann_bc.md)
- [Convergence Results](poisson_convergence_results.md)
