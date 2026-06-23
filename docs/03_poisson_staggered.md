# Staggered (Mixed) Poisson Equation Solver

## Problem Description

Solve the 2D Poisson equation using a **mixed first-order formulation** on a staggered DMStag grid:

$$-\Delta u = f \quad \text{in } \Omega = [0,1]^2, \qquad u = 0 \text{ on } \partial\Omega$$

The second-order equation is split into a first-order system with three unknown fields:

| Field | Location | Count |
|-------|----------|-------|
| $u$ | Element centers | $N_x \times N_y$ |
| $g_x = \partial u/\partial x$ | Left faces | $(N_x+1) \times N_y$ |
| $g_y = \partial u/\partial y$ | Down faces | $N_x \times (N_y+1)$ |

### First-Order System

$$\begin{cases}
g_x - \dfrac{\partial u}{\partial x} = 0 & \text{(on left faces)} \\[8pt]
g_y - \dfrac{\partial u}{\partial y} = 0 & \text{(on down faces)} \\[8pt]
-\left(\dfrac{\partial g_x}{\partial x} + \dfrac{\partial g_y}{\partial y}\right) = f & \text{(on elements)}
\end{cases}$$

In block matrix form:

$$\begin{bmatrix} I & 0 & -G_x \\ 0 & I & -G_y \\ -D_x & -D_y & 0 \end{bmatrix}
\begin{bmatrix} g_x \\ g_y \\ u \end{bmatrix} =
\begin{bmatrix} 0 \\ 0 \\ f \end{bmatrix}$$

The Schur complement $S = D_x G_x + D_y G_y$ recovers the standard 5-point Laplacian.

### Manufactured Solution

$$u_{\text{exact}} = \sin(\pi x)\sin(\pi y), \quad
g_x^{\text{exact}} = \pi\cos(\pi x)\sin(\pi y), \quad
g_y^{\text{exact}} = \pi\sin(\pi x)\cos(\pi y)$$

---

## Numerical Method

- **Spatial discretization**: Central differences on the staggered grid
  - Interior faces: $g_x = (u_R - u_L)/h$, distance $h$
  - Boundary faces: one-sided with distance $h/2$, e.g. $g_x(0) = 2u_0/h$
- **Solver**: GMRES (saddle-point system — $u$-block has zero diagonal)
- **Preconditioner**: None (identity blocks on $g_x,g_y$ make the system well-conditioned; typically converges in 2 iterations)

---

## Convergence

| $N$ | $h$ | $\|u - u_{\text{exact}}\|_{L^2}$ | Rate |
|:---:|:---:|:---:|:---:|
| 16 | 0.0625 | 1.609×10⁻³ | — |
| 32 | 0.03125 | 4.018×10⁻⁴ | 2.00 |
| 64 | 0.015625 | 1.004×10⁻⁴ | 2.00 |
| 128 | 0.0078125 | 2.510×10⁻⁵ | 2.00 |

**✅ Second-order accuracy confirmed** (rate = 2.00). Results are identical to cell-centered Poisson — the mixed formulation is algebraically equivalent but uses a different matrix structure and solver strategy.

---

## Usage

```bash
./poisson_staggered_2D -nx 64 -ny 64 -poisson_check_error

# MPI (no preconditioner needed)
mpirun -np 4 ./poisson_staggered_2D -nx 64 -ny 64 -poisson_check_error

# CTest
ctest -R poisson_staggered_convergence -V
```

| Option | Description | Default |
|--------|-------------|---------|
| `-nx`, `-ny` | Grid cells | 64 |
| `-poisson_check_error` | Compute L² error | off |
| `-convergence_test` | CONVERGENCE output | off |

### Why a Mixed Formulation?

This first-order system is the natural "staggered" discretization — it separates gradient and divergence operators onto different grid locations (faces vs. elements). It serves as a building block for the **Stokes projection method** solver.

## Source

`src/poisson/poisson_staggered_2D.cpp`
