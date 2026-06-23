# Cell-Centered Poisson Equation Solver (DMStag)

## Problem Description

Solve the 2D Poisson equation on the unit square with homogeneous Dirichlet boundary conditions:

$$-\Delta u = f \quad \text{in } \Omega = [0,1] \times [0,1], \qquad u = 0 \text{ on } \partial\Omega$$

### Manufactured Solution

$$u_{\text{exact}}(x,y) = \sin(\pi x) \sin(\pi y), \qquad f(x,y) = 2\pi^2 \sin(\pi x) \sin(\pi y)$$

---

## Numerical Method

- **Grid**: DMStag element-centered, $N_x \times N_y$ cells on $[0,1]^2$
- **Spacing**: $h_x = 1/N_x$, $h_y = 1/N_y$, cell centers at $x_i = (i+0.5)h_x$, $y_j = (j+0.5)h_y$
- **Stencil**: 5-point at cell centers
- **Boundary**: Ghost-cell reflection — $u_{\text{ghost}} = -u_{\text{interior}}$ enforces $u=0$ at the physical boundary
- **Solver**: Red-Black Gauss-Seidel (matrix-free), tolerance $10^{-8}$

$$\frac{4u_{i,j} - u_{i+1,j} - u_{i-1,j} - u_{i,j+1} - u_{i,j-1}}{h^2} = f_{i,j}$$

---

## Convergence

| $N$ | $h$ | $\|u - u_{\text{exact}}\|_{L^2}$ | Rate |
|:---:|:---:|:---:|:---:|
| 16 | 0.0625 | 1.609×10⁻³ | — |
| 32 | 0.03125 | 4.018×10⁻⁴ | 2.00 |
| 64 | 0.015625 | 1.004×10⁻⁴ | 2.00 |
| 128 | 0.0078125 | 2.510×10⁻⁵ | 2.00 |

**✅ Second-order accuracy confirmed** (rate = 2.00).

---

## Usage

```bash
./poisson_DMStag_2D -nx 64 -ny 64 -poisson_check_error
mpirun -np 4 ./poisson_DMStag_2D -nx 64 -ny 64 -convergence_test
```

## Source

`src/1_poisson/poisson_DMStag_2D.cpp`
