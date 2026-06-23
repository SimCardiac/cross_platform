# Cell-Centered Heat Equation Solver (Item 5)

## Problem

$$\frac{\partial u}{\partial t} = \alpha \Delta u, \quad u=0 \text{ on } \partial\Omega, \quad \alpha=0.1$$

## Manufactured Solution

From `SINPI_SCALAR_2D` in `analytical/unsteady.h`:

$$u_{\text{exact}}(x,y,t) = e^{-2\pi^2\alpha t} \sin(\pi x) \sin(\pi y)$$

Satisfies $u=0$ on all boundaries exactly. Source term $f=0$.

## Numerical Method

| Component | Detail |
|-----------|--------|
| **Grid** | DMStag element-centered, $N\times N$ cells |
| **Time** | Implicit Euler, $\Delta t = h^2$ (coupled refinement) |
| **Space** | 5-point Helmholtz stencil |
| **BC** | Ghost-cell reflection ($u_{\text{ghost}}=-u_{\text{interior}}$) |
| **Solver** | Red-Black Gauss-Seidel, tol $10^{-8}$ |

## Convergence

| $N$ | $h$ | $\|u-u_{\text{exact}}\|_{L^2}$ | Rate |
|:---:|:---:|:---:|:---:|
| 16 | 0.0625 | 3.114×10⁻⁴ | — |
| 32 | 0.03125 | 7.822×10⁻⁵ | 1.99 |
| 64 | 0.015625 | 1.973×10⁻⁵ | 1.99 |
| 128 | 0.0078125 | 4.964×10⁻⁶ | 1.99 |

$T_{\text{final}}=0.05$. **✅ Second-order confirmed.**

## Usage

```bash
./heat_center_2D -nx 64 -ny 64 -T 0.05 -heat_check_error
mpirun -np 4 ./heat_center_2D -nx 64 -ny 64 -T 0.05 -heat_check_error
```

## Source

`src/2_heat/heat_center_2D.cpp`
