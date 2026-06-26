# Cell-Centered Poisson Solver

Steady-state Poisson on a **cell-centered** (DMStag) grid. Ghost-cell method for boundary conditions.

## 1. Governing Equation

$-\Delta u = f$ on $[0,1]^2$. Same 3 MMS as [vertex Poisson](01_poisson_vertex_centered.md).

## 2. Discretization

- **Grid**: $N_x \times N_y$ cells, centers at $((i+0.5)h, (j+0.5)h)$, $h = 1/N$
- **Stencil**: 5-point Laplacian at cell centers
- **Dirichlet BC**: Ghost cell $u_{-1} = 2g - u_0$ → diagonal $+1/h^2$, RHS $+2g/h^2$
- **Neumann BC**: Ghost cell $u_{-1} = u_0 + hg$ → diagonal $-1/h^2$, RHS $+g/h$
- **All-Neumann**: Pin cell $(0,0)$ to exact value
- **Solver**: GMRES + ILU, tolerance $10^{-10}$

## 3. Convergence Results

`-mms sinpi -bc_type DDDD`:

```
$ ./poisson_DMStag_2D -nx 16 -ny 16 -mms sinpi -bc_type DDDD -poisson_check_error
Cell Poisson: N=16x16 h=0.0625
||u-u_ex||=0.00160948

$ ./poisson_DMStag_2D -nx 32 -ny 32 -mms sinpi -bc_type DDDD -poisson_check_error
Cell Poisson: N=32x32 h=0.03125
||u-u_ex||=0.000401789

$ ./poisson_DMStag_2D -nx 64 -ny 64 -mms sinpi -bc_type DDDD -poisson_check_error
Cell Poisson: N=64x64 h=0.015625
||u-u_ex||=0.000100411

$ ./poisson_DMStag_2D -nx 128 -ny 128 -mms sinpi -bc_type DDDD -poisson_check_error
Cell Poisson: N=128x128 h=0.0078125
||u-u_ex||=2.51004e-05
```

| $N$ | $h$ | $\|u-u_{\text{ex}}\|_{L^2}$ | Rate |
|:---:|:---:|:---:|:---:|
| 16 | 0.0625 | $1.609\times10^{-3}$ | — |
| 32 | 0.03125 | $4.018\times10^{-4}$ | 2.00 |
| 64 | 0.015625 | $1.004\times10^{-4}$ | 2.00 |
| 128 | 0.0078125 | $2.510\times10^{-5}$ | 2.00 |

All 18 BC combinations verified at second order.

## 4. Usage

Same CLI as [vertex Poisson](01_poisson_vertex_centered.md).

## 5. Source Code

- **Solver**: `src/1_poisson/poisson_DMStag_2D.cpp`
