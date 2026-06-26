# Cell-Centered Heat Solver

Heat equation on a **cell-centered** (DMStag) grid. Implicit Euler with ghost-cell BC treatment.

## 1. Governing Equation

Same PDE and MMS as [vertex heat](04_heat_vertex_centered.md). $\alpha = 0.1$, $\Delta t = h^2$, $T = 0.05$.

## 2. Discretization

- **Grid**: $N_x \times N_y$ cells, $h = 1/N$
- **Matrix**: $I - \alpha\Delta t\Delta_h$ — ghost-cell stencil + identity
- **Dirichlet BC**: Ghost $u_{-1}=2g-u_0$ → diagonal $+\alpha\Delta t/h^2$, RHS $+2\alpha\Delta t\,g/h^2$
- **Neumann BC**: Ghost $u_{-1}=u_0+hg$ → diagonal $-\alpha\Delta t/h^2$, RHS $+\alpha\Delta t\,g/h$
- **Solver**: GMRES + ILU, tolerance $10^{-12}$

## 3. Convergence Results

`-mms sinpi -bc_type DDDD`:

```
$ ./heat_center_2D -nx 16 -ny 16 -mms sinpi -bc_type DDDD -heat_check_error
Cell Heat: N=16x16 h=0.0625 dt=0.00384615 steps=13 T=0.05
||u-u_ex||=0.000311352

$ ./heat_center_2D -nx 32 -ny 32 -mms sinpi -bc_type DDDD -heat_check_error
Cell Heat: N=32x32 h=0.03125 dt=0.000961538 steps=52 T=0.05
||u-u_ex||=7.82144e-05

$ ./heat_center_2D -nx 64 -ny 64 -mms sinpi -bc_type DDDD -heat_check_error
Cell Heat: N=64x64 h=0.015625 dt=0.000243902 steps=205 T=0.05
||u-u_ex||=1.97324e-05

$ ./heat_center_2D -nx 128 -ny 128 -mms sinpi -bc_type DDDD -heat_check_error
Cell Heat: N=128x128 h=0.0078125 dt=6.09756e-05 steps=820 T=0.05
||u-u_ex||=4.9346e-06
```

| $N$ | $h$ | $\Delta t$ | Steps | $\|u-u_{\text{ex}}\|_{L^2}$ | Rate |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 16 | 0.0625 | $3.85\times10^{-3}$ | 13 | $3.114\times10^{-4}$ | — |
| 32 | 0.03125 | $9.62\times10^{-4}$ | 52 | $7.821\times10^{-5}$ | 1.99 |
| 64 | 0.015625 | $2.44\times10^{-4}$ | 205 | $1.973\times10^{-5}$ | 1.99 |
| 128 | 0.0078125 | $6.10\times10^{-5}$ | 820 | $4.935\times10^{-6}$ | 2.00 |

All 18 BC combinations verified at second order.

## 4. Usage

Same CLI as [vertex heat](04_heat_vertex_centered.md).

## 5. Source Code

- **Solver**: `src/2_heat/heat_center_2D.cpp`
