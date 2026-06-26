# Staggered Poisson Solver

Mixed first-order formulation on a **staggered** (MAC) grid. DOFs: $u$ (cell centers), $g_x$ (x-faces), $g_y$ (y-faces).

## 1. Governing Equation

First-order system equivalent to $-\Delta u = f$:

$$g_x = \frac{\partial u}{\partial x},\quad g_y = \frac{\partial u}{\partial y},\quad -\left(\frac{\partial g_x}{\partial x} + \frac{\partial g_y}{\partial y}\right) = f$$

Same 3 MMS as [vertex Poisson](01_poisson_vertex_centered.md).

## 2. Discretization

- **Grid**: $N_x \times N_y$ cells, $g_x$ at $x = i h$, $g_y$ at $y = j h$, $u$ at cell centers, $h = 1/N$
- **System**: Saddle-point matrix $\begin{bmatrix} I & -G \\ -D & 0 \end{bmatrix} \begin{bmatrix} g \\ u \end{bmatrix} = \begin{bmatrix} b_g \\ f \end{bmatrix}$
- **Dirichlet BC**: $g_x - \frac{2}{h}u_0 = -\frac{2}{h}g$ at boundary faces (one-sided gradient)
- **Neumann BC**: $g_x = \mp g$ at boundary faces (identity row, no $u$ coupling)
- **All-Neumann**: Pin $u(0,0)$ with identity row
- **Solver**: PCLU (direct solver — saddle-point needs it), $10^{-8}$ regularization on $u$-block

## 3. Convergence Results

`-mms sinpi -bc_type DDDD`:

```
$ ./poisson_staggered_2D -nx 16 -ny 16 -mms sinpi -bc_type DDDD -poisson_check_error
Stagg Poisson: N=16x16 h=0.0625 DOFs=800
||u-u_ex||=0.00160949

$ ./poisson_staggered_2D -nx 32 -ny 32 -mms sinpi -bc_type DDDD -poisson_check_error
Stagg Poisson: N=32x32 h=0.03125 DOFs=3136
||u-u_ex||=0.000401792

$ ./poisson_staggered_2D -nx 64 -ny 64 -mms sinpi -bc_type DDDD -poisson_check_error
Stagg Poisson: N=64x64 h=0.015625 DOFs=12416
||u-u_ex||=0.00010041

$ ./poisson_staggered_2D -nx 128 -ny 128 -mms sinpi -bc_type DDDD -poisson_check_error
Stagg Poisson: N=128x128 h=0.0078125 DOFs=49408
||u-u_ex||=2.51584e-05
```

| $N$ | $h$ | DOFs | $\|u-u_{\text{ex}}\|_{L^2}$ | Rate |
|:---:|:---:|:---:|:---:|:---:|
| 16 | 0.0625 | 800 | $1.609\times10^{-3}$ | — |
| 32 | 0.03125 | 3,136 | $4.018\times10^{-4}$ | 2.00 |
| 64 | 0.015625 | 12,416 | $1.004\times10^{-4}$ | 2.00 |
| 128 | 0.0078125 | 49,408 | $2.516\times10^{-5}$ | 2.00 |

All 18 BC combinations verified at second order. Results match cell-centered Poisson (algebraically equivalent formulations).

## 4. Usage

Same CLI. Default solver is PCLU; override with `-pc_type`.

## 5. Source Code

- **Solver**: `src/1_poisson/poisson_staggered_2D.cpp`
