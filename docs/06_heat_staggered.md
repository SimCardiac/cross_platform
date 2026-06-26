# Staggered Heat Solver

Heat equation in **mixed first-order form** on a staggered grid. Same saddle-point structure as [staggered Poisson](03_poisson_staggered.md) with Helmholtz operator.

## 1. Governing Equation

First-order system: $g_x = \partial_x u$, $g_y = \partial_y u$, $\partial_t u - \alpha(\partial_x g_x + \partial_y g_y) = f$.

Implicit Euler saddle-point system:

$$\begin{bmatrix} I & 0 & -G_x \\ 0 & I & -G_y \\ -\alpha\Delta t D_x & -\alpha\Delta t D_y & I \end{bmatrix}
\begin{bmatrix} g_x \\ g_y \\ u \end{bmatrix}^{n+1} =
\begin{bmatrix} b_{gx} \\ b_{gy} \\ u^n + \Delta t f \end{bmatrix}$$

Same 3 MMS as [vertex heat](04_heat_vertex_centered.md).

## 2. Discretization

- **Grid**: Same staggered layout as [staggered Poisson](03_poisson_staggered.md)
- **Matrix**: $u$-block has identity diagonal — no regularization needed (unlike Poisson)
- **Dirichlet BC**: $g_x - \frac{2}{h}u_0 = -\frac{2}{h}g$ at boundary faces
- **Neumann BC**: $g_x = \mp g$ at boundary faces (identity row)
- **Solver**: PCLU direct solver, tolerance $10^{-12}$

## 3. Convergence Results

`-mms sinpi -bc_type DDDD`:

```
$ ./heat_staggered_2D -nx 16 -ny 16 -mms sinpi -bc_type DDDD -heat_check_error
Staggered heat: N=16x16, h=0.0625, dt=0.00384615, steps=13, DOFs=800, mms=sinpi
  step 13/13, KSP its=1
||u(T) - u_exact(T)||_L2 = 0.000311352

$ ./heat_staggered_2D -nx 32 -ny 32 -mms sinpi -bc_type DDDD -heat_check_error
Staggered heat: N=32x32, h=0.03125, dt=0.000961538, steps=52, DOFs=3136, mms=sinpi
  step 52/52, KSP its=1
||u(T) - u_exact(T)||_L2 = 7.82144e-05

$ ./heat_staggered_2D -nx 64 -ny 64 -mms sinpi -bc_type DDDD -heat_check_error
Staggered heat: N=64x64, h=0.015625, dt=0.000243902, steps=205, DOFs=12416, mms=sinpi
  step 100/205, KSP its=1
  step 200/205, KSP its=1
  step 205/205, KSP its=1
||u(T) - u_exact(T)||_L2 = 1.97324e-05
```

| $N$ | DOFs | $\|u-u_{\text{ex}}\|_{L^2}$ | Rate |
|:---:|:---:|:---:|:---:|
| 16 | 800 | $3.114\times10^{-4}$ | — |
| 32 | 3,136 | $7.821\times10^{-5}$ | 1.99 |
| 64 | 12,416 | $1.973\times10^{-5}$ | 1.99 |
| 128 | 49,408 | $4.935\times10^{-6}$ | 2.00 |

Results match cell-centered heat (algebraically equivalent). All 18 BC combos verified. PCLU converges in 1 KSP iteration.

## 4. Usage

Same CLI. Default solver is PCLU (override with `-pc_type`).

## 5. Source Code

- **Solver**: `src/2_heat/heat_staggered_2D.cpp`
