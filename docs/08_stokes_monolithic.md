# Stokes Monolithic Solver

Steady/unsteady Stokes via a **monolithic saddle-point** system. All DOFs $(u, v, p)$ solved simultaneously.

## 1. Governing Equation

Same Stokes equations as [projection method](07_stokes_projection.md). Implicit Euler saddle-point:

$$\begin{bmatrix} I - \nu\Delta t\Delta & 0 & \Delta t\,G_x \\ 0 & I - \nu\Delta t\Delta & \Delta t\,G_y \\ -D_x & -D_y & \varepsilon I \end{bmatrix}
\begin{bmatrix} u \\ v \\ p \end{bmatrix}^{n+1} =
\begin{bmatrix} u^n \\ v^n \\ 0 \end{bmatrix}$$

$\varepsilon = 10^{-8}$ on pressure diagonal (stabilization). $\Delta t = h^2$, $T_{\text{final}} = 0.01$.

## 2. Numerical Method

- **Grid**: Same staggered layout — $u$ (x-faces), $v$ (y-faces), $p$ (cell centers)
- **Pressure pinning**: $p(0,0) = 0$ removes constant nullspace
- **Matrix**: Manual assembly, DOFs $= (N_x+1)N_y + N_x(N_y+1) + N_x N_y$
- **Solver**: GMRES + ILU, tolerance $10^{-10}$

## 3. Convergence Results

```
$ ./stokes_monolithic_2D -nx 16 -ny 16 -stokes_check_error
Monolithic Stokes (pin): N=16x16 h=0.0625 dt=0.00333333 steps=3 DOFs=800
  done 3 steps
||u-u_ex||=0.000446155 ||v-v_ex||=0.00044602 ||p-p_ex||=0.00189402

$ ./stokes_monolithic_2D -nx 32 -ny 32 -stokes_check_error
Monolithic Stokes (pin): N=32x32 h=0.03125 dt=0.000909091 steps=11 DOFs=3136
  done 11 steps
||u-u_ex||=5.46558e-05 ||v-v_ex||=5.94786e-05 ||p-p_ex||=0.00099311

$ ./stokes_monolithic_2D -nx 64 -ny 64 -stokes_check_error
Monolithic Stokes (pin): N=64x64 h=0.015625 dt=0.000243902 steps=41 DOFs=12416
  done 41 steps
||u-u_ex||=7.00555e-05 ||v-v_ex||=6.75935e-05 ||p-p_ex||=0.00102508
```

| $N$ | DOFs | $\|u-u_{\text{ex}}\|$ | $\|v-v_{\text{ex}}\|$ | $\|p-p_{\text{ex}}\|$ | Rate (u) |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 16 | 800 | $4.462\times10^{-4}$ | $4.460\times10^{-4}$ | $1.894\times10^{-3}$ | — |
| 32 | 3,136 | $5.466\times10^{-5}$ | $5.948\times10^{-5}$ | $9.931\times10^{-4}$ | 3.03 |
| 64 | 12,416 | $7.006\times10^{-5}$ | $6.759\times10^{-5}$ | $1.025\times10^{-3}$ | −0.36 |

Velocity converges from N=16 to N=32 (~3rd order), then stagnates at N=64 due to ILU on the saddle-point system. Pressure error stagnates throughout. A proper saddle-point preconditioner (PCD, LSC, FieldSplit) is needed for optimal convergence.

## 4. Usage

Same CLI. Current implementation uses manual matrix assembly.

## 5. Source Code

- **Solver**: `src/3_stokes/stokes_monolithic_2D.cpp`
