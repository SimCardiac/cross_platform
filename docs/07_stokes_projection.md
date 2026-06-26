# Stokes Projection Solver

Unsteady Stokes flow via **Chorin's projection method** on a staggered (MAC) grid. Fractional-step: viscous → pressure Poisson → correction.

## 1. Governing Equation

$$\frac{\partial \mathbf{u}}{\partial t} = \nu\Delta\mathbf{u} - \nabla p + \mathbf{f}, \quad \nabla\cdot\mathbf{u} = 0$$

**Manufactured solution (built-in):**

$$u = e^{-2\pi^2\nu t}\sin(\pi x)\cos(\pi y),\quad v = -e^{-2\pi^2\nu t}\cos(\pi x)\sin(\pi y),\quad p = 0$$

$\nu = 0.1$, homogeneous Dirichlet on velocity, $T_{\text{final}} = 0.01$.

## 2. Numerical Method

- **Grid**: $u$ at x-faces, $v$ at y-faces, $p$ at cell centers, $h = 1/N$
- **Step 1 (Viscous)**: $(I - \nu\Delta t\Delta)\mathbf{u}^* = \mathbf{u}^n$ — Red-Black Gauss-Seidel
- **Step 2 (Projection)**: $\Delta p = \frac{1}{\Delta t}\nabla\cdot\mathbf{u}^*$
- **Step 3 (Correction)**: $\mathbf{u}^{n+1} = \mathbf{u}^* - \Delta t\,\nabla p$
- **Time**: $\Delta t = h^2$

## 3. Convergence Results

```
$ ./stokes_projection_2D -nx 16 -ny 16 -stokes_check_error
Stokes: N=16x16 h=0.0625 dt=0.00333333 steps=3
  step 3/3
||u-u_ex||=6.25464e-05
||v-v_ex||=6.25464e-05

$ ./stokes_projection_2D -nx 32 -ny 32 -stokes_check_error
Stokes: N=32x32 h=0.03125 dt=0.000909091 steps=11
  step 11/11
||u-u_ex||=1.6428e-05
||v-v_ex||=1.6428e-05

$ ./stokes_projection_2D -nx 64 -ny 64 -stokes_check_error
Stokes: N=64x64 h=0.015625 dt=0.000243902 steps=41
  step 41/41
||u-u_ex||=4.27063e-06
||v-v_ex||=4.27063e-06

$ ./stokes_projection_2D -nx 128 -ny 128 -stokes_check_error
Stokes: N=128x128 h=0.0078125 dt=6.09756e-05 steps=164
  step 100/164
  step 164/164
||u-u_ex||=1.068e-06
||v-v_ex||=1.068e-06
```

| $N$ | $h$ | Steps | $\|u-u_{\text{ex}}\|_{L^2}$ | $\|v-v_{\text{ex}}\|_{L^2}$ | Rate |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 16 | 0.0625 | 3 | $6.255\times10^{-5}$ | $6.255\times10^{-5}$ | — |
| 32 | 0.03125 | 11 | $1.643\times10^{-5}$ | $1.643\times10^{-5}$ | 1.93 |
| 64 | 0.015625 | 41 | $4.271\times10^{-6}$ | $4.271\times10^{-6}$ | 1.94 |
| 128 | 0.0078125 | 164 | $1.068\times10^{-6}$ | $1.068\times10^{-6}$ | 2.00 |

Second-order convergence verified for velocity. $u$ and $v$ errors identical due to symmetry of the MMS. Pressure identically zero ($p=0$ MMS).

## 4. Usage

| Flag | Description | Default |
|------|-------------|---------|
| `-nx N`, `-ny N` | Grid size | 32 |
| `-T T` | Final time | 0.01 |
| `-stokes_check_error` | Print L² error | off |

## 5. Source Code

- **Solver**: `src/3_stokes/stokes_projection_2D.cpp`
