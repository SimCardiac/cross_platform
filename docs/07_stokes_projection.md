# 7. Stokes Projection Method (Staggered Grid)

## 1. Problem Description

Unsteady Stokes on $[0,1]^2$, homogeneous Dirichlet BC:

$$\frac{\partial \mathbf{u}}{\partial t} - \nu\Delta\mathbf{u} + \nabla p = 0,\quad \nabla\cdot\mathbf{u} = 0,\quad \nu=0.1$$

**Manufactured solution** (divergence-free, BC-compatible, $f_x=f_y=0$):

$$u = \sin(\pi x)\cos(\pi y)\,e^{-2\pi^2\nu t},\quad v = -\cos(\pi x)\sin(\pi y)\,e^{-2\pi^2\nu t},\quad p = 0$$

## 2. Numerical Method

| Component | Detail |
|-----------|--------|
| **Grid** | DMStag: u (LEFT faces), v (DOWN faces), p (ELEMENT centers) |
| **Time** | Implicit Euler, $\Delta t = h^2$ |
| **Algorithm** | Chorin-Temam projection (4 sub-steps per step) |

**Per time step:**
1. **Helmholtz u** — matrix-free Red-Black GS on left faces: $(I-\nu\Delta t\Delta_h)\mathbf{u}^*=\mathbf{u}^n$
2. **Helmholtz v** — same on down faces
3. **Pressure Poisson** — GS on elements with Neumann BC: $\Delta_h p = \frac{1}{\Delta t}\nabla\cdot\mathbf{u}^*$
4. **Correction** — $\mathbf{u}^{n+1} = \mathbf{u}^* - \Delta t\nabla_h p$

## 3. Convergence Results

$T_{\text{final}}=0.01$, $\Delta t = h^2$:

| $N\times N$ | $h$ | $\|u-u_{\text{exact}}\|_{L^2}$ | Rate |
|:---:|:---:|:---:|:---:|
| 16 | 0.0625 | 6.25×10⁻⁵ | — |
| 32 | 0.03125 | 1.64×10⁻⁵ | 1.93 |
| 64 | 0.015625 | 4.27×10⁻⁶ | 1.94 |
| 128 | 0.0078125 | 1.07×10⁻⁶ | 2.00 |

**✅ Rate → 2.00 asymptotically.**

> **Reproduce:**
> ```bash
> for n in 16 32 64 128; do
>   ./stokes_projection_2D -nx $n -ny $n -T 0.01 -stokes_check_error -convergence_test
> done
> ```

## 4. Usage

```bash
./stokes_projection_2D -nx 64 -ny 64 -T 0.01 -stokes_check_error
ctest -R stokes_projection_convergence -V
```

| Option | Default | Description |
|--------|---------|-------------|
| `-nx`, `-ny` | 32 | Grid cells |
| `-T` | 0.01 | Final time |
| `-stokes_check_error` | off | Compute L² error |

## 5. Source Code

`src/3_stokes/stokes_projection_2D.cpp`
