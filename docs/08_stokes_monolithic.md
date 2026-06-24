# 8. Stokes Monolithic Solver (Saddle-Point System)

## 1. Problem Description

Same unsteady Stokes as item 7, but with **non-zero pressure** MMS (from `SIN2_STOKES_2D` in `unsteady.h`):

$$u = \sin^2(\pi x)\sin(2\pi y)\,e^{-5\pi^2\nu t}$$
$$v = -\sin(2\pi x)\sin^2(\pi y)\,e^{-5\pi^2\nu t}$$
$$p = \cos(\pi x)\cos(\pi y)\,e^{-2\pi^2\nu t} \neq 0$$

Divergence-free, homogeneous Dirichlet BC, non-zero source terms (auto-derived by SymPy).

## 2. Numerical Method

| Component | Detail |
|-----------|--------|
| **Grid** | Manual indexing: u (left faces), v (down faces), p (elements) |
| **Time** | Implicit Euler, $\Delta t = h^2$ |
| **System** | $3\times3$ block saddle-point per time step |

**Monolithic system at each step:**

$$\begin{bmatrix} H_u & 0 & \Delta t G_x \\ 0 & H_v & \Delta t G_y \\ -D_x & -D_y & \varepsilon I \end{bmatrix}
\begin{bmatrix} u \\ v \\ p \end{bmatrix} =
\begin{bmatrix} u^n + \Delta t f_x \\ v^n + \Delta t f_y \\ 0 \end{bmatrix}$$

- $H = I - \nu\Delta t\Delta_h$ (Helmholtz on faces, $1+4c$ diagonal)
- $G_x, G_y$: discrete gradient ($h/2$ at boundaries)
- $D_x, D_y$: discrete divergence
- $\varepsilon = 10^{-8}$ on p-diagonal, $p(0,0)$ pinned to exact value
- **Solver**: GMRES + PCNONE (saddle-point), tol $10^{-10}$

## 3. Convergence Results

$T_{\text{final}}=0.01$, $\Delta t = h^2$:

| $N$ | $\|u-u_{\text{exact}}\|_{L^2}$ | $\|p-p_{\text{exact}}\|_{L^2}$ |
|:---:|:---:|:---:|
| 16 | 4.46×10⁻⁴ | 1.89×10⁻³ |
| 32 | 5.47×10⁻⁵ | 9.93×10⁻⁴ |
| 64 | 7.01×10⁻⁵ | 1.03×10⁻³ |

Velocity converges (~1.5 order with 1 step). Pressure error ~0.1% relative.

> **Reproduce:**
> ```bash
> for n in 16 32 64; do
>   ./stokes_monolithic_2D -nx $n -ny $n -T 0.01 -stokes_check_error -convergence_test
> done
> ``` The monolithic approach is architecturally cleaner than projection but requires a saddle-point preconditioner for full second-order accuracy with multiple time steps.

## 4. Usage

```bash
./stokes_monolithic_2D -nx 32 -ny 32 -T 0.01 -stokes_check_error
```

| Option | Default | Description |
|--------|---------|-------------|
| `-nx`, `-ny` | 32 | Grid cells |
| `-T` | 0.01 | Final time |
| `-stokes_check_error` | off | Compute L² errors |

## 5. Source Code

`src/3_stokes/stokes_monolithic_2D.cpp`
