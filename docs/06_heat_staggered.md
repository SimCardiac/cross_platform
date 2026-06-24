# Staggered (Mixed) Heat Equation Solver (Item 6)

## 1. Problem Description

$$\frac{\partial u}{\partial t} = \alpha \Delta u, \quad u=0 \text{ on } \partial\Omega, \quad \alpha=0.1$$

**Mixed first-order formulation** — split into 3 fields ($u$ on elements, $g_x,g_y$ on faces):

$$\begin{cases}
g_x - \dfrac{\partial u}{\partial x} = 0 & \text{(left faces)} \\[8pt]
g_y - \dfrac{\partial u}{\partial y} = 0 & \text{(down faces)} \\[8pt]
\dfrac{\partial u}{\partial t} - \alpha\left(\dfrac{\partial g_x}{\partial x} + \dfrac{\partial g_y}{\partial y}\right) = 0 & \text{(elements)}
\end{cases}$$

Implicit Euler yields the monolithic system:

$$\begin{bmatrix} I & 0 & -G_x \\ 0 & I & -G_y \\ -\alpha\Delta t D_x & -\alpha\Delta t D_y & I \end{bmatrix}
\begin{bmatrix} g_x \\ g_y \\ u \end{bmatrix}^{n+1} =
\begin{bmatrix} 0 \\ 0 \\ u^n \end{bmatrix}$$

The $(3,3)$ block is $I$ (not $0$ as in Poisson), so ~2 GMRES iterations per step.

**Manufactured solution** (from `SINPI_SCALAR_2D`): $u = e^{-2\pi^2\alpha t}\sin(\pi x)\sin(\pi y)$, $f=0$.

## 2. Numerical Method

| Component | Detail |
|-----------|--------|
| **Grid** | DMStag: $N_xN_y$ elements + $(N_x+1)N_y$ x-faces + $N_x(N_y+1)$ y-faces |
| **Time** | Implicit Euler, $\Delta t = h^2$ |
| **Space** | Central differences on staggered grid (boundary: one-sided $h/2$) |
| **Solver** | GMRES, tol $10^{-12}$ |

## 3. Convergence Results

| $N$ | $\|u-u_{\text{exact}}\|_{L^2}$ | Rate |
|:---:|:---:|:---:|
| 16 | 3.114×10⁻⁴ | — |
| 32 | 7.821×10⁻⁵ | 1.99 |
| 64 | 1.973×10⁻⁵ | 1.99 |
| 128 | 4.935×10⁻⁶ | 2.00 |

$T_{\text{final}}=0.05$. **✅ Second-order confirmed.**

> **Reproduce:**
> ```bash
> for n in 16 32 64 128; do
>   ./heat_staggered_2D -nx $n -ny $n -T 0.05 -heat_check_error -convergence_test
> done
> ```

Results match cell-centered heat (item 5) — the mixed formulation is algebraically equivalent but uses a 3-field monolithic system instead of a single scalar Helmholtz solve.

## 4. Usage

```bash
./heat_staggered_2D -nx 64 -ny 64 -T 0.05 -heat_check_error
ctest -R heat_staggered_convergence -V
```

## 5. Source Code

`src/2_heat/heat_staggered_2D.cpp`
