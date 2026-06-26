# Vertex-Centered Heat Solver

Time-dependent heat equation on a **vertex-centered** (DMDA) grid. Implicit Euler with $\Delta t = h^2$.

## 1. Governing Equation

$$\frac{\partial u}{\partial t} = \alpha\Delta u + f, \quad \alpha = 0.1$$

Implicit Euler: $(I - \alpha\Delta t\Delta) u^{n+1} = u^n + \Delta t\,f^{n+1}$

**Manufactured solutions (time-dependent):**

| Name | $u(x,y,t)$ | $f$ |
|------|-----------|-----|
| `sinpi` | $e^{-2\pi^2\alpha t}\sin(\pi x)\sin(\pi y)$ | $0$ |
| `poly2` | $x^2 + y^2$ (steady) | $-4\alpha$ |
| `cospi` | $e^{-2\pi^2\alpha t}\cos(\pi x)\cos(\pi y)$ | $0$ |

## 2. Discretization

- **Grid**: $N_x \times N_y$ vertices, $h = 1/(N-1)$
- **Time**: $\Delta t = h^2$, $T_{\text{final}} = 0.05$
- **Matrix**: Helmholtz $I - \alpha\Delta t\Delta_h$ — 5-point stencil + identity
- **Dirichlet BC**: Identity row, RHS $= u_{\text{exact}}(x,y,t^{n+1})$
- **Neumann BC**: Same stencil as Poisson, RHS $= u^n + \Delta t\,f + 2\alpha\Delta t\,g/h$
- **Solver**: GMRES + ILU, tolerance $10^{-12}$

$I - \alpha\Delta t\Delta$ is positive definite — no pin needed for all-Neumann.

## 3. Convergence Results

`-mms sinpi -bc_type DDDD`, $T = 0.05$:

```
$ ./heat_vertex_2D -nx 16 -ny 16 -mms sinpi -bc_type DDDD -heat_check_error
Vertex Heat: N=16x16 h=0.0666667 dt=0.00416667 steps=12 T=0.05
||u-u_ex||=0.000344859

$ ./heat_vertex_2D -nx 32 -ny 32 -mms sinpi -bc_type DDDD -heat_check_error
Vertex Heat: N=32x32 h=0.0322581 dt=0.00102041 steps=49 T=0.05
||u-u_ex||=8.31501e-05

$ ./heat_vertex_2D -nx 64 -ny 64 -mms sinpi -bc_type DDDD -heat_check_error
Vertex Heat: N=64x64 h=0.015873 dt=0.000251256 steps=199 T=0.05
||u-u_ex||=2.03436e-05

$ ./heat_vertex_2D -nx 128 -ny 128 -mms sinpi -bc_type DDDD -heat_check_error
Vertex Heat: N=128x128 h=0.00787402 dt=6.19579e-05 steps=807 T=0.05
||u-u_ex||=5.01342e-06
```

| $N$ | $h$ | $\Delta t$ | Steps | $\|u-u_{\text{ex}}\|_{L^2}$ | Rate |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 16 | 0.0667 | $4.17\times10^{-3}$ | 12 | $3.449\times10^{-4}$ | — |
| 32 | 0.0323 | $1.02\times10^{-3}$ | 49 | $8.315\times10^{-5}$ | 2.05 |
| 64 | 0.0159 | $2.51\times10^{-4}$ | 199 | $2.034\times10^{-5}$ | 2.03 |
| 128 | 0.00787 | $6.20\times10^{-5}$ | 807 | $5.013\times10^{-6}$ | 2.02 |

All 18 BC combinations verified at second order. Temporal + spatial errors both $\mathcal{O}(h^2)$ under coupled refinement.

## 4. Usage

| Flag | Description | Default |
|------|-------------|---------|
| `-nx N`, `-ny N` | Grid points | 64 |
| `-T T` | Final time | 0.05 |
| `-mms NAME` | `sinpi`, `poly2`, `cospi` | `sinpi` |
| `-bc_type LRTB` | BC types | auto |
| `-heat_check_error` | Print L² error | off |

## 5. Source Code

- **Solver**: `src/2_heat/heat_vertex_2D.cpp`
