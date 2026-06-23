# Vertex-Centered Heat Equation Solver

## Problem Description

Solve the 2D heat (diffusion) equation with homogeneous Dirichlet BC:

$$\frac{\partial u}{\partial t} = \alpha \Delta u, \quad u = 0 \text{ on } \partial\Omega$$

- Domain: $\Omega = [0,1]^2$
- Diffusivity: $\alpha = 0.1$
- Initial condition: $u(x,y,0) = \sin(\pi x)\sin(\pi y)$

### Exact Solution

$$u_{\text{exact}}(x,y,t) = e^{-2\pi^2\alpha t} \sin(\pi x) \sin(\pi y)$$

The source term $f = 0$ — the homogeneous heat equation with this initial condition has a known closed-form solution.

---

## Numerical Method

### Time Discretization — Implicit Euler

$$\frac{u^{n+1} - u^n}{\Delta t} = \alpha \Delta_h u^{n+1}$$

Rearranged as a Helmholtz equation:

$$(I - \alpha\Delta t \Delta_h) u^{n+1} = u^n$$

### Spatial Discretization

- **Grid**: DMDA vertex-centered, $N_x \times N_y$ vertices
- **Spacing**: $h_x = 1/(N_x-1)$, $h_y = 1/(N_y-1)$
- **Stencil**: 5-point Helmholtz operator $(I - \alpha\Delta t \Delta_h)$
- **Boundary**: Identity rows at boundary vertices ($u = 0$)
- **Solver**: GMRES + ILU (serial) or Block Jacobi (MPI)

### Time Step Selection

For convergence tests, we use **coupled space-time refinement**:

$$\Delta t = h^2$$

This ensures the temporal error is $\mathcal{O}(h^2)$, matching the spatial accuracy.

---

## Convergence

| $N$ | $h$ | $\Delta t$ | Steps | $\|u - u_{\text{exact}}\|_{L^2}$ | Rate |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 16 | 0.0667 | 4.17×10⁻³ | 12 | 3.449×10⁻⁴ | — |
| 32 | 0.0323 | 1.02×10⁻³ | 49 | 8.315×10⁻⁵ | 2.05 |
| 64 | 0.0159 | 2.51×10⁻⁴ | 199 | 2.034×10⁻⁵ | 2.03 |
| 128 | 0.00787 | 6.20×10⁻⁵ | 807 | 5.013×10⁻⁶ | 2.02 |

$T_{\text{final}} = 0.05$, rate asymptotically approaches 2.0.

**✅ Second-order convergence verified** (implicit Euler is $\mathcal{O}(\Delta t)$, but with $\Delta t \propto h^2$ the overall scheme achieves $\mathcal{O}(h^2)$).

---

## Usage

```bash
# Default: 64×64, T=0.05, dt=h²
./heat_vertex_2D -heat_check_error

# Custom grid
./heat_vertex_2D -nx 32 -ny 32 -T 0.05 -heat_check_error

# Convergence test
./heat_vertex_2D -nx 64 -ny 64 -T 0.05 -heat_check_error -convergence_test
```

| Option | Description | Default |
|--------|-------------|---------|
| `-nx`, `-ny` | Grid vertices | 64 |
| `-T` | Final time | 0.05 |
| `-heat_check_error` | Compute L² error | off |
| `-convergence_test` | CONVERGENCE output | off |

## Source

`src/heat/heat_vertex_2D.cpp`
