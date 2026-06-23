# Stokes Projection Method (Item 7)

## Problem

Unsteady Stokes equations with homogeneous Dirichlet BC:

$$\frac{\partial u}{\partial t} - \nu\Delta u + \frac{\partial p}{\partial x} = 0, \quad
\frac{\partial v}{\partial t} - \nu\Delta v + \frac{\partial p}{\partial y} = 0, \quad
\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0$$

## Manufactured Solution

Divergence-free, Dirichlet-compatible:

$$u = \sin(\pi x)\cos(\pi y)e^{-2\pi^2\nu t}, \quad
v = -\cos(\pi x)\sin(\pi y)e^{-2\pi^2\nu t}, \quad p = 0$$

$\nu=0.1$, $f_x=f_y=0$.

## Algorithm (Projection Method)

Per time step on staggered grid (DMStag: $p$ at elements, $u$ at left faces, $v$ at down faces):

1. **Helmholtz** (GS): $(I-\nu\Delta t\Delta_h)u^* = u^n$, $(I-\nu\Delta t\Delta_h)v^* = v^n$
2. **Pressure Poisson** (GS): $\Delta_h p = \frac{1}{\Delta t}\nabla\cdot\mathbf{u}^*$, Neumann BC $\partial p/\partial n=0$
3. **Correction**: $u^{n+1}=u^*-\Delta t\nabla_h p$, $v^{n+1}=v^*-\Delta t\nabla_h p$

## Convergence

| $N$ | $h$ | $\|u-u_{\text{exact}}\|_{L^2}$ | Rate |
|:---:|:---:|:---:|:---:|
| 16 | 0.0625 | 6.255×10⁻⁵ | — |
| 32 | 0.03125 | 1.643×10⁻⁵ | 1.93 |
| 64 | 0.015625 | 4.271×10⁻⁶ | 1.94 |
| 128 | 0.0078125 | 1.068×10⁻⁶ | 2.00 |

$T=0.01$, $\Delta t=h^2$. Rate approaches 2.0 asymptotically.

**✅ Second-order accuracy confirmed.**

## Usage

```bash
./stokes_projection_2D -nx 64 -ny 64 -T 0.01 -stokes_check_error
ctest -R stokes_convergence -V
```

## Source

`src/3_stokes/stokes_projection_2D.cpp`
