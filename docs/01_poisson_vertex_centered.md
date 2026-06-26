# Vertex-Centered Poisson Solver

Steady-state Poisson equation on a **vertex-centered** (DMDA) grid. Supports Dirichlet and Neumann BCs with 3 manufactured solutions.

## 1. Governing Equation

$$-\Delta u = f \quad \text{in } \Omega = [0,1]^2$$

**Manufactured solutions (MMS):**

| Name | $u(x,y)$ | $f = -\Delta u$ | Default BC |
|------|----------|----------------|------------|
| `sinpi` | $\sin(\pi x)\sin(\pi y)$ | $2\pi^2\sin(\pi x)\sin(\pi y)$ | DDDD |
| `poly2` | $x^2 + y^2$ | $-4$ | DDDD |
| `cospi` | $\cos(\pi x)\cos(\pi y)$ | $2\pi^2\cos(\pi x)\cos(\pi y)$ | NNNN |

Boundary conditions set via `-bc_type LRTB` (e.g. `-bc_type DNDN` = Dirichlet left/right, Neumann bottom/top).

## 2. Discretization

- **Grid**: $N_x \times N_y$ vertices at $(i h_x, j h_y)$, $h_x = 1/(N_x-1)$, $h_y = 1/(N_y-1)$
- **Stencil**: 5-point Laplacian, interior diagonal $2(h_x^{-2}+h_y^{-2})$
- **Dirichlet BC**: Identity row at boundary vertices
- **Neumann BC**: Ghost-cell reflection with $2/h^2$ coefficients, RHS $= f + 2g/h$
- **All-Neumann**: Pin corner $(0,0)$ to exact value (removes constant nullspace)
- **Solver**: GMRES + ILU, tolerance $10^{-10}$

## 3. Convergence Results

`-mms sinpi -bc_type DDDD` — second-order verified:

```
$ ./poisson_vertex_2D -nx 16 -ny 16 -mms sinpi -bc_type DDDD -poisson_check_error
Vertex Poisson: N=16x16 h=0.0666667
||u-u_ex||=0.00183172

$ ./poisson_vertex_2D -nx 32 -ny 32 -mms sinpi -bc_type DDDD -poisson_check_error
Vertex Poisson: N=32x32 h=0.0322581
||u-u_ex||=0.000428142

$ ./poisson_vertex_2D -nx 64 -ny 64 -mms sinpi -bc_type DDDD -poisson_check_error
Vertex Poisson: N=64x64 h=0.015873
||u-u_ex||=0.000103624

$ ./poisson_vertex_2D -nx 128 -ny 128 -mms sinpi -bc_type DDDD -poisson_check_error
Vertex Poisson: N=128x128 h=0.00787402
||u-u_ex||=2.54973e-05
```

| $N$ | $h$ | $\|u-u_{\text{ex}}\|_{L^2}$ | Rate |
|:---:|:---:|:---:|:---:|
| 16 | 0.0667 | $1.832\times10^{-3}$ | — |
| 32 | 0.0323 | $4.281\times10^{-4}$ | 2.00 |
| 64 | 0.0159 | $1.036\times10^{-4}$ | 2.00 |
| 128 | 0.00787 | $2.550\times10^{-5}$ | 2.00 |

All 3 MMS × 6 BC configurations converge at second order. `poly2` achieves machine precision ($\sim10^{-11}$).

**Loop for batch testing:**
```bash
for n in 16 32 64 128; do
  ./poisson_vertex_2D -nx $n -ny $n -mms sinpi -poisson_check_error
done
```

**Other BC examples:**
```
$ ./poisson_vertex_2D -nx 64 -ny 64 -mms cospi -bc_type NNNN -poisson_check_error
Vertex Poisson: N=64x64 h=0.015873
||u-u_ex||=0.000236118

$ ./poisson_vertex_2D -nx 64 -ny 64 -mms sinpi -bc_type DNDN -poisson_check_error
Vertex Poisson: N=64x64 h=0.015873
||u-u_ex||=0.000182514
```

## 4. Usage

| Flag | Description | Default |
|------|-------------|---------|
| `-nx N`, `-ny N` | Grid points | 64 |
| `-mms NAME` | `sinpi`, `poly2`, `cospi` | `sinpi` |
| `-bc_type LRTB` | D=Dirichlet, N=Neumann | auto |
| `-poisson_check_error` | Print L² error | off |
| `-convergence_test` | Machine-readable output | off |

## 5. Source Code

- **Solver**: `src/1_poisson/poisson_vertex_2D.cpp`
- **BC/MMS**: `src/common/boundary.h`, `src/common/mms.h`
