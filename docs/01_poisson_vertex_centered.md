# Vertex-Centered Poisson Equation Solver

## Problem Description

Solve the 2D Poisson equation on the unit square with homogeneous Dirichlet boundary conditions:

$$-\Delta u = f \quad \text{in } \Omega = [0,1] \times [0,1]$$
$$u = 0 \quad \text{on } \partial\Omega$$

### Manufactured Solution

To verify accuracy, we use the method of manufactured solutions (MMS):

$$u_{\text{exact}}(x,y) = \sin(\pi x) \sin(\pi y)$$

The corresponding source term is:

$$f(x,y) = -\Delta u_{\text{exact}} = 2\pi^2 \sin(\pi x) \sin(\pi y)$$

This solution satisfies $u = 0$ on all boundaries exactly.

---

## Numerical Method

### Spatial Discretization

- **Grid**: DMDA vertex-centered structured grid with $N_x \times N_y$ points
- **Grid spacing**: $h_x = \frac{1}{N_x-1}$, $h_y = \frac{1}{N_y-1}$
- **DOF location**: Vertices (grid points) — unknowns live at $(x_i, y_j) = (i h_x, j h_y)$
- **Stencil**: Standard 5-point finite difference

For an interior vertex $(i,j)$, the discrete Laplacian is:

$$-\Delta_h u_{i,j} = \frac{-u_{i+1,j} + 2u_{i,j} - u_{i-1,j}}{h_x^2} + \frac{-u_{i,j+1} + 2u_{i,j} - u_{i,j-1}}{h_y^2}$$

For $h_x = h_y = h$:

$$-\Delta_h u_{i,j} = \frac{4u_{i,j} - u_{i+1,j} - u_{i-1,j} - u_{i,j+1} - u_{i,j-1}}{h^2}$$

### Boundary Conditions

Dirichlet BC $u = 0$ is enforced by setting the identity row at boundary vertices:

$$u_{i,j} = 0 \quad \text{for } i = 0, i = N_x-1, j = 0, j = N_y-1$$

### Linear Solver

- **Solver**: GMRES (Generalized Minimal Residual)
- **Preconditioner**: Block Jacobi with ILU on each block (for MPI), or ILU (for serial)
- **Tolerance**: $\|r\|_2 < 10^{-12}$ for convergence tests

---

## Convergence Results

### Spatial Convergence (Second-Order Accuracy)

| $N_x \times N_y$ | $h$ | $\|u - u_{\text{exact}}\|_{L^2}$ | Rate | KSP Iterations |
|:---:|:---:|:---:|:---:|:---:|
| 16×16 | 0.06667 | 1.832×10⁻³ | — | 16 |
| 32×32 | 0.03226 | 4.281×10⁻⁴ | 2.10 | 29 |
| 64×64 | 0.01587 | 1.036×10⁻⁴ | 2.05 | 59 |
| 128×128 | 0.00787 | 2.550×10⁻⁵ | 2.02 | 143 |
| 256×256 | 0.00392 | 6.324×10⁻⁶ | 2.01 | 405 |

**Convergence rate** is computed as:

$$\text{rate} = \frac{\log(e_h / e_{h/2})}{\log(2)}$$

where $e_h$ is the $L^2$ error at grid spacing $h$.

**✅ Second-order convergence verified.** The rate asymptotically approaches 2.0.

### Error Norm Definition

The $L^2$ error is computed as the discrete $L^2$ norm:

$$\|u - u_{\text{exact}}\|_{L^2} = \sqrt{h_x h_y \sum_{i,j} \big(u_{i,j} - u_{\text{exact}}(x_i, y_j)\big)^2}$$

---

## Comparison: Vertex-Centered vs Cell-Centered

| Property | Vertex-Centered | Cell-Centered (DMStag) |
|----------|:---------------:|:----------------------:|
| **PETSc DM** | DMDA | DMStag |
| **DOF location** | Grid vertices | Cell centers |
| **Boundary nodes** | Explicit DOFs | Ghost cells |
| **BC enforcement** | Identity row in matrix | $u_{\text{ghost}} = -u_{\text{interior}}$ |
| **Grid points** | $N_x \times N_y$ | $N_x \times N_y$ cells |
| **Spacing** | $h = 1/(N-1)$ | $h = 1/N$ |
| **Stencil** | 5-point at vertices | 5-point at cell centers |
| **Accuracy** | $\mathcal{O}(h^2)$ | $\mathcal{O}(h^2)$ |

---

## Usage

### Basic Run

```bash
# Solve with default 64×64 grid
./poisson_vertex_2D

# Specify grid and compute error
./poisson_vertex_2D -nx 128 -ny 128 -poisson_check_error

# Run convergence test
./poisson_vertex_2D -nx 64 -ny 64 -poisson_check_error -convergence_test
```

### MPI Parallel Run

```bash
mpirun -np 4 ./poisson_vertex_2D -nx 128 -ny 128 -poisson_check_error \
    -ksp_type gmres -pc_type bjacobi -sub_pc_type ilu
```

### CTest

```bash
# Run the vertex-centered Poisson convergence test
ctest -R poisson_vertex_convergence -V
```

### Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `-nx <N>` | Grid points in x-direction | 64 |
| `-ny <N>` | Grid points in y-direction | 64 |
| `-poisson_check_error` | Compute and print L² error | off |
| `-convergence_test` | Output parseable CONVERGENCE line | off |
| `-ksp_type <type>` | KSP solver (gmres, cg, etc.) | gmres |
| `-pc_type <type>` | Preconditioner (ilu, bjacobi, etc.) | ilu |

---

## Algorithm

```
1. Create DMDA grid with Nx × Ny vertices on [0,1]²
2. Assemble 5-point Laplacian matrix A
   - Interior: -Δ_h stencil
   - Boundary: identity row (Dirichlet BC)
3. Fill RHS vector b: f(x_i, y_j) at interior, 0 at boundaries
4. Solve Au = b with KSP GMRES
5. Compute L² error: ||u - u_exact||
```

---

## Related Solvers

- **Cell-Centered Poisson**: `src/1_poisson/poisson_DMStag_2D.cpp` — DMStag element-centered
- **Staggered Poisson**: Planned — DMStag face-centered for mixed formulations
- **Heat Equation**: `src/2_heat/heat_DMStag_2D.cpp` — time-dependent, implicit Euler
- **Stokes Equation**: `src/stokes/stokes_DMStag_2D.cpp` — saddle-point system
