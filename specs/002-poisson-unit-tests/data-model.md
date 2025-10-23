# Data Model for Poisson Tests

## Entities

- TestCase
  - id: string
  - description: string
  - grid: { nx: int, ny: int }
  - domain: { x0: double, x1: double, y0: double, y1: double }
  - pair: string  # identifier for (u_exact, f)
  - bc: { west: string, east: string, south: string, north: string }  # e.g., Dirichlet/Neumann
  - tol: { l2: double, linf: double }

- Pair
  - id: string
  - u_exact(x,y): function
  - f(x,y): function
  - derivs: optional function pointers for boundary derivatives if needed

- Boundary
  - type: enum { Dirichlet, Neumann }
  - value(x,y): function  # for Dirichlet u_exact
  - flux(x,y): function   # for Neumann ∂u/∂n

## Notes

- Domain default: unit square [0,1] × [0,1]
- Grid default: powers of two for convenience; tests may use (16,16), (32,32), (64,64)
- Tolerances default: l2=1e-10, linf=1e-8; can be overridden per TestCase via CLI
- Pair identifiers: poly2, sinpi, cospi
