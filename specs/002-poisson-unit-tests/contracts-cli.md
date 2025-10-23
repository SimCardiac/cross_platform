# CLI Contract: Poisson Test Runner

We will drive tests using a single binary (e.g., ex_poisson_center or a small wrapper) with CLI flags. Multiple (u_exact, f) pairs are supported in one run via repeatable flags.

## Flags

- --grid nx ny
  - Required. Example: --grid 32 32
- --domain x0 x1 y0 y1
  - Optional. Default: 0 1 0 1
- --pair ID
  - Repeatable. Selects a manufactured solution pair. Accepted: poly2 | sinpi | cospi
  - Example: --pair sinpi --pair poly2
- --bc west east south north
  - Optional. Default: Dirichlet Dirichlet Dirichlet Dirichlet
  - Each is Dirichlet|Neumann
- --tol l2 linf
  - Optional. Default: 1e-10 1e-8
- --refine N
  - Optional. If provided, perform N successive uniform refinements and report observed order
- --output path
  - Optional. Where to write summary JSON (per-pair residuals, errors, and pass/fail)

## Behavior

- For each --pair occurrence, run a full solve with specified grid and BCs
- Report per-pair metrics: residual L2, error L2/Linf, and status
- If --refine is set, repeat solves on 2x refinements and compute slope of log(error) vs log(h)
- Exit code: 0 if all pairs pass; non-zero if any fail
