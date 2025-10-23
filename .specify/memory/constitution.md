<!--
Sync Impact Report
- Version change: 1.0.0 → 1.1.0 (MINOR: added phased roadmap, expanded equation coverage and accuracy requirements, clarified coupling method)
- Modified/Added Principles:
	- Added V. Phased Delivery & Equation Coverage
- Added sections:
	- Roadmap & Phasing
- Removed sections:
	- None
- Templates requiring updates:
	- .specify/templates/spec-template.md ✅ (alignment verified; no structural change required)
	- .specify/templates/plan-template.md ✅ (Constitution Check gate remains valid)
	- .specify/templates/tasks-template.md ✅ (phasing already supported via phases; no change required)
	- .specify/templates/commands/* ⚠ pending (no files found; nothing to update)
- Follow-up TODOs:
	- TODO(RATIFICATION_DATE_HISTORY): If historical ratification differs from 2025-10-23, update accordingly.
-->

# Incompressible Navier-Stokes Solver Constitution

## Core Principles

### I. PETSc API Compliance
- All numerical operations MUST use PETSc data structures (Vec, Mat, DM, KSP, SNES)
- Direct manipulation of underlying arrays ONLY when PETSc provides no alternative
- Proper PETSc object lifecycle: Create → Configure → Use → Destroy
- Error handling via PETSc macros (CHKERRQ, PetscCall) is MANDATORY
- All PETSc objects must be properly destroyed to prevent memory leaks

### II. Modern C++ Standards
- Minimum C++17 standard; prefer C++20 features when available
- Use RAII for resource management (PETSc object wrappers encouraged)
- Prefer `std::unique_ptr` and `std::shared_ptr` for non-PETSc resources
- Use `constexpr` for compile-time computations when possible
- STL containers and algorithms over raw arrays/loops where performance permits

### III. Template Programming Strategy
- Use templates for dimension-agnostic code (2D/3D solvers from same codebase)
- Template parameters for: spatial dimensions, discretization schemes, boundary conditions
- Compile-time polymorphism via CRTP (Curiously Recurring Template Pattern) preferred over runtime polymorphism
- Template specialization for performance-critical sections ONLY
- Concepts (C++20) or SFINAE to constrain template parameters

### IV. Numerical Correctness
- Divergence-free velocity field enforcement (incompressibility constraint)
- Timestep stability: CFL condition checks for explicit schemes
- Validation against analytical solutions for test cases
- Convergence testing: spatial and temporal; equation coverage MUST include Poisson, Heat, and Incompressible flows

### V. Phased Delivery & Equation Coverage
- Delivery MUST follow phases: Phase 1 Poisson → Phase 2 Heat → Phase 3 Incompressible.
- Progression to a phase REQUIRES all prior phase acceptance criteria to pass.
- Pressure–velocity coupling MUST be via a projection method with both first- and second-order temporal accuracy variants.
- Each phase MUST include unit tests, verification (e.g., MMS), and documented exit criteria.

## Technical Requirements

### Discretization Standards
- Spatial: Finite Difference
- Temporal: Implicit Backward Euler; Crank–Nicolson. The convective term MAY be explicit (upwind, PPM, or WENO scheme).
- Pressure–velocity coupling: Projection method (e.g., Chorin first-order and a second-order variant such as incremental pressure-correction).
- Staggered grid (DMStag) preferred for incompressible flow

### Solver Configuration
- Linear solvers: self-coded matrix-free solver; projection method as a preconditioner is allowed
- Solver tolerances and max iterations configurable via command line
- Convergence monitoring with PETSc's built-in monitors

### Boundary Conditions
- Support: Dirichlet, Neumann
- Template-based BC application for extensibility
- BC specification via runtime configuration (PETSc options or config files)
## Roadmap & Phasing

- Phase 1 — Poisson:
	- Discrete operator identities, boundary-condition verification, residual and error convergence on refined grids.
	- Exit criteria: residual norms below target tolerance; observed convergence rates within defined bounds; docs updated.

- Phase 2 — Heat:
	- Temporal accuracy verification for first-order (Backward Euler) and second-order (Crank–Nicolson) using manufactured solutions.
	- Exit criteria: measured temporal orders ≈1 and ≈2 respectively across time-step halvings; stability checks documented.

- Phase 3 — Incompressible:
	- Projection method tests (first- and second-order variants), divergence reduction, boundary handling, pressure nullspace robustness.
	- Exit criteria: post-projection divergence below tolerance; temporal accuracy verified; comprehensive documentation complete.


### Code Organization
- Separate modules: Grid management, Discretization, Time integration, Solvers, I/O
- Header-only templates in `include/`, implementations in `src/`
- Example programs demonstrating: lid-driven cavity, channel flow, Taylor-Green vortex

## Quality Assurance

### Testing Requirements
- Unit tests for: grid operations, discretization stencils, BC application
- Integration tests for: pressure Poisson solve, momentum solve, full timestep
- Verification tests: method of manufactured solutions (MMS)
- Phase-specific tests:
	- Poisson: operator consistency and solver convergence on refined grids
	- Heat: first- and second-order temporal accuracy on manufactured solutions
	- Incompressible: projection method divergence reduction and temporal accuracy
- Performance regression tests: track solve times, memory usage

### Documentation Standards
- Doxygen comments for all public APIs
- Mathematical formulation documented in separate markdown files
- README with: build instructions, dependencies, quick start examples
- Inline comments explaining non-obvious numerical techniques
- Comprehensive docs MUST include:
	- Covered equations (Poisson, Heat, Incompressible) and discretization/time-integration choices
	- Projection method variants (1st/2nd order), assumptions, and limitations
	- Phase exit criteria, test tolerances, and how to extend tests and benchmarks

## Governance

- Constitution defines minimum requirements for correctness and maintainability
- Performance optimizations must not compromise numerical accuracy
- All code changes must pass: compilation, unit tests, verification tests
- Breaking changes to public API require version bump and migration guide
- PETSc version compatibility: support last 2 major releases minimum

**Version**: 1.0.0 | **Ratified**: 2025-10-23 | **Last Amended**: 2025-10-23
**Version**: 1.1.0 | **Ratified**: 2025-10-23 | **Last Amended**: 2025-10-23
