# Feature Specification: Unit Test Suite for Poisson, Heat, and Incompressible Flow

**Feature Branch**: `001-add-unit-tests`  
**Created**: 2025-10-23  
**Status**: Draft  
**Input**: User description: "Add unit tests for incompressible Navier-Stokes solver"

## User Scenarios & Testing *(mandatory)*

<!--
  IMPORTANT: User stories should be PRIORITIZED as user journeys ordered by importance.
  Each user story/journey must be INDEPENDENTLY TESTABLE - meaning if you implement just ONE of them,
  you should still have a viable MVP (Minimum Viable Product) that delivers value.
  
  Assign priorities (P1, P2, P3, etc.) to each story, where P1 is the most critical.
  Think of each story as a standalone slice of functionality that can be:
  - Developed independently
  - Tested independently
  - Deployed independently
  - Demonstrated to users independently
-->

### User Story 1 - Run all tests locally (Priority: P1)

As a developer, I can run the full unit test suite with a single command and immediately see clear, actionable failure messages, so I can safely refactor core solver code.

**Why this priority**: Enables safe iteration and reduces regressions in core numerical routines.

**Independent Test**: Running the test command completes successfully and returns a pass/fail summary without any additional setup beyond standard project instructions.

**Acceptance Scenarios**:

1. **Given** a clean checkout, **When** the developer runs the test entrypoint, **Then** all unit tests execute and a clear summary is printed with pass/fail counts.
2. **Given** an intentionally introduced error in a core routine, **When** tests run, **Then** at least one test fails with a message pointing to the violated property and a minimal reproducer.


### User Story 2 - Validate numerical properties (Priority: P1)

As a maintainer, I can verify key numerical properties (e.g., divergence-free projection, boundary condition application, discrete operator consistency) are enforced via unit tests, so changes don’t silently degrade solution quality.

**Why this priority**: Protects correctness of physical constraints central to incompressible flow.

**Independent Test**: Individual tests assert norms/identities against thresholds and pass independently of other tests.

**Acceptance Scenarios**:

1. **Given** a manufactured velocity-pressure field, **When** the projection step is applied, **Then** the divergence norm falls below a defined threshold.
2. **Given** canonical boundary setups (no-slip, free-slip, periodic), **When** discrete operators are evaluated near boundaries, **Then** the results satisfy their boundary definitions within tolerance.


### User Story 3 - Catch regressions automatically (Priority: P2)

As a reviewer, I rely on automated tests to block merging changes that break correctness or assumptions, so review time focuses on design rather than manual verification.

**Why this priority**: Ensures continuous quality and reduces manual verification burden.

**Independent Test**: A change that alters a discrete operator or solver contract triggers one or more failing tests until fixed.

**Acceptance Scenarios**:

1. **Given** a change to a discrete operator, **When** tests run, **Then** any deviation from expected identities (symmetry, consistency) causes a failure.
2. **Given** changes to time-stepping or coupling, **When** a stability-related test executes, **Then** it fails if a defined stability condition is violated on a standard setup.



### User Story 4 - Step-by-step phased delivery (Priority: P1)

As a project owner, I want the test suite to be delivered in phases — first Poisson, then Heat, then Incompressible — so quality is built incrementally and each phase has clear exit criteria before moving on.

**Why this priority**: De-risks delivery by validating foundations (elliptic and parabolic) before the full coupled system.

**Independent Test**: Each phase can be run and validated independently; the next phase cannot pass CI until the previous phase’s exit criteria are met.

**Acceptance Scenarios**:

1. Given Phase 1 (Poisson) tests pass, When Phase 2 (Heat) tests are introduced, Then Phase 1 remains green and Phase 2 must pass in isolation before Phase 3 starts.
2. Given Phase 2 (Heat) tests pass, When Phase 3 (Incompressible) tests are introduced, Then a projection-method suite with first- and second-order temporal accuracy tests must pass before feature completion.

### Edge Cases

- Zero or uniform velocity field inputs (should remain divergence-free and unchanged after projection)
- Very high or very low viscosity (diffusion-dominated vs. advection-dominated regimes)
- Very small and relatively large time steps at stability boundaries for explicit components
- Periodic boundaries and mixed boundary conditions on different faces
- Pressure nullspace handling (constant offset invariance in pressure-related validations)
- Thermal boundary combinations for Heat tests (Dirichlet/Neumann mixes on orthogonal faces)
## Requirements *(mandatory)*

### Scope Boundaries

- In scope: Unit tests for Poisson, Heat, and Incompressible (projection method) covering discrete operators, boundary condition application, projection step, time-stepping primitives, and basic I/O invariants (e.g., shape/metadata)
- Out of scope: Full end-to-end long transient simulations, performance benchmarking, and visualization validation beyond smoke checks
### Assumptions

- Dimensionless test problems unless a scale is specified; default tolerance references unit scales
- Deterministic test inputs; acceptable numerical variance handled via explicit tolerances
- Test fixtures provide minimal domains and canonical boundary condition sets
### Functional Requirements

- FR-001: Provide a single entrypoint to run all unit tests and report a concise summary.
- FR-002: Include tests for discrete gradient, divergence, and Laplacian operators verifying consistency and expected identities on minimal grids.
- FR-003: Include tests that validate boundary conditions (no-slip, free-slip, inflow/outflow, periodic) on representative faces and corners.
- FR-004: Provide a projection-step test that reduces velocity-field divergence below a specified threshold on standard setups.
- FR-005: Include a pressure Poisson validation that achieves a residual norm below a specified threshold for representative right-hand sides.
- FR-006: Include time-stepping tests (e.g., diffusion-only and advection-only toy problems) that verify stability/accuracy bounds on standard setups.
- FR-007: Provide manufactured/analytical solutions or discrete identities enabling numerical error measurement with clear pass/fail thresholds.
- FR-008: Ensure tests complete quickly (target under 3 minutes on a typical developer laptop for the full suite).
- FR-009: Ensure tests are deterministic across supported platforms within stated tolerances and do not require internet or external services.
- FR-010: Provide minimal documentation in the repository explaining how to run tests and interpret failures.
- FR-011: Provide Poisson equation tests validating discrete operator consistency and solver convergence with residual norms below thresholds on refined grids.
- FR-012: Provide Heat equation tests validating time-integration accuracy with first-order (e.g., backward Euler) and second-order (e.g., Crank–Nicolson) methods using manufactured solutions and reporting observed order.
- FR-013: Provide Incompressible Navier–Stokes tests based on a projection method, including both first- and second-order temporal accuracy variants, verifying divergence reduction and temporal convergence against manufactured solutions.
- FR-014: Enforce phased delivery: Phase 1 (Poisson) must pass all acceptance criteria before Phase 2 (Heat) merges; Phase 2 must pass before Phase 3 (Incompressible) merges.
- FR-015: Provide comprehensive documentation covering test scope, assumptions, fixture domains, boundary conditions, expected orders, tolerances, and how to extend the suite.

### Acceptance Criteria

- AC-001 (FR-001): Running the documented test command executes all tests and prints a summary with exit status reflecting pass/fail.
- AC-002 (FR-002): On minimal domains, discrete operator tests pass with relative error ≤ 1e-8 (or absolute ≤ 1e-12 when norms are small).
- AC-003 (FR-003): Boundary-condition tests pass with error within 1e-8 relative to expected face/corner values on canonical cases.
- AC-004 (FR-004): Post-projection velocity divergence L2 norm ≤ 1e-8 for unit-scale manufactured inputs.
- AC-005 (FR-005): Pressure equation residual L2 norm ≤ 1e-10 for representative right-hand sides on minimal grids.
- AC-006 (FR-006): Stability/accuracy tests pass with errors within defined bounds for standard time steps (documented in the test case).
- AC-007 (FR-007): Tests using manufactured/analytical references document the formula and report measured order when applicable.
- AC-008 (FR-008): End-to-end test runtime ≤ 3 minutes on a contemporary laptop as measured during review.
- AC-009 (FR-009): Re-running the suite three times consecutively yields identical pass/fail outcomes and errors within tolerance windows.
- AC-010 (FR-010): A short README section explains how to run tests and how to diagnose failures.
- AC-011 (FR-011): On grid refinement, Poisson solver residuals converge below 1e-10 and observed rate matches discrete operator expectations within 10%.
- AC-012 (FR-012): Heat equation temporal error shows ≈1st order for BE and ≈2nd order for CN on manufactured solutions over at least two time-step halvings.
- AC-013 (FR-013): Projection tests achieve velocity divergence L2 ≤ 1e-8 post-projection; temporal convergence matches 1st/2nd order for respective variants.
- AC-014 (FR-014): CI blocks merging of Phase 2 until all Phase 1 criteria pass; similarly blocks Phase 3 until all Phase 2 criteria pass (documented gating logic).
- AC-015 (FR-015): Documentation includes sectioned coverage for Poisson, Heat, and Incompressible tests with explicit tolerances, norms, and expected orders; reviewers can locate and follow steps without external context.

### Key Entities *(include if feature involves data)*

- Test Case: Minimal, self-contained scenario with domain size, initial/boundary data, and expected quantitative outcomes.
- Reference Solution: Analytical/manufactured expressions or discrete identities used to compute expected results and error bounds.
- Configuration: Parameters (domain size, time step, boundary types) defining a test’s conditions and tolerances.


## Success Criteria *(mandatory)*

### Measurable Outcomes

