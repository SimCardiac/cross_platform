# Feature Specification: Poisson Equation Unit Tests

**Feature Branch**: `002-poisson-unit-tests`  
**Created**: 2025-10-23  
**Status**: Draft  
**Input**: User description: "Unit test for the Poisson equation based on src/poisson_DMStag_2D.cpp"

## User Scenarios & Testing (mandatory)

### User Story 1 - Run Poisson tests locally (Priority: P1)

As a developer, I can run the Poisson unit tests with a single command and get a clear summary, so I can quickly verify changes.

**Why this priority**: Ensures rapid, reliable feedback when modifying core elliptic routines.

**Independent Test**: Running the test entrypoint executes Poisson tests only and prints pass/fail with counts; no extra setup beyond project standard.

**Acceptance Scenarios**:

1. Given a clean checkout, When the developer runs the test entrypoint, Then all Poisson tests execute and a summary is printed with pass/fail counts.
2. Given an intentional change that breaks a stencil identity, When tests run, Then at least one test fails with an actionable message.

---

### User Story 2 - Verify operator consistency and convergence (Priority: P1)

As a maintainer, I can validate discrete operator consistency and observe expected convergence on grid refinement, so I trust the numerical foundation.

**Why this priority**: Operator correctness and convergence are critical to solution accuracy.

**Independent Test**: Manufactured-solution tests compute errors and compare observed order with expected rate; passes independently of other tests.

**Acceptance Scenarios**:

1. Given a manufactured solution on multiple grid resolutions, When errors are computed, Then the observed rate matches expected order within 10%.
2. Given a constant RHS and appropriate BCs, When the solution is computed, Then residual norms drop below the target tolerance.

---

### User Story 3 - Boundary condition validation (Priority: P2)

As a reviewer, I can confirm the implementation respects Dirichlet and Neumann boundary conditions on representative faces and corners, so edge behaviors are correct.

**Why this priority**: Boundary handling errors commonly cause subtle inaccuracies.

**Independent Test**: Boundary-only test fixtures verify values/fluxes against expectations within tolerance.

**Acceptance Scenarios**:

1. Given Dirichlet BCs on all sides, When the solver runs, Then boundary values match prescribed values within tolerance.
2. Given Neumann BCs on selected sides, When the solver runs, Then normal fluxes match prescribed values within tolerance and solution is unique up to a constant when applicable.

### Edge Cases

- Homogeneous RHS with Dirichlet BCs (zero solution)
- Constant RHS and mixed BCs (sanity check for offsets/fluxes)
- Small domains (2x2, 3x3) to catch indexing/stagger issues
- Pressure/nullspace-like invariances (constant offsets where applicable)
- Highly anisotropic grid spacing (if configurable) impacting error rates

## Requirements (mandatory)

### Scope Boundaries

- In scope: Unit tests for Poisson equation covering discrete operator identities, boundary enforcement, solver residuals, and grid-refinement convergence using manufactured solutions.
- Out of scope: Full transient simulations, performance benchmarking, visualization validation.

### Assumptions

- Existing Poisson implementation is available in the repository.
- Deterministic tests; acceptable numerical variance handled via explicit tolerances.
- Minimal domains and canonical boundary sets are sufficient to expose defects.

### Functional Requirements

- FR-001: Provide a single entrypoint to run only the Poisson unit tests and report a concise summary.
- FR-002: Include tests verifying discrete Laplacian/operator identities on minimal grids.
- FR-003: Include grid-refinement tests using a manufactured solution measuring observed order against expected.
- FR-004: Validate Dirichlet and Neumann boundary conditions on representative faces/corners within tolerance.
- FR-005: Require solver residual norms to drop below a specified tolerance on representative RHS/BC setups.
- FR-006: Ensure tests complete quickly (target ≤ 60 seconds on a typical developer laptop for the Poisson subset).
- FR-007: Ensure tests are deterministic across supported environments without external services.
- FR-008: Provide brief documentation describing how to run Poisson tests and interpret failures.

### Acceptance Criteria

- AC-001 (FR-001): Running the documented command executes only Poisson tests and prints a summary with exit status reflecting pass/fail.
- AC-002 (FR-002): Operator identity tests pass with relative error ≤ 1e-8 (or absolute ≤ 1e-12 for near-zero norms).
- AC-003 (FR-003): Observed convergence rate on grid refinement matches expected order within 10% on manufactured solutions.
- AC-004 (FR-004): Boundary-condition tests pass with errors within 1e-8 relative to prescribed values/fluxes.
- AC-005 (FR-005): Solver residual L2 norm ≤ 1e-10 on representative RHS/BCs for minimal grids.
- AC-006 (FR-006): Poisson subset test runtime ≤ 60 seconds on a contemporary laptop during review.
- AC-007 (FR-007): Re-running the Poisson tests three times yields identical pass/fail outcomes within tolerance windows.
- AC-008 (FR-008): A short README or doc section explains running Poisson tests and diagnosing common failures.

### Key Entities (if data involved)

- Test Case: Minimal scenario with domain, RHS, BC specification, and expected quantitative outcomes.
- Reference Solution: Manufactured expressions/identities used to compute expected results and error bounds.
- Configuration: Parameters (domain size, boundary types) defining a test’s conditions and tolerances.

## Success Criteria (mandatory)

### Measurable Outcomes

- SC-001: Developers can run Poisson tests locally and see a clear summary in ≤ 60 seconds.
- SC-002: Operator consistency and convergence properties are validated by tests that fail when tolerances/orders are not met.
- SC-003: Zero intermittent failures in three consecutive runs on supported platforms using the documented process.
- SC-004: Reviewers rely on the suite to detect regressions; seeded defects are caught during review.
