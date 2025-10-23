# Implementation Plan: Poisson Equation Unit Tests

**Branch**: `002-poisson-unit-tests` | **Date**: 2025-10-23 | **Spec**: /Users/pengfei/Documents/GitHub/cross_platform/specs/002-poisson-unit-tests/spec.md
**Input**: Feature specification from `/specs/002-poisson-unit-tests/spec.md`

**Note**: This template is filled in by the `/speckit.plan` command. See `.specify/templates/commands/plan.md` for the execution workflow.

## Summary

Deliver a focused Poisson equation unit test suite based on the existing implementation (`src/poisson_DMStag_2D.cpp`).
Tests validate discrete operator identities, boundary condition enforcement, solver residuals, and grid-refinement
convergence via manufactured solutions. Support multiple RHS/exact-solution pairs and accept all parameters via
command line flags. Build runs by sourcing `~/.bashrc` first, then using CMake + Make (out-of-source `build/`).

## Technical Context

<!--
  ACTION REQUIRED: Replace the content in this section with the technical details
  for the project. The structure here is presented in advisory capacity to guide
  the iteration process.
-->

**Language/Version**: C++17 (C++20 if toolchain allows)  
**Primary Dependencies**: PETSc (Vec/Mat/DM/KSP), CMake, Make  
**Storage**: N/A  
**Testing**: NEEDS CLARIFICATION (CTest + lightweight CLI tests vs. GoogleTest)  
**Target Platform**: macOS, Linux  
**Project Type**: Single project (CLI + tests)  
**Performance Goals**: Poisson test subset completes ≤ 60 seconds locally  
**Constraints**: Deterministic tests; no network; CLI-only parameterization; source `~/.bashrc` before build  
**Scale/Scope**: Unit-test level only (operators, BCs, residuals, refinement convergence)

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

- PETSc API Compliance: PASS (Poisson tests leverage existing PETSc-based implementation)
- Modern C++ Standards: PASS (C++17 minimum; RAII encouraged)
- Template Programming Strategy: PASS/NA (Poisson tests do not preclude templates; extensible design kept)
- Numerical Correctness: PASS (MMS-based convergence, residual thresholds, BC validation)
- Phased Delivery & Coverage: PASS (Phase 1 is Poisson; aligns with constitution roadmap)
- Discretization/Time: PASS (Elliptic; discrete Laplacian identities validated)
- Solver Configuration: PASS (Residual tolerance and iteration limits configurable)
- Boundary Conditions: PASS (Dirichlet/Neumann validated)
- QA & Docs: PASS (Unit + verification tests; quickstart/docs to be added)

## Project Structure

### Documentation (this feature)

```text
specs/[###-feature]/
├── plan.md              # This file (/speckit.plan command output)
├── research.md          # Phase 0 output (/speckit.plan command)
├── data-model.md        # Phase 1 output (/speckit.plan command)
├── quickstart.md        # Phase 1 output (/speckit.plan command)
├── contracts/           # Phase 1 output (/speckit.plan command)
└── tasks.md             # Phase 2 output (/speckit.tasks command - NOT created by /speckit.plan)
```

### Source Code (repository root)
<!--
  ACTION REQUIRED: Replace the placeholder tree below with the concrete layout
  for this feature. Delete unused options and expand the chosen structure with
  real paths (e.g., apps/admin, packages/something). The delivered plan must
  not include Option labels.
-->

```text
src/
├── poisson_DMStag_2D.cpp           # existing implementation (reference)
└── ...

tests/
└── unit/
  └── poisson/
    ├── test_poisson_cli.sh     # optional CLI smoke test
    └── (CTest or gtest sources TBD)

specs/002-poisson-unit-tests/
├── spec.md
├── plan.md
├── research.md
├── data-model.md
├── quickstart.md
└── contracts/
  └── cli.md
```

**Structure Decision**: Single-project structure; tests kept under `tests/unit/poisson/`.
Contracts expressed as CLI parameter spec (`specs/.../contracts/cli.md`).

## Complexity Tracking

> **Fill ONLY if Constitution Check has violations that must be justified**

| Violation | Why Needed | Simpler Alternative Rejected Because |
|-----------|------------|-------------------------------------|
| [e.g., 4th project] | [current need] | [why 3 projects insufficient] |
| [e.g., Repository pattern] | [specific problem] | [why direct DB access insufficient] |

---

## Phase 0: Outline & Research

Unknowns to resolve:
- Testing framework choice: CTest-only vs. GoogleTest (NEEDS CLARIFICATION)
- Exact CLI param schema for multiple RHS/exact pairs
- Default tolerance values and expected orders

Research tasks:
- Research CTest integration for C++ unit tests and CLI wrapping
- Define 2–3 manufactured solutions with analytic u and corresponding f for Poisson
- Establish CLI flags for: grid size, case selection, tolerances, BC types

Output: research.md consolidating decisions, rationale, and alternatives.

## Phase 1: Design & Contracts

Artifacts to generate:
- data-model.md: Entities (TestCase, Boundary, Tolerances, Pair)
- contracts/cli.md: CLI schema for parameters and multiple pairs handling
- quickstart.md: Build/run instructions including sourcing `~/.bashrc`, CMake + Make, sample commands
- Agent context update (executed)

Re-evaluate Constitution Check after Phase 1; all gates expected to remain PASS.
