---

description: "Task list for Poisson equation unit tests feature"
---

# Tasks: Poisson Equation Unit Tests

**Input**: Design documents from `/specs/002-poisson-unit-tests/`
**Prerequisites**: plan.md (required), spec.md (required for user stories), research.md, data-model.md, contracts/

**Tests**: This feature explicitly requests unit tests; test tasks are included per user story.

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

## Path Conventions

- Single project: `src/`, `tests/` at repository root
- Tests organized under `tests/unit/poisson/`

---

## Phase 1: Setup (Shared Infrastructure)

Purpose: Prepare testing infrastructure and project scaffolding needed by all stories.

- [x] T001 Enable CTest in top-level CMake at `CMakeLists.txt` (add `enable_testing()` and basic config)
- [x] T002 [P] Create tests directory structure `tests/unit/poisson/` with `.gitkeep`
- [x] T003 [P] Add `tests/unit/poisson/CMakeLists.txt` to register CLI-based tests via `add_test(...)`
- [x] T004 Add a top-level include for `tests/unit/poisson/CMakeLists.txt` from root `CMakeLists.txt`

---

## Phase 2: Foundational (Blocking Prerequisites)

Purpose: Core building blocks that MUST be complete before any user story work.

- [x] T005 Create manufactured solutions module `src/poisson_mms.hpp` and `src/poisson_mms.cpp` (pairs: poly2, sinpi, cospi)
- [x] T006 [P] Create error norms utilities `src/poisson_metrics.hpp` and `src/poisson_metrics.cpp` (L2/Linf, residual)
- [x] T007 Implement CLI test runner `src/poisson_tests_cli.cpp` with flags per `specs/002-poisson-unit-tests/contracts/cli.md`
- [x] T008 Wire PETSc/KSP setup in runner `src/poisson_tests_cli.cpp` (tolerances, max iters, residual capture)
- [x] T009 Refactor or wrap `src/poisson_DMStag_2D.cpp` usage: expose a callable solve path `run_poisson_case(params)` (new function)
- [x] T010 [P] Add nullspace handling for pure Neumann in runner (PETSc `MatNullSpace`) in `src/poisson_tests_cli.cpp`
- [ ] T011 [P] Add optional JSON summary output `--output` in `src/poisson_tests_cli.cpp`
- [x] T012 Register a smoke test in CTest `tests/unit/poisson/CMakeLists.txt` invoking runner with `--grid 8 8 --pair poly2`

Checkpoint: Foundation ready — user stories can proceed.

---

## Phase 3: User Story 1 - Run Poisson tests locally (Priority: P1) 🎯 MVP

Goal: A single command runs only Poisson tests and prints a concise summary with correct exit code.

Independent Test: From `build/`, `ctest -R poisson -j` executes Poisson tests only and shows pass/fail counts; exit code reflects failures.

### Tests for User Story 1 (Requested)

- [x] T013 [P] [US1] Add CTest label `poisson` and test definitions in `tests/unit/poisson/CMakeLists.txt`
- [x] T014 [P] [US1] Create CLI smoke test script `tests/unit/poisson/test_poisson_cli.sh` (invokes runner with one pair)

### Implementation for User Story 1

- [x] T015 [US1] Ensure runner prints per-case summary and final aggregated status in `src/poisson_tests_cli.cpp`
- [x] T016 [US1] Document run command in `specs/002-poisson-unit-tests/quickstart.md` (ctest and direct runner)
- [x] T017 [US1] Add CTest entry to run the shell smoke test `tests/unit/poisson/test_poisson_cli.sh`

Checkpoint: US1 independently runnable and verifiable.

---

## Phase 4: User Story 2 - Verify operator consistency and convergence (Priority: P1)

Goal: Validate discrete operator identities and expected convergence on refinement using manufactured solutions.

Independent Test: Running convergence tests over multiple grids produces observed order within 10% of expected; identity tests meet error tolerances.

### Tests for User Story 2 (Requested)

- [ ] T018 [P] [US2] Add identity test CTest entry (e.g., constant RHS) in `tests/unit/poisson/CMakeLists.txt` invoking runner
- [ ] T019 [P] [US2] Add refinement test CTest entry (e.g., sinpi pair, grids 16/32/64) in `tests/unit/poisson/CMakeLists.txt`

### Implementation for User Story 2

- [ ] T020 [P] [US2] Implement identity check path in `src/poisson_tests_cli.cpp` (reports relative error ≤ 1e-8)
- [ ] T021 [US2] Implement refinement loop `--refine N` in `src/poisson_tests_cli.cpp` and compute slope of log(error) vs log(h)`
- [ ] T022 [US2] Integrate MMS pairs usage from `src/poisson_mms.cpp` into solve/evaluation pipeline
- [ ] T023 [US2] Ensure residual threshold check (≤ 1e-10) and fail-fast behavior in `src/poisson_tests_cli.cpp`

Checkpoint: US2 independently validates consistency and convergence.

---

## Phase 5: User Story 3 - Boundary condition validation (Priority: P2)

Goal: Confirm Dirichlet/Neumann boundary conditions are enforced correctly on faces and corners.

Independent Test: Boundary-only fixtures compare values/fluxes against expectations within tolerance; pure Neumann uniqueness handled.

### Tests for User Story 3 (Requested)

- [ ] T024 [P] [US3] Add Dirichlet-all-sides test in `tests/unit/poisson/CMakeLists.txt` invoking runner with `--bc Dirichlet Dirichlet Dirichlet Dirichlet`
- [ ] T025 [P] [US3] Add mixed Neumann/Dirichlet test in `tests/unit/poisson/CMakeLists.txt` invoking runner with `--bc Dirichlet Neumann Dirichlet Neumann`

### Implementation for User Story 3

- [ ] T026 [P] [US3] Add BC specification parsing and application in `src/poisson_tests_cli.cpp`
- [ ] T027 [US3] Implement boundary value/flux checks vs MMS in `src/poisson_tests_cli.cpp` (≤ 1e-8 relative)
- [ ] T028 [US3] Ensure nullspace-aware validation for Neumann cases (compare up to constant) in `src/poisson_tests_cli.cpp`

Checkpoint: US3 independently validates BC enforcement.

---

## Phase N: Polish & Cross-Cutting Concerns

Purpose: Finishings, docs, and robustness improvements across stories.

- [x] T029 [P] Add developer docs for CLI flags in `specs/002-poisson-unit-tests/contracts/cli.md` (ensure synced with implementation)
- [ ] T030 Code cleanup and comments in `src/poisson/tests/poisson_tests_cli.cpp`, `src/poisson/poisson_mms.*`, `src/poisson/poisson_metrics.*`
- [ ] T031 [P] Add parallel CTest runs config and timeouts in `tests/unit/poisson/CMakeLists.txt`
- [ ] T032 Validate `specs/002-poisson-unit-tests/quickstart.md` by a fresh build and run

**Code Organization (Completed):**
- [x] Reorganized `src/` structure: Created `src/poisson/` module
- [x] Moved test infrastructure to `src/poisson/tests/`
- [x] Updated CMakeLists.txt paths
- [x] Created `src/README.md` and `PROJECT_STRUCTURE.md` documentation

---

## Dependencies & Execution Order

### Phase Dependencies

- Setup (Phase 1): No dependencies — can start immediately
- Foundational (Phase 2): Depends on Setup completion — BLOCKS all user stories
- User Stories (Phase 3+): Depend on Foundational completion; US1 and US2 can proceed in parallel after Phase 2; US3 can also run in parallel, but benefits from US2 utilities
- Polish (Final Phase): Depends on selected user stories being complete

### User Story Dependencies

- User Story 1 (P1): Starts after Phase 2; independent
- User Story 2 (P1): Starts after Phase 2; independent of US1 (shares foundation)
- User Story 3 (P2): Starts after Phase 2; independent of US1; may reuse US2 utilities but not strictly required

### Within Each User Story

- Tests (included) should be written and verified against expected failures before full implementation
- Implement CLI and computation paths, then wire into CTest
- Ensure each story’s criteria can be validated on its own

### Parallel Opportunities

- T002, T003 in Phase 1 can run in parallel
- In Phase 2: T006, T010, T011 can run in parallel; T012 can proceed once T007 is ready
- User stories: US1, US2, US3 can be developed by different contributors in parallel once Phase 2 completes

---

## Parallel Example: User Story 2

```
# Parallelize identity and refinement test wiring (after runner exists):
Task: "Add identity test CTest entry ..." (T018)
Task: "Add refinement test CTest entry ..." (T019)

# Meanwhile implement computation paths:
Task: "Implement identity check path ..." (T020)
Task: "Implement refinement loop ..." (T021)
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup
2. Complete Phase 2: Foundational (BLOCKS all stories)
3. Complete Phase 3: User Story 1
4. STOP and VALIDATE: `ctest -R poisson` returns summary and correct exit code

### Incremental Delivery

1. Setup + Foundational → base ready
2. Add US1 → test independently → demo
3. Add US2 → test independently → demo
4. Add US3 → test independently → demo

---

# Report

- Generated file: `/Users/pengfei/Documents/GitHub/cross_platform/specs/002-poisson-unit-tests/tasks.md`
- Total task count: 32
- Task count per user story:
  - US1: 3 implementation + 2 tests = 5 (T013–T017)
  - US2: 4 implementation + 2 tests = 6 (T018–T023)
  - US3: 3 implementation + 2 tests = 5 (T024–T028)
  - Setup + Foundational + Polish: 16 (T001–T012, T029–T032)
- Parallel opportunities identified: Phase 1 (T002, T003), Phase 2 (T006, T010, T011), US2 tests (T018–T019), US1 tests (T013–T014), US3 tests (T024–T025)
- Independent test criteria:
  - US1: `ctest -R poisson` prints summary; exit status reflects failures
  - US2: Identity error thresholds; convergence rate within 10% of expected
  - US3: Boundary values/fluxes within tolerance; Neumann uniqueness handled
- Suggested MVP scope: User Story 1 (US1) only
- Format validation: All tasks follow required checklist format (checkbox + TaskID + optional [P] + optional [US#] + file path)
