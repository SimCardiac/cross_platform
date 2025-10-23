# Project Structure: cross_platform

## Repository Layout

```
cross_platform/
├── .github/
│   ├── copilot-instructions.md         # AI assistant context
│   └── prompts/                        # Workflow prompts
│       ├── speckit.specify.prompt.md
│       ├── speckit.constitution.prompt.md
│       ├── speckit.plan.prompt.md
│       ├── speckit.tasks.prompt.md
│       └── speckit.implement.prompt.md
│
├── .specify/
│   ├── memory/
│   │   └── constitution.md             # Project governance & standards
│   ├── scripts/bash/                   # Workflow automation scripts
│   └── templates/                      # Document templates
│
├── specs/                              # Feature specifications
│   ├── 001-add-unit-tests/            # Broad test suite spec
│   └── 002-poisson-unit-tests/        # Current feature (Poisson tests)
│       ├── spec.md                    # Requirements & user stories
│       ├── plan.md                    # Implementation plan
│       ├── tasks.md                   # Task breakdown
│       ├── research.md                # Technical decisions
│       ├── data-model.md              # Data entities
│       ├── quickstart.md              # Usage guide
│       ├── contracts/
│       │   └── cli.md                 # CLI interface spec
│       └── checklists/
│           └── requirements.md        # Quality checklist
│
├── src/                                # Source code (see src/README.md)
│   ├── poisson/                       # Poisson module (NEW)
│   │   ├── poisson_mms.{hpp,cpp}
│   │   ├── poisson_metrics.{hpp,cpp}
│   │   ├── poisson_solver.{hpp,cpp}
│   │   └── tests/
│   │       └── poisson_tests_cli.cpp
│   ├── io/                            # I/O utilities
│   ├── dmstag/                        # DMStag helpers
│   ├── examples/                      # Example programs
│   ├── wyh/                           # Specialized solvers
│   ├── space_filling_curve/
│   ├── poisson_DMStag_2D.cpp         # Original solvers
│   ├── poisson_DM_2D.cpp
│   ├── heat_DMStag_2D.cpp
│   ├── steady_stokes_DMStag_2D.cpp
│   └── ...
│
├── tests/                              # Test suite
│   └── unit/
│       └── poisson/                   # Poisson unit tests
│           ├── CMakeLists.txt         # CTest configuration
│           ├── test_poisson_cli.sh    # Shell smoke test
│           └── .gitkeep
│
├── build/                              # CMake build directory (gitignored)
│   ├── poisson_tests_cli              # Test executable
│   └── ...
│
├── docs/                               # Documentation
│   ├── heat_convergence_results.md
│   └── heat_equation.md
│
├── CMakeLists.txt                      # Main build configuration
├── .gitignore                          # Git ignore patterns
├── spack.yaml                          # Spack environment
├── readme.md
└── ...
```

## Key Improvements in This Feature

### 1. Organized Test Infrastructure

**Before:**
- All files flat in `src/`
- No clear separation of concerns

**After:**
```
src/poisson/                 # Module boundary
├── poisson_mms.*           # Manufactured solutions
├── poisson_metrics.*       # Error computation
├── poisson_solver.*        # Solver wrapper
└── tests/                  # Test executables
    └── poisson_tests_cli.cpp
```

### 2. CTest Integration

**New files:**
- `tests/unit/poisson/CMakeLists.txt` - Registers 4 CTest suites
- `tests/unit/poisson/test_poisson_cli.sh` - Shell-based smoke tests

**Usage:**
```bash
cd build
ctest -R poisson           # Run all Poisson tests
ctest -R poisson -V        # Verbose output
```

### 3. Specification-Driven Development

**Workflow:**
1. `specs/002-poisson-unit-tests/spec.md` - What to build
2. `specs/002-poisson-unit-tests/plan.md` - How to build it
3. `specs/002-poisson-unit-tests/tasks.md` - Step-by-step tasks
4. Implementation in `src/poisson/`
5. Tests in `tests/unit/poisson/`

## Build Targets

| Target | Type | Location | Dependencies |
|--------|------|----------|--------------|
| `poisson_mms` | Library | `src/poisson/` | - |
| `poisson_metrics` | Library | `src/poisson/` | PETSc |
| `poisson_solver` | Library | `src/poisson/` | PETSc |
| `poisson_tests_cli` | Executable | `src/poisson/tests/` | All above libs |
| `ex_poisson_stagger` | Executable | `src/` | PETSc |
| `ex_poisson_center` | Executable | `src/` | PETSc |

## Next Steps for Organization

As the test suite expands to Heat and Stokes:

```
src/
├── common/              # Shared utilities (future)
│   ├── metrics.*       # Generic error metrics
│   └── mms.*           # Generic MMS framework
├── poisson/            # Existing
├── heat/               # Future: Heat equation
│   ├── heat_solver.*
│   └── tests/
└── stokes/             # Future: Stokes flow
    ├── stokes_solver.*
    └── tests/

tests/unit/
├── poisson/            # Existing
├── heat/               # Future
└── stokes/             # Future
```

## Documentation References

- Project governance: `.specify/memory/constitution.md`
- Poisson spec: `specs/002-poisson-unit-tests/spec.md`
- Source layout: `src/README.md`
- Build guide: `specs/002-poisson-unit-tests/quickstart.md`
