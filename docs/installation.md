# Installation Guide

This guide covers setting up the development environment on **Apple Silicon (M2) macOS** and building the project from source.

---

## Prerequisites

- macOS running on Apple Silicon (M1/M2/M3/M4)
- [Homebrew](https://brew.sh) package manager
- Internet connection for downloading dependencies

---

## 1. Install Host Compiler (GCC via Homebrew)

Apple's built-in Clang lacks Fortran support needed by PETSc, and does not support OpenMP. A host GCC compiler is required as a bootstrap toolchain.

```bash
brew install gcc@14
```

This installs:
| Binary | Path |
|--------|------|
| C compiler | `/opt/homebrew/bin/gcc-14` |
| C++ compiler | `/opt/homebrew/bin/g++-14` |
| Fortran compiler | `/opt/homebrew/bin/gfortran-14` |

> **Why not build GCC directly with Spack?** Spack needs a working C/C++/Fortran host compiler to bootstrap. Homebrew provides a pre-built GCC that serves this role reliably.

---

## 2. Install Spack

Download and set up the Spack package manager:

```bash
wget https://github.com/spack/spack/releases/download/v1.0.2/spack-1.0.2.tar.gz
tar -xzf spack-1.0.2.tar.gz && rm spack-1.0.2.tar.gz
mv spack-1.0.2 $HOME/spack
```

Add Spack to your shell (append to `~/.zshrc` for persistence):

```bash
. $HOME/spack/share/spack/setup-env.sh
```

---

## 3. Register Compilers with Spack

Tell Spack about the Homebrew GCC compilers. Create or edit `~/.spack/packages.yaml`:

```yaml
packages:
  gcc:
    externals:
    - spec: gcc@14.3.0 languages:='c,c++,fortran'
      prefix: /opt/homebrew
      extra_attributes:
        compilers:
          c: /opt/homebrew/bin/gcc-14
          cxx: /opt/homebrew/bin/g++-14
          fortran: /opt/homebrew/bin/gfortran-14
```

Then register them:

```bash
spack compiler find
```

Verify with:

```bash
spack compiler list
```

---

## 4. Create Spack Environment

Create a Spack environment for this project and add the required dependencies:

```bash
spack env create cross
spack env activate cross
spack add openmpi@5.0.8 %gcc@14
spack add openblas@0.3.27 %gcc@14
spack add kokkos@4.6.01 +openmp %gcc@14
spack add petsc@3.23.4 +kokkos +openmp %gcc@14
```

Alternatively, use the provided `spack.yaml` in the project root:

```bash
spack env create cross spack.yaml
spack env activate cross
```

The environment file (`spack.yaml`) specifies:

```yaml
spack:
  specs:
  - openmpi@5.0.8
  - openblas@0.3.27
  - kokkos@4.6.01 +openmp
  - petsc@3.23.4 +kokkos +openmp
  - py-pip
  view: true
  concretizer:
    unify: true
```

> **Note on OpenBLAS**: Version `0.3.27` is a known-stable version. Newer versions may fail to build on macOS.

Install all packages:

```bash
spack install
```

> ⚠️ This step compiles PETSc, Kokkos, OpenMPI, and OpenBLAS from source. It may take 30–60 minutes depending on your machine.

---

## 5. Activate the Environment

Every time you open a new terminal, source Spack and activate the environment:

```bash
. $HOME/spack/share/spack/setup-env.sh
spack env activate cross
```

The `spack view` creates a merged directory of all installed packages so CMake can discover them automatically.

---

## 6. Generate Analytical Solution Headers

Run the Python code generator to produce `src/analytical/unsteady.h`:

```bash
cd /path/to/cross_platform
python tools/unsteady.py
```

---

## 7. Build with CMake

```bash
mkdir -p build && cd build

# Clean any stale cache, then configure with GCC-14
rm -f CMakeCache.txt
CC=/opt/homebrew/bin/gcc-14 \
CXX=/opt/homebrew/bin/g++-14 \
FC=/opt/homebrew/bin/gfortran-14 \
cmake ..

# Build (using all CPU cores)
make -j$(sysctl -n hw.ncpu)
```

> **Why specify compilers manually?** CMake defaults to the system compiler (`/usr/bin/gcc`, which is Apple Clang). Clang lacks OpenMP, causing `find_package(Kokkos)` to fail. Explicitly setting `CC`/`CXX`/`FC` ensures the correct GCC toolchain is used.

---

## 8. Run Tests

Run convergence verification tests:

```bash
# Poisson equation
./poisson_DMStag_2D -poisson_check_error

# Heat equation
./heat_DMStag_2D -heat_check_error

# Stokes equation
./stokes_DMStag_2D -stokes_check_error

# Multi-process MPI
mpirun -np 2 ./poisson_DMStag_2D -ksp_monitor_short -ksp_converged_reason
```

Run CTest test suites:

```bash
ctest -R heat_dmstag_convergence -V
ctest -R heat_dmstag_edge_convergence -V
ctest -R heat_dmstag_CN_temporal_convergence -V
ctest -R heat_dmstag_temporal_convergence -V
```

---

## Troubleshooting

### Spack 
```
SystemError: buffer overflow
```
**Fix**: Do not use the latest version of Python.
```bash
export SPACK_PYTHON=/opt/homebrew/opt/python@3.12/bin/python3.12
```
### CMake cannot find OpenMP

```
Could NOT find OpenMP_CXX (missing: OpenMP_CXX_FLAGS OpenMP_CXX_LIB_NAMES)
```

**Fix**: Ensure `CC` and `CXX` point to Homebrew GCC when running `cmake`:

```bash
CC=/opt/homebrew/bin/gcc-14 CXX=/opt/homebrew/bin/g++-14 cmake ..
```

### OpenBLAS fails to build

Use the known-stable version:

```bash
spack add openblas@0.3.27 %gcc@14
```

### PETSc requires Kokkos with OpenMP

Both Kokkos and PETSc must be built with `+openmp` enabled:

```bash
spack add kokkos@4.6.01 +openmp %gcc@14
spack add petsc@3.23.4 +kokkos +openmp %gcc@14
```

### GPU Acceleration

Kokkos Metal backend is not yet available. GPU acceleration on Apple Silicon will be supported once upstream Kokkos adds Metal support or via MoltenVK translation layer.