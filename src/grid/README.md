# Space-Filling Curve Grid Ordering

A C++ library for mapping structured grid indices using space-filling curves, designed to be compatible with PETSc DMStag/DMDA.

## Features

- **Three ordering strategies:**
  - **XYZ**: Natural row-major ordering (i fastest, then j, then k)
  - **Morton (Z-order)**: Bit-interleaving for better cache locality
  - **Hilbert**: Optimal spatial locality with continuous curve

- **Unified interface:** All orderings provide the same API for bidirectional mapping between (i,j,k) coordinates and linear indices

- **DMStag/DMDA compatible:** Adapter class for seamless integration with PETSc's distributed structured grids

- **Template-based design:** Easy to extend with custom ordering strategies

## API Overview

### Basic Usage

```cpp
#include "grid/space_filling_order.hpp"
using namespace grid;

// Create ordering for 16x16 grid
auto order = create_ordering(OrderingType::Hilbert, 16, 16);

// Convert spatial coordinates to linear index
int64_t idx = order->to_linear(3, 5, 0);  // i=3, j=5, k=0

// Convert linear index back to coordinates
int i, j, k;
order->to_ijk(idx, i, j, k);
```

### With PETSc DMStag/DMDA

```cpp
#include "grid/space_filling_order.hpp"
#include <petscdmda.h>

// Create DMDA
DM da;
DMDACreate2d(comm, ..., &da);

// Create space-filling curve ordering
auto sfc = create_ordering(OrderingType::Hilbert, mx, my);

// Map local grid points to SFC indices
DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
for (int j = ys; j < ys + ym; j++) {
  for (int i = xs; i < xs + xm; i++) {
    int64_t sfc_idx = sfc->to_linear(i, j, 0);
    // Use sfc_idx for data access, I/O ordering, etc.
  }
}
```

### With Multiple DOF (DMStag-style)

```cpp
// Create adapter for grid with 3 DOF per point
auto order = create_ordering(OrderingType::Morton, 8, 8);
DMStagSFCAdapter adapter(std::move(order), 3);

// Map (i,j,k,component) to linear index
int64_t idx = adapter.to_linear(2, 3, 0, 1);  // point (2,3), component 1

// Reverse mapping
int i, j, k, c;
adapter.to_ijkc(idx, i, j, k, c);
```

## Requirements

- Grid dimensions must be powers of 2 for Morton and Hilbert orderings
- Hilbert ordering requires square (2D) or cubic (3D) grids
- 3D Hilbert is not yet implemented (use Morton or XYZ for 3D)

## Performance Characteristics

| Ordering | Spatial Locality | Computation Cost | Use Case |
|----------|------------------|------------------|----------|
| XYZ      | Poor (jumps at row boundaries) | O(1) | Simple, no constraints |
| Morton   | Good (Z-pattern) | O(log n) | General purpose, 2D/3D |
| Hilbert  | Excellent (continuous) | O(log n) | Best locality, 2D grids |

**Benchmark (16x16 grid, average Manhattan distance between consecutive indices):**
- XYZ: 1.004
- Morton: 1.498
- Hilbert: 1.250

Lower values indicate better cache locality.

## Build

```bash
cd build
cmake ..
make grid_sfc ex_sfc ex_sfc_dmda
```

## Examples

### Basic test
```bash
./ex_sfc
```
Shows ordering visualization, DMStag adapter usage, and locality comparison.

### PETSc integration
```bash
mpirun -np 4 ./ex_sfc_dmda
```
Demonstrates using space-filling curves with PETSc DMDA for domain decomposition.

## Use Cases

1. **I/O optimization**: Write grid data in SFC order for better disk locality
2. **Load balancing**: Partition along SFC for better subdomain connectivity
3. **Cache efficiency**: Access data in SFC order to improve cache hit rate
4. **Visualization**: Output data in SFC order for progressive rendering

## Integration with DMStag

DMStag uses staggered grid with multiple DOF per location. This library provides:

```cpp
// Example: 2D staggered grid with cell centers + edges
// DMStag manages ghost cells and MPI communication
// Space-filling curves provide ordering for owned cells

DMStagGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL, NULL, NULL, NULL);

auto sfc = create_ordering(OrderingType::Hilbert, mx, my);
DMStagSFCAdapter adapter(std::move(sfc), dof_per_cell);

// Access pattern optimized for SFC
for (int64_t sfc_idx = 0; sfc_idx < adapter.size() / dof_per_cell; ++sfc_idx) {
  int i, j, k;
  adapter.ordering()->to_ijk(sfc_idx, i, j, k);
  
  if (i >= xs && i < xs + xm && j >= ys && j < ys + ym) {
    // This is a local point, process all DOF
    for (int c = 0; c < dof_per_cell; ++c) {
      int64_t linear_idx = adapter.to_linear(i, j, k, c);
      // Access data[linear_idx]
    }
  }
}
```

## Extending

To add a new ordering strategy:

1. Derive from `SpaceFillingOrder`
2. Implement `to_linear()` and `to_ijk()`
3. Add to `create_ordering()` factory
4. Update `OrderingType` enum

## References

- Morton (Z-order): [Wikipedia](https://en.wikipedia.org/wiki/Z-order_curve)
- Hilbert curve: [Wikipedia](https://en.wikipedia.org/wiki/Hilbert_curve)
- PETSc DMStag: [PETSc Manual](https://petsc.org/release/manual/dmstag/)
