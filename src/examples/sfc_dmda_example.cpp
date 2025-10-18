#include "../grid/space_filling_order.hpp"
#include <petscdm.h>
#include <petscdmda.h>
#include <petscvec.h>
#include <iostream>
#include <vector>

using namespace grid;

/**
 * @brief Example: Use space-filling curve ordering with PETSc DMStag/DMDA
 * 
 * This demonstrates how to:
 * 1. Create a DMStag/DMDA for structured grid
 * 2. Use space-filling curve ordering to reorder data
 * 3. Maintain compatibility with PETSc's distributed arrays
 */

PetscErrorCode demonstrate_sfc_with_dmda(MPI_Comm comm, OrderingType order_type) {
  PetscFunctionBeginUser;
  DM             da;
  Vec            global, local;
  PetscScalar  **array;
  PetscInt       i, j, xs, ys, xm, ym, mx = 16, my = 16;
  PetscMPIInt    rank;
  
  PetscCallMPI(MPI_Comm_rank(comm, &rank));
  
  // Create 2D DMDA
  PetscCall(DMDACreate2d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                         DMDA_STENCIL_BOX, mx, my,
                         PETSC_DECIDE, PETSC_DECIDE,
                         1, 1, NULL, NULL, &da));
  PetscCall(DMSetFromOptions(da));
  PetscCall(DMSetUp(da));
  
  // Get local portion
  PetscCall(DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL));
  
  PetscCall(PetscPrintf(comm, "\n=== Testing %s ordering with DMDA ===\n",
                        order_type == OrderingType::XYZ ? "XYZ" :
                        order_type == OrderingType::Morton ? "Morton" : "Hilbert"));
  PetscCall(PetscSynchronizedPrintf(comm, 
      "[%d] Local subdomain: (%d,%d) size (%d,%d)\n",
      rank, xs, ys, xm, ym));
  PetscCall(PetscSynchronizedFlush(comm, PETSC_STDOUT));
  
  // Create global vector
  PetscCall(DMCreateGlobalVector(da, &global));
  PetscCall(DMCreateLocalVector(da, &local));
  
  // Initialize data using space-filling curve ordering
  PetscCall(DMDAVecGetArray(da, local, &array));
  
  try {
    // Create space-filling order for the global grid
    auto sfc_order = create_ordering(order_type, mx, my);
    
    // Fill local array with SFC index as value
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        int64_t sfc_idx = sfc_order->to_linear(i, j, 0);
        array[j][i] = (PetscScalar)sfc_idx;
      }
    }
    
    // Compute statistics for this ordering
    std::vector<int64_t> local_indices;
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        local_indices.push_back(sfc_order->to_linear(i, j, 0));
      }
    }
    
    // Check if indices are contiguous (good for memory access)
    bool contiguous = true;
    for (size_t k = 1; k < local_indices.size(); k++) {
      if (local_indices[k] != local_indices[k-1] + 1) {
        contiguous = false;
        break;
      }
    }
    
    PetscCall(PetscSynchronizedPrintf(comm,
        "[%d] SFC index range: [%lld, %lld], contiguous: %s\n",
        rank, 
        (long long)*std::min_element(local_indices.begin(), local_indices.end()),
        (long long)*std::max_element(local_indices.begin(), local_indices.end()),
        contiguous ? "yes" : "no"));
    PetscCall(PetscSynchronizedFlush(comm, PETSC_STDOUT));
    
  } catch (const std::exception& e) {
    PetscCall(PetscPrintf(comm, "Error: %s\n", e.what()));
  }
  
  PetscCall(DMDAVecRestoreArray(da, local, &array));
  
  // Local to global
  PetscCall(DMLocalToGlobal(da, local, INSERT_VALUES, global));
  
  // View vector (only for small grids)
  if (mx <= 8 && my <= 8) {
    PetscCall(PetscPrintf(comm, "\nGlobal vector (SFC indices):\n"));
    PetscCall(VecView(global, PETSC_VIEWER_STDOUT_WORLD));
  }
  
  // Cleanup
  PetscCall(VecDestroy(&local));
  PetscCall(VecDestroy(&global));
  PetscCall(DMDestroy(&da));
  
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc, char **argv) {
  PetscCall(PetscInitialize(&argc, &argv, NULL, 
      "Space-filling curve ordering with PETSc DMStag/DMDA\n"));
  
  MPI_Comm comm = PETSC_COMM_WORLD;
  PetscMPIInt size;
  PetscCallMPI(MPI_Comm_size(comm, &size));
  
  PetscCall(PetscPrintf(comm, 
      "Running with %d MPI ranks\n", size));
  
  // Test all three orderings
  PetscCall(demonstrate_sfc_with_dmda(comm, OrderingType::XYZ));
  PetscCall(demonstrate_sfc_with_dmda(comm, OrderingType::Morton));
  PetscCall(demonstrate_sfc_with_dmda(comm, OrderingType::Hilbert));
  
  PetscCall(PetscPrintf(comm, "\n"));
  PetscCall(PetscPrintf(comm, "=== Summary ===\n"));
  PetscCall(PetscPrintf(comm, "XYZ ordering: Natural row-major, may not be contiguous across ranks\n"));
  PetscCall(PetscPrintf(comm, "Morton ordering: Z-curve, better locality than XYZ\n"));
  PetscCall(PetscPrintf(comm, "Hilbert ordering: Best spatial locality, indices more contiguous\n"));
  PetscCall(PetscPrintf(comm, "\nFor optimal performance, partition domain along the space-filling curve.\n"));
  
  PetscCall(PetscFinalize());
  return 0;
}
