#include <iostream>
#include <vector>

#include <Kokkos_Core.hpp>
#include <petsc.h>
#include <petscsys.h>
#include <petscversion.h>


#include <petscdm.h>
#include <petscdmda.h>

int main_1(int argc, char *argv[])
{
  DM              da, daX, daY;
  DMDALocalInfo   info;
  MPI_Comm        commX, commY;
  Vec             basisX, basisY;
  PetscScalar   **arrayX, **arrayY;
  const PetscInt *lx, *ly;
  PetscInt        M = 3, N = 3;
  PetscInt        p     = 1;
  PetscInt        numGP = 3;
  PetscInt        dof   = 2 * (p + 1) * numGP;
  PetscMPIInt     rank, subsize, subrank;
  
  static char help[] = "Tests DMDAVecGetArrayDOF()\n";
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, 0, help));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
  /* Create 2D DMDA */
  PetscCall(DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, M, N, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &da));
  PetscCall(DMSetFromOptions(da));
  PetscCall(DMSetUp(da));
  /* Create 1D DMDAs along two directions. */
  PetscCall(DMDAGetOwnershipRanges(da, &lx, &ly, NULL));
  PetscCall(DMDAGetLocalInfo(da, &info));
  /* Partitioning in the X direction makes a subcomm extending in the Y direction and vice-versa. */
  PetscCall(DMDAGetProcessorSubsets(da, DM_X, &commY));
  PetscCall(DMDAGetProcessorSubsets(da, DM_Y, &commX));
  PetscCallMPI(MPI_Comm_size(commX, &subsize));
  PetscCallMPI(MPI_Comm_rank(commX, &subrank));
  PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[%d]X subrank: %d subsize: %d\n", rank, subrank, subsize));
  PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));
  PetscCallMPI(MPI_Comm_size(commY, &subsize));
  PetscCallMPI(MPI_Comm_rank(commY, &subrank));
  PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[%d]Y subrank: %d subsize: %d\n", rank, subrank, subsize));
  PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));
  PetscCall(DMDACreate1d(commX, DM_BOUNDARY_NONE, info.mx, dof, 1, lx, &daX));
  PetscCall(DMSetUp(daX));
  PetscCall(DMDACreate1d(commY, DM_BOUNDARY_NONE, info.my, dof, 1, ly, &daY));
  PetscCall(DMSetUp(daY));
  /* Create 1D vectors for basis functions */
  PetscCall(DMGetGlobalVector(daX, &basisX));
  PetscCall(DMGetGlobalVector(daY, &basisY));
  /* Extract basis functions */
  PetscCall(DMDAVecGetArrayDOF(daX, basisX, &arrayX));
  PetscCall(DMDAVecGetArrayDOF(daY, basisY, &arrayY));
  /*arrayX[i][ndof]; */
  /*arrayY[j][ndof]; */
  PetscCall(DMDAVecRestoreArrayDOF(daX, basisX, &arrayX));
  PetscCall(DMDAVecRestoreArrayDOF(daY, basisY, &arrayY));
  /* Return basis vectors */
  PetscCall(DMRestoreGlobalVector(daX, &basisX));
  PetscCall(DMRestoreGlobalVector(daY, &basisY));
  /* Cleanup */
  PetscCallMPI(MPI_Comm_free(&commX));
  PetscCallMPI(MPI_Comm_free(&commY));
  PetscCall(DMDestroy(&daX));
  PetscCall(DMDestroy(&daY));
  PetscCall(DMDestroy(&da));
  PetscCall(PetscFinalize());
  return 0;
}

/*TEST

   test:
      nsize: 2

TEST*/


int main(int argc, char **argv) {

    main_1(argc, argv);
  PetscErrorCode ierr = PetscInitialize(&argc, &argv, NULL, NULL);
  if (ierr) {
    std::cerr << "Failed to initialize PETSc" << std::endl;
    return ierr;
  }

  char version[128];
  ierr = PetscGetVersion(version, sizeof(version));
  if (!ierr) {
    PetscPrintf(PETSC_COMM_WORLD, "PETSc successfully initialized.\n");
  }

  // Finalize PETSc
  ierr = PetscFinalize();
  if (ierr) {
    std::cerr << "Failed to finalize PETSc" << std::endl;
    return ierr;
  }

  // Initialize Kokkos
  Kokkos::initialize(argc, argv);
  {
    Kokkos::print_configuration(std::cout);

    const int N = 1000;
    using view_t = Kokkos::View<double *>;
    view_t a("a", N), b("b", N), c("c", N);

    // Fill a and b
    Kokkos::parallel_for(
        "init", N, KOKKOS_LAMBDA(const int i) {
          a(i) = 1.0;
          b(i) = 2.0;
        });

    // c = a + b
    Kokkos::parallel_for(
        "add", N, KOKKOS_LAMBDA(const int i) { c(i) = a(i) + b(i); });

    // Reduce sum of c on host
    double sum = 0.0;
    Kokkos::parallel_reduce(
        "sum", N, KOKKOS_LAMBDA(const int i, double &lsum) { lsum += c(i); },
        sum);

    std::cout << "Kokkos computed sum = " << sum << " (expected " << 3.0 * N
              << ")\n";
  }
  Kokkos::finalize();


  return 0;
}
