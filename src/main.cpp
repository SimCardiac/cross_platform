#include <iostream>
#include <vector>

#include <Kokkos_Core.hpp>
#include <petsc.h>
#include <petscsys.h>
#include <petscversion.h>

int main(int argc, char **argv) {
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
