#include <mpi.h>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>

#include "../io/vti.hpp"

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Global grid: Gx x Gy points
  const int Gx = 64, Gy = 64;

  // 1D block decomposition along x for simplicity
  int px = size; // processes along x
  int py = 1;    // processes along y
  if (Gy % py != 0 || Gx % px != 0) {
    if (rank == 0) std::cerr << "Grid not divisible by process grid" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int rx = rank % px;
  int ry = rank / px;

  int local_nx = Gx / px;
  int local_ny = Gy / py;
  int x0 = rx * local_nx;
  int y0 = ry * local_ny;

  std::vector<double> data(local_nx * local_ny);
  for (int j = 0; j < local_ny; ++j) {
    for (int i = 0; i < local_nx; ++i) {
      double x = (x0 + i) / double(Gx - 1);
      double y = (y0 + j) / double(Gy - 1);
      data[j * local_nx + i] = std::sin(2 * M_PI * x) * std::cos(2 * M_PI * y);
    }
  }

  io::VTIOptions opts;
  opts.arrayName = "phi";
  opts.pointData = true;
  opts.spacing[0] = opts.spacing[1] = 1.0;
  opts.spacing[2] = 1.0;

  io::write_vti_parallel_2d("output", MPI_COMM_WORLD,
                            Gx, Gy,
                            x0, y0,
                            local_nx, local_ny,
                            data.data(), opts);

  if (rank == 0) {
    std::cout << "Wrote output.pvti and per-rank output_*.vti" << std::endl;
  }

  MPI_Finalize();
  return 0;
}
