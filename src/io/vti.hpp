#pragma once
#include <mpi.h>
#include <string>
#include <vector>

namespace io {

struct VTIOptions {
  // If true, data is written under <PointData>, otherwise <CellData>
  bool pointData = true;
  // Name of the single scalar array
  std::string arrayName = "data";
  // Origin and spacing
  double origin[3]  = {0.0, 0.0, 0.0};
  double spacing[3] = {1.0, 1.0, 1.0};
    // Number of ghost layers to add at piece boundaries (for point-data tiling)
  // Many VTK readers request 1-column overlap between pieces; set to 1 to avoid UpdateExtent gaps.
  int ghostLevel = 1;
};

// Write a 2D parallel VTI consisting of one piece per MPI rank and a master PVTI on rank 0.
// All extents are for point-data indexing: [iMin, iMax, jMin, jMax, kMin, kMax]
// global_nx, global_ny are the global point counts in i/j. Local tile starts at (local_x0, local_y0)
// with size (local_nx, local_ny). local_data has local_nx*local_ny entries in i-fastest order.
// Files are written as: <stem>_<rank>.vti and master <stem>.pvti in the current working directory
void write_vti_parallel_2d(const std::string &stem,
                           MPI_Comm comm,
                           int global_nx, int global_ny,
                           int local_x0, int local_y0,
                           int local_nx, int local_ny,
                           const double *local_data,
                           const VTIOptions &opts = {});

} // namespace io
