#include "../grid/space_filling_order.hpp"
#include <iostream>
#include <iomanip>
#include <vector>

using namespace grid;

void print_ordering_2d(const SpaceFillingOrder* order, const std::string& name) {
  int nx, ny, nz;
  order->get_dims(nx, ny, nz);
  
  std::cout << "\n" << name << " Ordering (" << nx << "x" << ny << "):\n";
  std::cout << "Linear index -> (i, j) coordinates:\n";
  
  for (int64_t idx = 0; idx < std::min((int64_t)16, order->size()); ++idx) {
    int i, j, k;
    order->to_ijk(idx, i, j, k);
    std::cout << "  " << std::setw(3) << idx << " -> (" 
              << i << ", " << j << ")\n";
  }
  
  // Show grid visualization
  std::cout << "\nGrid visualization (showing linear index at each point):\n";
  std::vector<std::vector<int>> grid(ny, std::vector<int>(nx, -1));
  
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      grid[ny-1-j][i] = order->to_linear(i, j, 0);
    }
  }
  
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      std::cout << std::setw(4) << grid[j][i];
    }
    std::cout << "\n";
  }
}

void test_dmstag_adapter() {
  std::cout << "\n" << std::string(60, '=') << "\n";
  std::cout << "DMStag-compatible adapter test\n";
  std::cout << std::string(60, '=') << "\n";
  
  const int nx = 4, ny = 4;
  const int dof = 3; // 3 DOF per grid point (e.g., velocity u, v, pressure p)
  
  // Create adapter with Hilbert ordering
  auto order = create_ordering(OrderingType::Hilbert, nx, ny);
  DMStagSFCAdapter adapter(std::move(order), dof);
  
  std::cout << "Grid: " << nx << "x" << ny << " with " << dof << " DOF per point\n";
  std::cout << "Total size: " << adapter.size() << " elements\n\n";
  
  std::cout << "Sample mappings (i,j,k,c) -> linear index:\n";
  for (int j = 0; j < 2; ++j) {
    for (int i = 0; i < 2; ++i) {
      for (int c = 0; c < dof; ++c) {
        int64_t idx = adapter.to_linear(i, j, 0, c);
        std::cout << "  (" << i << "," << j << ",0," << c << ") -> " << idx << "\n";
      }
    }
  }
  
  std::cout << "\nReverse mapping (first 12 linear indices):\n";
  for (int64_t idx = 0; idx < 12; ++idx) {
    int i, j, k, c;
    adapter.to_ijkc(idx, i, j, k, c);
    std::cout << "  " << idx << " -> (" << i << "," << j << "," << k << "," << c << ")\n";
  }
}

void compare_orderings() {
  std::cout << "\n" << std::string(60, '=') << "\n";
  std::cout << "Comparing different space-filling curves\n";
  std::cout << std::string(60, '=') << "\n";
  
  const int n = 4;
  
  auto xyz = create_ordering(OrderingType::XYZ, n, n);
  auto morton = create_ordering(OrderingType::Morton, n, n);
  auto hilbert = create_ordering(OrderingType::Hilbert, n, n);
  
  print_ordering_2d(xyz.get(), "XYZ (Natural)");
  print_ordering_2d(morton.get(), "Morton (Z-order)");
  print_ordering_2d(hilbert.get(), "Hilbert");
}

void benchmark_locality() {
  std::cout << "\n" << std::string(60, '=') << "\n";
  std::cout << "Spatial locality comparison\n";
  std::cout << std::string(60, '=') << "\n";
  
  const int n = 16;
  
  auto xyz = create_ordering(OrderingType::XYZ, n, n);
  auto morton = create_ordering(OrderingType::Morton, n, n);
  auto hilbert = create_ordering(OrderingType::Hilbert, n, n);
  
  // Measure average Manhattan distance between consecutive indices
  auto calc_avg_distance = [n](const SpaceFillingOrder* order) -> double {
    double sum = 0.0;
    for (int64_t idx = 0; idx < order->size() - 1; ++idx) {
      int i1, j1, k1, i2, j2, k2;
      order->to_ijk(idx, i1, j1, k1);
      order->to_ijk(idx + 1, i2, j2, k2);
      sum += std::abs(i2 - i1) + std::abs(j2 - j1);
    }
    return sum / (order->size() - 1);
  };
  
  std::cout << "Average Manhattan distance between consecutive indices:\n";
  std::cout << "  XYZ:     " << std::fixed << std::setprecision(3) 
            << calc_avg_distance(xyz.get()) << "\n";
  std::cout << "  Morton:  " << calc_avg_distance(morton.get()) << "\n";
  std::cout << "  Hilbert: " << calc_avg_distance(hilbert.get()) << "\n";
  std::cout << "\n(Lower is better for cache locality)\n";
}

int main() {
  try {
    compare_orderings();
    test_dmstag_adapter();
    benchmark_locality();
    
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "All tests completed successfully!\n";
    std::cout << std::string(60, '=') << "\n";
    
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  
  return 0;
}
