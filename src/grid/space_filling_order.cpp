#include "space_filling_order.hpp"
#include <cmath>
#include <stdexcept>

namespace grid {

// ============================================================================
// Factory
// ============================================================================
std::unique_ptr<SpaceFillingOrder> create_ordering(
    OrderingType type, int nx, int ny, int nz) {
  switch (type) {
  case OrderingType::XYZ:
    return std::make_unique<XYZOrder>(nx, ny, nz);
  case OrderingType::Morton:
    return std::make_unique<MortonOrder>(nx, ny, nz);
  case OrderingType::Hilbert:
    return std::make_unique<HilbertOrder>(nx, ny, nz);
  default:
    throw std::runtime_error("Unknown ordering type");
  }
}

// ============================================================================
// Morton Order Implementation
// ============================================================================

MortonOrder::MortonOrder(int nx, int ny, int nz)
    : nx_(nx), ny_(ny), nz_(nz) {
  // For Morton code, dimensions should be powers of 2
  auto is_pow2 = [](int n) {
    return n > 0 && (n & (n - 1)) == 0;
  };
  if (!is_pow2(nx) || !is_pow2(ny) || (nz > 1 && !is_pow2(nz))) {
    throw std::invalid_argument(
        "Morton ordering requires power-of-2 dimensions");
  }
}

uint32_t MortonOrder::part1by1(uint32_t n) {
  n &= 0x0000ffff;
  n = (n | (n << 8)) & 0x00FF00FF;
  n = (n | (n << 4)) & 0x0F0F0F0F;
  n = (n | (n << 2)) & 0x33333333;
  n = (n | (n << 1)) & 0x55555555;
  return n;
}

uint32_t MortonOrder::part1by2(uint32_t n) {
  n &= 0x000003ff;
  n = (n | (n << 16)) & 0xFF0000FF;
  n = (n | (n << 8))  & 0x0F00F00F;
  n = (n | (n << 4))  & 0xC30C30C3;
  n = (n | (n << 2))  & 0x49249249;
  return n;
}

uint32_t MortonOrder::compact1by1(uint32_t n) {
  n &= 0x55555555;
  n = (n ^ (n >> 1))  & 0x33333333;
  n = (n ^ (n >> 2))  & 0x0F0F0F0F;
  n = (n ^ (n >> 4))  & 0x00FF00FF;
  n = (n ^ (n >> 8))  & 0x0000FFFF;
  return n;
}

uint32_t MortonOrder::compact1by2(uint32_t n) {
  n &= 0x49249249;
  n = (n ^ (n >> 2))  & 0xC30C30C3;
  n = (n ^ (n >> 4))  & 0x0F00F00F;
  n = (n ^ (n >> 8))  & 0xFF0000FF;
  n = (n ^ (n >> 16)) & 0x000003FF;
  return n;
}

int64_t MortonOrder::to_linear(int i, int j, int k) const {
  if (nz_ == 1) {
    // 2D Morton
    return (int64_t)(part1by1(i) | (part1by1(j) << 1));
  } else {
    // 3D Morton
    return (int64_t)(part1by2(i) | (part1by2(j) << 1) | (part1by2(k) << 2));
  }
}

void MortonOrder::to_ijk(int64_t idx, int &i, int &j, int &k) const {
  if (nz_ == 1) {
    // 2D Morton
    i = compact1by1(idx);
    j = compact1by1(idx >> 1);
    k = 0;
  } else {
    // 3D Morton
    i = compact1by2(idx);
    j = compact1by2(idx >> 1);
    k = compact1by2(idx >> 2);
  }
}

// ============================================================================
// Hilbert Order Implementation
// ============================================================================

HilbertOrder::HilbertOrder(int nx, int ny, int nz)
    : nx_(nx), ny_(ny), nz_(nz) {
  // Check power of 2 and square/cube
  auto log2_check = [](int n) -> int {
    if (n <= 0 || (n & (n - 1)) != 0) return -1;
    int log = 0;
    while ((1 << log) < n) ++log;
    return log;
  };

  int logx = log2_check(nx);
  int logy = log2_check(ny);
  int logz = log2_check(nz);

  if (logx < 0 || logy < 0 || (nz > 1 && logz < 0)) {
    throw std::invalid_argument(
        "Hilbert ordering requires power-of-2 dimensions");
  }

  if (nz == 1) {
    if (nx != ny) {
      throw std::invalid_argument(
          "2D Hilbert ordering requires square grid (nx == ny)");
    }
    order_ = logx;
  } else {
    if (nx != ny || ny != nz) {
      throw std::invalid_argument(
          "3D Hilbert ordering requires cubic grid (nx == ny == nz)");
    }
    order_ = logx;
  }
}

void HilbertOrder::rot(int n, int &x, int &y, int rx, int ry) {
  if (ry == 0) {
    if (rx == 1) {
      x = n - 1 - x;
      y = n - 1 - y;
    }
    std::swap(x, y);
  }
}

void HilbertOrder::xy2d_2d(int x, int y, int64_t &d) const {
  d = 0;
  int s = nx_ / 2;
  while (s > 0) {
    int rx = (x & s) > 0 ? 1 : 0;
    int ry = (y & s) > 0 ? 1 : 0;
    d += (int64_t)s * s * ((3 * rx) ^ ry);
    rot(s, x, y, rx, ry);
    s /= 2;
  }
}

void HilbertOrder::d2xy_2d(int64_t d, int &x, int &y) const {
  x = y = 0;
  int s = 1;
  while (s < nx_) {
    int rx = 1 & (d / 2);
    int ry = 1 & (d ^ rx);
    rot(s, x, y, rx, ry);
    x += s * rx;
    y += s * ry;
    d /= 4;
    s *= 2;
  }
}

int64_t HilbertOrder::to_linear(int i, int j, int k) const {
  if (nz_ == 1) {
    int64_t d;
    xy2d_2d(i, j, d);
    return d;
  } else {
    // 3D Hilbert is more complex; for now throw or use fallback
    throw std::runtime_error("3D Hilbert not yet implemented, use 2D or XYZ/Morton");
  }
}

void HilbertOrder::to_ijk(int64_t idx, int &i, int &j, int &k) const {
  if (nz_ == 1) {
    d2xy_2d(idx, i, j);
    k = 0;
  } else {
    throw std::runtime_error("3D Hilbert not yet implemented");
  }
}

} // namespace grid
