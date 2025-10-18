#pragma once

#include <cstdint>
#include <memory>
#include <vector>
#include <stdexcept>

namespace grid {

/**
 * @brief Space-filling curve ordering strategies for structured grids
 * 
 * This provides a unified interface for mapping between:
 * - (i, j, k) spatial coordinates
 * - linear index along a space-filling curve
 * 
 * Compatible with PETSc DMStag local/global indexing patterns.
 */
enum class OrderingType {
  XYZ,      ///< Natural row-major ordering (i fastest, then j, then k)
  Morton,   ///< Morton/Z-order curve (bit-interleaving)
  Hilbert   ///< Hilbert space-filling curve
};

/**
 * @brief Abstract base class for space-filling curve ordering
 */
class SpaceFillingOrder {
public:
  virtual ~SpaceFillingOrder() = default;

  /**
   * @brief Convert (i,j,k) coordinates to linear index
   * @param i x-coordinate [0, nx)
   * @param j y-coordinate [0, ny)
   * @param k z-coordinate [0, nz)
   * @return linear index in space-filling order
   */
  virtual int64_t to_linear(int i, int j, int k = 0) const = 0;

  /**
   * @brief Convert linear index to (i,j,k) coordinates
   * @param idx linear index
   * @param[out] i x-coordinate
   * @param[out] j y-coordinate
   * @param[out] k z-coordinate
   */
  virtual void to_ijk(int64_t idx, int &i, int &j, int &k) const = 0;

  /**
   * @brief Get grid dimensions
   */
  virtual void get_dims(int &nx, int &ny, int &nz) const = 0;

  /**
   * @brief Total number of grid points
   */
  virtual int64_t size() const = 0;

  /**
   * @brief Get ordering type
   */
  virtual OrderingType type() const = 0;
};

/**
 * @brief Factory function to create space-filling order
 * @param type Ordering strategy
 * @param nx Grid size in x (must be power of 2 for Morton/Hilbert)
 * @param ny Grid size in y (must be power of 2 for Morton/Hilbert)
 * @param nz Grid size in z (must be power of 2 for Morton/Hilbert, use 1 for 2D)
 */
std::unique_ptr<SpaceFillingOrder> create_ordering(
    OrderingType type, int nx, int ny, int nz = 1);

// ============================================================================
// XYZ Natural Ordering (i fastest)
// ============================================================================
class XYZOrder : public SpaceFillingOrder {
public:
  XYZOrder(int nx, int ny, int nz = 1)
      : nx_(nx), ny_(ny), nz_(nz) {}

  int64_t to_linear(int i, int j, int k = 0) const override {
    return (int64_t)k * ny_ * nx_ + (int64_t)j * nx_ + i;
  }

  void to_ijk(int64_t idx, int &i, int &j, int &k) const override {
    k = idx / (ny_ * nx_);
    int64_t rem = idx % (ny_ * nx_);
    j = rem / nx_;
    i = rem % nx_;
  }

  void get_dims(int &nx, int &ny, int &nz) const override {
    nx = nx_; ny = ny_; nz = nz_;
  }

  int64_t size() const override {
    return (int64_t)nx_ * ny_ * nz_;
  }

  OrderingType type() const override { return OrderingType::XYZ; }

private:
  int nx_, ny_, nz_;
};

// ============================================================================
// Morton (Z-order) Ordering
// ============================================================================
class MortonOrder : public SpaceFillingOrder {
public:
  MortonOrder(int nx, int ny, int nz = 1);

  int64_t to_linear(int i, int j, int k = 0) const override;
  void to_ijk(int64_t idx, int &i, int &j, int &k) const override;

  void get_dims(int &nx, int &ny, int &nz) const override {
    nx = nx_; ny = ny_; nz = nz_;
  }

  int64_t size() const override {
    return (int64_t)nx_ * ny_ * nz_;
  }

  OrderingType type() const override { return OrderingType::Morton; }

private:
  int nx_, ny_, nz_;

  // Helper: spread bits by inserting zeros between them
  static uint32_t part1by1(uint32_t n);   // for 2D
  static uint32_t part1by2(uint32_t n);   // for 3D
  static uint32_t compact1by1(uint32_t n);
  static uint32_t compact1by2(uint32_t n);
};

// ============================================================================
// Hilbert Ordering
// ============================================================================
class HilbertOrder : public SpaceFillingOrder {
public:
  HilbertOrder(int nx, int ny, int nz = 1);

  int64_t to_linear(int i, int j, int k = 0) const override;
  void to_ijk(int64_t idx, int &i, int &j, int &k) const override;

  void get_dims(int &nx, int &ny, int &nz) const override {
    nx = nx_; ny = ny_; nz = nz_;
  }

  int64_t size() const override {
    return (int64_t)nx_ * ny_ * nz_;
  }

  OrderingType type() const override { return OrderingType::Hilbert; }

private:
  int nx_, ny_, nz_;
  int order_;  // log2(nx), assuming nx=ny=nz=2^order

  // 2D Hilbert helpers
  static void rot(int n, int &x, int &y, int rx, int ry);
  void xy2d_2d(int x, int y, int64_t &d) const;
  void d2xy_2d(int64_t d, int &x, int &y) const;
};

// ============================================================================
// DMStag-compatible adapter
// ============================================================================

/**
 * @brief Adapter to use space-filling ordering with PETSc DMStag local arrays
 * 
 * DMStag uses (i,j,k,c) indexing where c is the component/DOF.
 * This adapter provides mapping between DMStag's local array layout and
 * space-filling curve ordering for spatial indices.
 */
class DMStagSFCAdapter {
public:
  /**
   * @param order Space-filling order for spatial indices
   * @param dof Number of DOF per grid point
   */
  DMStagSFCAdapter(std::unique_ptr<SpaceFillingOrder> order, int dof = 1)
      : order_(std::move(order)), dof_(dof) {}

  /**
   * @brief Map (i,j,k,component) to linear index in SFC-ordered array
   */
  int64_t to_linear(int i, int j, int k, int c) const {
    int64_t spatial_idx = order_->to_linear(i, j, k);
    return spatial_idx * dof_ + c;
  }

  /**
   * @brief Map linear index to (i,j,k,component)
   */
  void to_ijkc(int64_t idx, int &i, int &j, int &k, int &c) const {
    c = idx % dof_;
    int64_t spatial_idx = idx / dof_;
    order_->to_ijk(spatial_idx, i, j, k);
  }

  /**
   * @brief Total size including all DOF
   */
  int64_t size() const {
    return order_->size() * dof_;
  }

  /**
   * @brief Get underlying ordering
   */
  const SpaceFillingOrder* ordering() const { return order_.get(); }

  int dof() const { return dof_; }

private:
  std::unique_ptr<SpaceFillingOrder> order_;
  int dof_;
};

} // namespace grid
