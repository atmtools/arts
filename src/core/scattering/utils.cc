#include "utils.h"

namespace scattering {

namespace {
std::span<const Numeric> as_span(const Vector& x) { return {x.data_handle(), x.size()}; }

//! All the weight on the one point of a one-point grid
lagrange_interp::lag_t<1> one_point_lag() {
  lagrange_interp::lag_t<1> lag;
  lag.indx = {0, 0};
  lag.data = {1.0, 0.0};
  return lag;
}
}  // namespace

lagrange_interp::lag_t<1> linear_lag(const Vector& grid, Numeric x) {
  if (grid.empty()) arts::user_error("Cannot interpolate on an empty grid");
  if (grid.size() == 1) return one_point_lag();
  if (grid[0] < grid[1]) return lagrange_interp::lag_t<1>(as_span(grid), x, lagrange_interp::ascending_grid_t{});
  return lagrange_interp::lag_t<1>(as_span(grid), x, lagrange_interp::descending_grid_t{});
}

std::vector<lagrange_interp::lag_t<1>> linear_lags(const Vector& grid, const Vector& new_grid) {
  if (grid.empty()) arts::user_error("Cannot interpolate on an empty grid");
  if (grid.size() == 1) return std::vector<lagrange_interp::lag_t<1>>(new_grid.size(), one_point_lag());
  // An extrapolation limit of 0 checks nothing: the weights extrapolate linearly beyond the ends, without limit
  return lagrange_interp::make_lags<1, lagrange_interp::grid_identity>(
      as_span(grid), as_span(new_grid), 0.0, "scattering data grid");
}

lagrange_interp::lag_t<1> clamped_lag(const Vector& grid, Numeric x) {
  if (grid.empty()) arts::user_error("Cannot interpolate on an empty grid");
  const auto [low, upp] = std::minmax(grid.front(), grid.back());
  return linear_lag(grid, std::clamp(x, low, upp));
}

std::vector<lagrange_interp::lag_t<0>> identity_lags(Index n) {
  std::vector<lagrange_interp::lag_t<0>> out(n);
  for (Index i = 0; i < n; ++i) {
    out[i].indx = {i};
    out[i].data = {1.0};
  }
  return out;
}

std::vector<lagrange_interp::lag_t<0>> selection_lags(const std::vector<Index>& indices) {
  std::vector<lagrange_interp::lag_t<0>> out(indices.size());
  for (Size i = 0; i < indices.size(); ++i) {
    out[i].indx = {indices[i]};
    out[i].data = {1.0};
  }
  return out;
}

/**! Find bin index for a given value.
     *
     * Returns the 0-based bin index of a given value.
     *
     * @param boundaries: A vector containing the n + 1 bin-boundaries of the n bins in increasing order.
     * @param value: The value to sort into the bins.
     *
     * @return  The index of the bin or:
     *     - -1 if value is smaller than any of the bin boundaries
     *     - n + 1 if value is larger than any of the bin boundaries.
     */
Index digitize(const Vector& boundaries, Numeric value) {
  // Handle boundary cases
  if (value <= boundaries.front()) { return -1; }
  if (value >= boundaries.back()) { return boundaries.size(); }

  // Find the interval using std::lower_bound
  auto  it    = std::lower_bound(boundaries.begin(), boundaries.end(), value);
  Index index = std::distance(boundaries.begin(), it);
  return index - 1;
}

}  // namespace scattering
