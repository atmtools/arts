#ifndef UTILS_H_
#define UTILS_H_
#include <lagrange_interp.h>
#include <matpack.h>

#include <algorithm>
#include <iostream>
#include <vector>

namespace scattering {
/** The linear Lagrange weights of x on grid
 *
 * Beyond the ends of grid they extrapolate linearly, without limit; a grid of
 * one point gives its one value everywhere.
 */
lagrange_interp::lag_t<1> linear_lag(const Vector& grid, Numeric x);

//! linear_lag of every point of new_grid
std::vector<lagrange_interp::lag_t<1>> linear_lags(const Vector& grid, const Vector& new_grid);

//! The linear Lagrange weights of x on grid, clamped to its ends: beyond them, the value at the nearest end
lagrange_interp::lag_t<1> clamped_lag(const Vector& grid, Numeric x);

//! Lags that keep each of the n indices of an axis that is not interpolated
std::vector<lagrange_interp::lag_t<0>> identity_lags(Index n);

//! Lags that take the given index of an axis for each new point: the axis is selected, not interpolated
std::vector<lagrange_interp::lag_t<0>> selection_lags(const std::vector<Index>& indices);

Index digitize(const Vector& boundaries, Numeric value);

}  // namespace scattering

#endif  // UTILS_H_
