#include "antenna_pattern.h"

#include <arts_conversions.h>
#include <lagrange_interp.h>

#include <algorithm>
#include <boost/math/distributions/normal.hpp>
#include <format>
#include <limits>
#include <stdexcept>
#include <variant>

namespace sensor {
namespace {
constexpr Numeric response_eps = 64 * std::numeric_limits<Numeric>::epsilon();
using identity = lagrange_interp::identity;

Numeric scale(Numeric x) { return std::max<Numeric>(1.0, std::abs(x)); }

bool is_close(Numeric a, Numeric b) {
  return std::abs(a - b) <= response_eps * std::max(scale(a), scale(b));
}

bool is_within_support(const AscendingGrid& grid, Numeric x) {
  if (grid.empty()) return false;
  return x >= grid.front() - response_eps * scale(grid.front()) and
         x <= grid.back() + response_eps * scale(grid.back());
}

template <lagrange_interp::transformer Transform>
using LookupLag =
    std::variant<lagrange_interp::lag_t<0, Transform>,
                 lagrange_interp::lag_t<1, Transform>>;

template <lagrange_interp::transformer Transform, typename Grid>
LookupLag<Transform> lookup_lag(const Grid& grid, Numeric x) {
  if (grid.size() == 1) return grid.template lag<0, Transform>(x);
  return grid.template lag<1, Transform>(x);
}

template <typename ZenithGrid, typename AzimuthGrid>
Numeric lookup_interpolate_impl(const ZenithGrid& zenith_grid,
                                const AzimuthGrid& azimuth_grid,
                                ConstMatrixView lookup_response,
                                Numeric delta_zenith,
                                Numeric delta_azimuth) {
  return std::visit(
      [&lookup_response](const auto& zenith_lag, const auto& azimuth_lag) {
        return lagrange_interp::interp(
            lookup_response, zenith_lag, azimuth_lag);
      },
      lookup_lag<identity>(zenith_grid, delta_zenith),
  lookup_lag<identity>(azimuth_grid, delta_azimuth));
}

Numeric lookup_interpolate(const AscendingGrid& zenith_grid,
                           const AscendingGrid& azimuth_grid,
                           ConstMatrixView lookup_response,
                           Numeric delta_zenith,
                           Numeric delta_azimuth) {
  if (not is_within_support(zenith_grid, delta_zenith) or
      not is_within_support(azimuth_grid, delta_azimuth)) {
    return 0.0;
  }

  return std::visit(
      [&lookup_response](const auto& zenith_lag, const auto& azimuth_lag) {
        return lagrange_interp::interp(
            lookup_response, zenith_lag, azimuth_lag);
      },
      lookup_lag<identity>(zenith_grid, delta_zenith),
      lookup_lag<identity>(azimuth_grid, delta_azimuth));
}

Numeric gaussian_component(Numeric delta, Numeric sigma) {
  using gauss = boost::math::normal_distribution<Numeric>;
  using boost::math::pdf;

  static const gauss unit_normal(0.0, 1.0);
  static const Numeric unit_peak = pdf(unit_normal, 0.0);

  return pdf(unit_normal, delta / sigma) / unit_peak;
}

void assert_positive_width(Numeric value, std::string_view name) {
  if (value <= 0) {
    throw std::invalid_argument(
        std::format("{} must be positive.  Got: {}.", name, value));
  }
}
}  // namespace

AntennaPattern AntennaPattern::pencil() {
  AntennaPattern pattern;
  pattern.set_pencil_beam();
  return pattern;
}

AntennaPattern AntennaPattern::gaussian(Numeric zenith_std,
                                        Numeric azimuth_std) {
  AntennaPattern pattern;
  pattern.set_gaussian(zenith_std, azimuth_std);
  return pattern;
}

AntennaPattern AntennaPattern::gaussian_fwhm(Numeric zenith_fwhm,
                                             Numeric azimuth_fwhm) {
  AntennaPattern pattern;
  pattern.set_gaussian_fwhm(zenith_fwhm, azimuth_fwhm);
  return pattern;
}

AntennaPattern AntennaPattern::lookup(const AscendingGrid& zenith_grid,
                                      const AscendingGrid& azimuth_grid,
                                      const Matrix& response_lookup) {
  AntennaPattern pattern;
  pattern.set_lookup(zenith_grid, azimuth_grid, response_lookup);
  return pattern;
}

AntennaPattern AntennaPattern::lookup(const ZenGrid& zenith_grid,
                                      const AziGrid& azimuth_grid,
                                      const Matrix& response_lookup) {
  AntennaPattern pattern;
  pattern.set_lookup(zenith_grid, azimuth_grid, response_lookup);
  return pattern;
}

void AntennaPattern::set_pencil_beam() {
  type                = AntennaType::PencilBeam;
  sigma_zenith        = 0.0;
  sigma_azimuth       = 0.0;
  lookup_uses_angular_grids = false;
  lookup_zenith_grid  = AscendingGrid{};
  lookup_azimuth_grid = AscendingGrid{};
  lookup_response.resize(0, 0);
}

void AntennaPattern::set_gaussian(Numeric zenith_std, Numeric azimuth_std) {
  assert_positive_width(zenith_std, "Zenith standard deviation");
  assert_positive_width(azimuth_std, "Azimuth standard deviation");

  type                = AntennaType::Gaussian;
  sigma_zenith        = zenith_std;
  sigma_azimuth       = azimuth_std;
  lookup_uses_angular_grids = false;
  lookup_zenith_grid  = AscendingGrid{};
  lookup_azimuth_grid = AscendingGrid{};
  lookup_response.resize(0, 0);
}

void AntennaPattern::set_gaussian_fwhm(Numeric zenith_fwhm,
                                       Numeric azimuth_fwhm) {
  set_gaussian(Conversion::fwhm2std(zenith_fwhm),
               Conversion::fwhm2std(azimuth_fwhm));
}

void AntennaPattern::set_lookup(const AscendingGrid& zenith_grid,
                                const AscendingGrid& azimuth_grid,
                                ConstMatrixView response_lookup_) {
  if (zenith_grid.empty() or azimuth_grid.empty()) {
    throw std::invalid_argument("Lookup antenna grids must not be empty.");
  }

  if (response_lookup_.nrows() != static_cast<Index>(zenith_grid.size()) or
      response_lookup_.ncols() != static_cast<Index>(azimuth_grid.size())) {
    throw std::invalid_argument(std::format(
        "Lookup antenna response must have shape ({}, {}).  Got ({}, {}).",
        zenith_grid.size(),
        azimuth_grid.size(),
        response_lookup_.nrows(),
        response_lookup_.ncols()));
  }

  type                = AntennaType::Lookup;
  sigma_zenith        = 0.0;
  sigma_azimuth       = 0.0;
  lookup_uses_angular_grids = false;
  lookup_zenith_grid  = zenith_grid;
  lookup_azimuth_grid = azimuth_grid;
  lookup_response     = Matrix(response_lookup_);
}

void AntennaPattern::set_lookup(const ZenGrid& zenith_grid,
                                const AziGrid& azimuth_grid,
                                ConstMatrixView response_lookup_) {
  set_lookup(AscendingGrid{zenith_grid.vec()},
             AscendingGrid{azimuth_grid.vec()},
             response_lookup_);
  lookup_uses_angular_grids = true;
}

Numeric AntennaPattern::operator()(Numeric delta_zenith,
                                   Numeric delta_azimuth) const {
  switch (type) {
    case AntennaType::PencilBeam:
      return is_close(delta_zenith, 0.0) and is_close(delta_azimuth, 0.0) ? 1.0
                                                                          : 0.0;

    case AntennaType::Gaussian:
      return gaussian_component(delta_zenith, sigma_zenith) *
             gaussian_component(delta_azimuth, sigma_azimuth);

    case AntennaType::Lookup:
      if (lookup_uses_angular_grids) {
        if (not ZenGrid::is_valid(delta_zenith) or
            not is_within_support(lookup_zenith_grid, delta_zenith) or
            not is_within_support(lookup_azimuth_grid, delta_azimuth)) {
          return 0.0;
        }

        return lookup_interpolate_impl(ZenGrid{lookup_zenith_grid.vec()},
                                       AziGrid{lookup_azimuth_grid.vec()},
                                       lookup_response,
                                       delta_zenith,
                                       delta_azimuth);
      }

      return lookup_interpolate(lookup_zenith_grid,
                                lookup_azimuth_grid,
                                lookup_response,
                                delta_zenith,
                                delta_azimuth);

    default:
      throw std::invalid_argument(
          std::format("Unsupported antenna type: {}.", type));
  }
}

AntennaPatternGriddedField AntennaPattern::response(
    const AscendingGrid& zenith_grid, const AscendingGrid& azimuth_grid) const {
  Matrix data(zenith_grid.size(), azimuth_grid.size());

  for (Size i = 0; i < zenith_grid.size(); i++) {
    for (Size j = 0; j < azimuth_grid.size(); j++) {
      data[i, j] = (*this)(zenith_grid[i], azimuth_grid[j]);
    }
  }

  return {.data_name  = "antenna-pattern"s,
          .data       = std::move(data),
          .grid_names = std::array<String, 2>{"dzen"s, "dazi"s},
          .grids = std::array<AscendingGrid, 2>{zenith_grid, azimuth_grid}};
}

AntennaPatternGriddedField AntennaPattern::normalized_response(
    const AscendingGrid& zenith_grid, const AscendingGrid& azimuth_grid) const {
  auto out = response(zenith_grid, azimuth_grid);

  Numeric total = 0.0;
  for (Index i = 0; i < out.data.nrows(); i++) {
    for (Index j = 0; j < out.data.ncols(); j++) {
      total += out.data[i, j];
    }
  }

  if (total > 0) out.data /= total;

  out.data_name = "antenna-pattern-normalized"s;
  return out;
}

AntennaPatternGriddedField AntennaPattern::raw_sensor(
    const AscendingGrid& dzen_grid, const AscendingGrid& dazi_grid) const {
  auto out      = response(dzen_grid, dazi_grid);
  out.data_name = "antenna-raw-sensor"s;
  return out;
}

AntennaPatternGriddedField AntennaPattern::normalized_raw_sensor(
    const AscendingGrid& dzen_grid, const AscendingGrid& dazi_grid) const {
  auto out      = normalized_response(dzen_grid, dazi_grid);
  out.data_name = "antenna-raw-sensor-normalized"s;
  return out;
}

static_assert(AntennaPatternSelection<AntennaPattern>);
}  // namespace sensor