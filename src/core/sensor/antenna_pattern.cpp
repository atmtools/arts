#include "antenna_pattern.h"

#include <arts_conversions.h>
#include <geodetic.h>

#include <cmath>

namespace sensor {

namespace {
struct AntennaBasis {
  Vector3 v;
  Vector3 h;
  Vector3 k;
};

[[nodiscard]] AntennaBasis antenna_basis(Vector2 bore_los) {
  using Conversion::cosd, Conversion::sind;

  const Numeric cza = cosd(bore_los[0]);
  const Numeric sza = sind(bore_los[0]);
  const Numeric caa = cosd(bore_los[1]);
  const Numeric saa = sind(bore_los[1]);

  return {
      .v = {-cza * saa, -cza * caa, sza},
      .h = {caa, -saa, 0.0},
      .k = {sza * saa, sza * caa, cza},
  };
}

[[nodiscard]] Vector3 antenna_frame_los(Vector2 local_los) {
  using Conversion::cosd, Conversion::sind;

  const Numeric cza = cosd(local_los[0]);
  const Numeric sza = sind(local_los[0]);
  const Numeric caa = cosd(local_los[1]);
  const Numeric saa = sind(local_los[1]);

  return {-sza * caa, sza * saa, cza};
}

[[nodiscard]] AntennaPatternField make_gaussian_field(ZenGrid zen_grid,
                                                      AziGrid azi_grid,
                                                      Numeric zenith_std,
                                                      Numeric azimuth_std,
                                                      Stokvec weight) {
  ARTS_USER_ERROR_IF(zenith_std <= 0.0,
                     "Gaussian antenna zenith_std must be positive")
  ARTS_USER_ERROR_IF(azimuth_std <= 0.0,
                     "Gaussian antenna azimuth_std must be positive")

  AntennaPatternField out{
      .data_name  = "gaussian"s,
      .data       = StokvecMatrix(zen_grid.size(), azi_grid.size()),
      .grid_names = {"zenith"s, "azimuth"s},
      .grids      = {std::move(zen_grid), std::move(azi_grid)},
  };

  using Conversion::atan2d;

  for (Size izen = 0; izen < out.grid<0>().size(); ++izen) {
    for (Size iazi = 0; iazi < out.grid<1>().size(); ++iazi) {
      const Vector3 local =
          antenna_frame_los({out.grid<0>()[izen], out.grid<1>()[iazi]});

      if (local[2] <= 0.0) {
        out[izen, iazi] = {0.0, 0.0, 0.0, 0.0};
        continue;
      }

      const Numeric ant_zen = atan2d(local[0], local[2]);
      const Numeric ant_azi = atan2d(local[1], local[2]);
      const Numeric exponent =
          -0.5 * ((ant_zen / zenith_std) * (ant_zen / zenith_std) +
                  (ant_azi / azimuth_std) * (ant_azi / azimuth_std));

      out[izen, iazi] = std::exp(exponent) * weight;
    }
  }

  return out;
}
}  // namespace

PencilBeamAntenna::PencilBeamAntenna(Stokvec weight)
    : AntennaPattern({.data_name  = "pencil beam"s,
                      .data       = StokvecMatrix(1, 1, weight),
                      .grid_names = {"zenith"s, "azimuth"s},
                      .grids      = {Vector{0.0}, Vector{0.0}}}) {}

GaussianAntenna::GaussianAntenna(ZenGrid zen_grid,
                                 AziGrid azi_grid,
                                 Numeric zenith_std,
                                 Numeric azimuth_std,
                                 Stokvec weight)
    : AntennaPattern(make_gaussian_field(std::move(zen_grid),
                                         std::move(azi_grid),
                                         zenith_std,
                                         azimuth_std,
                                         weight)) {}

std::vector<std::pair<Stokvec, Vector2>> AntennaPattern::operator()(
    Vector2 bore_los) const {
  ARTS_USER_ERROR_IF(not data.ok(),
                     "SensorAntennaPattern data shape does not match its grids")

  const auto& zen_grid = data.grid<0>();
  const auto& azi_grid = data.grid<1>();

  std::vector<std::pair<Stokvec, Vector2>> out;
  out.reserve(static_cast<std::size_t>(zen_grid.size()) *
              static_cast<std::size_t>(azi_grid.size()));

  const auto basis = antenna_basis(bore_los);

  for (Size izen = 0; izen < zen_grid.size(); ++izen) {
    for (Size iazi = 0; iazi < azi_grid.size(); ++iazi) {
      const Vector3 local = antenna_frame_los({zen_grid[izen], azi_grid[iazi]});
      const Vector3 enu   = normalized(local[0] * basis.v + local[1] * basis.h +
                                       local[2] * basis.k);
      out.emplace_back(data[izen, iazi], enu2los(enu));
    }
  }

  return out;
}

static_assert(AntennaPatternSelection<AntennaPattern>);
static_assert(AntennaPatternSelection<PencilBeamAntenna>);
static_assert(AntennaPatternSelection<GaussianAntenna>);
}  // namespace sensor