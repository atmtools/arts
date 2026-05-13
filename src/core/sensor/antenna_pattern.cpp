#include "antenna_pattern.h"

#include <arts_constants.h>
#include <arts_conversions.h>
#include <geodetic.h>

#include <cmath>
#include <memory>

namespace sensor {

namespace {
constexpr Numeric gaussian_airy_hwhm_factor =
    3.8317059702075123156 / Constant::pi;

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

[[nodiscard]] Numeric gaussian_airy_std(Numeric frequency,
                                        Numeric aperture_diameter) {
  const Numeric wavelength = Constant::speed_of_light / frequency;
  const Numeric hwhm_deg = Conversion::rad2deg(gaussian_airy_hwhm_factor *
                                               wavelength / aperture_diameter);

  return Conversion::hwhm2std(hwhm_deg);
}

std::shared_ptr<const PosLosVector> make_single_poslos_grid(
    const Vector3& pos, const Vector2& los) {
  return std::make_shared<PosLosVector>(1, PosLos{.pos = pos, .los = los});
}
}  // namespace

PencilBeamAntenna::PencilBeamAntenna(Stokvec weight) : weight(weight) {}

Obsel PencilBeamAntenna::operator()(const Channel& channel,
                                    const Vector3& pos,
                                    const Vector2& bore_los) const {
  const auto& channel_weights = channel.weights();
  const auto freq_grid =
      std::make_shared<const AscendingGrid>(channel.freq_grid());
  SparseStokvecMatrix weight_matrix(1, channel_weights.size());

  if (not weight.is_zero()) {
    for (Size ifreq = 0; ifreq < channel_weights.size(); ++ifreq) {
      if (channel_weights[ifreq] == 0.0) continue;
      weight_matrix[0, ifreq] = channel_weights[ifreq] * weight;
    }
  }

  return {freq_grid,
          make_single_poslos_grid(pos, bore_los),
          std::move(weight_matrix)};
}

std::shared_ptr<const AntennaPattern> PencilBeamAntenna::clone() const {
  return std::make_shared<PencilBeamAntenna>(*this);
}

Obsel GriddedAntennaPattern::operator()(const Channel& channel,
                                        const Vector3& pos,
                                        const Vector2& bore_los) const {
  ARTS_USER_ERROR_IF(
      not data.ok(),
      "SensorGriddedAntennaPattern data shape does not match its grids")

  const auto& zen_grid        = data.grid<0>();
  const auto& azi_grid        = data.grid<1>();
  const auto& channel_weights = channel.weights();
  const auto freq_grid =
      std::make_shared<const AscendingGrid>(channel.freq_grid());

  const Size nsamples = zen_grid.size() * azi_grid.size();
  auto poslos_grid    = std::make_shared<PosLosVector>(nsamples);
  SparseStokvecMatrix weight_matrix(nsamples, channel_weights.size());

  const auto basis = antenna_basis(bore_los);
  Size isample     = 0;

  for (Size izen = 0; izen < zen_grid.size(); ++izen) {
    for (Size iazi = 0; iazi < azi_grid.size(); ++iazi) {
      const Vector3 local = antenna_frame_los({zen_grid[izen], azi_grid[iazi]});
      const Vector3 enu   = normalized(local[0] * basis.v + local[1] * basis.h +
                                       local[2] * basis.k);
      (*poslos_grid)[isample] = {.pos = pos, .los = enu2los(enu)};

      const auto& antenna_weight = data[izen, iazi];
      if (not antenna_weight.is_zero()) {
        for (Size ifreq = 0; ifreq < channel_weights.size(); ++ifreq) {
          if (channel_weights[ifreq] == 0.0) continue;
          weight_matrix[isample, ifreq] =
              channel_weights[ifreq] * antenna_weight;
        }
      }

      ++isample;
    }
  }

  return {freq_grid, poslos_grid, std::move(weight_matrix)};
}

std::shared_ptr<const AntennaPattern> GriddedAntennaPattern::clone() const {
  return std::make_shared<GriddedAntennaPattern>(*this);
}

GaussianAntenna::GaussianAntenna(ZenGrid zen_grid,
                                 AziGrid azi_grid,
                                 Numeric zenith_std,
                                 Numeric azimuth_std,
                                 Stokvec weight) {
  data = make_gaussian_field(std::move(zen_grid),
                             std::move(azi_grid),
                             zenith_std,
                             azimuth_std,
                             weight);
}

GaussianAntenna::GaussianAntenna(ZenGrid zen_grid,
                                 AziGrid azi_grid,
                                 Numeric std,
                                 Stokvec weight)
    : GaussianAntenna(
          std::move(zen_grid), std::move(azi_grid), std, std, weight) {}

std::shared_ptr<const AntennaPattern> GaussianAntenna::clone() const {
  return std::make_shared<GaussianAntenna>(*this);
}

GaussianAiryAntenna::GaussianAiryAntenna(ZenGrid zen_grid,
                                         AziGrid azi_grid,
                                         Numeric aperture_diameter,
                                         Stokvec weight)
    : zen_grid(std::move(zen_grid)),
      azi_grid(std::move(azi_grid)),
      aperture_diameter(aperture_diameter),
      weight(weight) {}

Obsel GaussianAiryAntenna::operator()(const Channel& channel,
                                      const Vector3& pos,
                                      const Vector2& bore_los) const {
  ARTS_USER_ERROR_IF(aperture_diameter <= 0.0,
                     "Gaussian Airy antenna aperture_diameter must be positive")

  const auto& channel_weights = channel.weights();
  const auto freq_grid =
      std::make_shared<const AscendingGrid>(channel.freq_grid());

  ARTS_USER_ERROR_IF(
      freq_grid->front() <= 0.0,
      "Gaussian Airy antenna requires strictly positive channel frequencies because SensorBuilder only provides the channel frequency grid")

  const Size nsamples = zen_grid.size() * azi_grid.size();
  auto poslos_grid    = std::make_shared<PosLosVector>(nsamples);
  SparseStokvecMatrix weight_matrix(nsamples, channel_weights.size());

  const auto basis = antenna_basis(bore_los);
  Size isample     = 0;

  for (double izen : zen_grid) {
    for (double iazi : azi_grid) {
      const Vector3 local = antenna_frame_los({izen, iazi});
      const Vector3 enu   = normalized(local[0] * basis.v + local[1] * basis.h +
                                       local[2] * basis.k);
      (*poslos_grid)[isample] = {.pos = pos, .los = enu2los(enu)};

      ++isample;
    }
  }

  if (not weight.is_zero()) {
    for (Size ifreq = 0; ifreq < channel_weights.size(); ++ifreq) {
      if (channel_weights[ifreq] == 0.0) continue;

      const Numeric airy_std =
          gaussian_airy_std((*freq_grid)[ifreq], aperture_diameter);
      const GaussianAntenna frequency_pattern{
          zen_grid, azi_grid, airy_std, weight};

      isample = 0;
      for (Size izen = 0; izen < zen_grid.size(); ++izen) {
        for (Size iazi = 0; iazi < azi_grid.size(); ++iazi) {
          const auto& antenna_weight = frequency_pattern.data[izen, iazi];
          if (not antenna_weight.is_zero()) {
            weight_matrix[isample, ifreq] =
                channel_weights[ifreq] * antenna_weight;
          }

          ++isample;
        }
      }
    }
  }

  return {freq_grid, poslos_grid, std::move(weight_matrix)};
}

std::shared_ptr<const AntennaPattern> GaussianAiryAntenna::clone() const {
  return std::make_shared<GaussianAiryAntenna>(*this);
}

static_assert(AntennaPatternSelection<AntennaPattern>);
static_assert(AntennaPatternSelection<GriddedAntennaPattern>);
static_assert(AntennaPatternSelection<PencilBeamAntenna>);
static_assert(AntennaPatternSelection<GaussianAntenna>);
static_assert(AntennaPatternSelection<GaussianAiryAntenna>);
}  // namespace sensor
