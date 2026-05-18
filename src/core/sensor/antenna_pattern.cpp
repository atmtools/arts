#include "antenna_pattern.h"

#include <arts_constants.h>
#include <arts_conversions.h>
#include <geodetic.h>

#include <cmath>
#include <limits>
#include <memory>
#include <vector>

namespace sensor {

namespace {
constexpr Numeric gaussian_airy_hwhm_factor =
    3.8317059702075123156 / Constant::pi;

struct AntennaBasis {
  Vector3 v;
  Vector3 h;
  Vector3 k;
};

struct AntennaGeometrySample {
  Size sample_index{};
  Numeric ant_zen{};
  Numeric ant_azi{};
  bool has_response{};
  bool is_representative{};
};

struct AntennaGeometryLayout {
  std::shared_ptr<PosLosVector> poslos_grid;
  std::vector<AntennaGeometrySample> samples;
};

constexpr Numeric azimuth_degenerate_tol = 1e-12;

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

[[nodiscard]] bool antenna_azimuth_is_degenerate(const Vector3& local) {
  return std::hypot(local[0], local[1]) <= azimuth_degenerate_tol;
}

template <typename ZenGridT, typename AziGridT>
[[nodiscard]] AntennaGeometryLayout make_antenna_geometry_layout(
    const ZenGridT& zen_grid,
    const AziGridT& azi_grid,
    const Vector3& pos,
    const Vector2& bore_los) {
  const auto basis = antenna_basis(bore_los);

  std::vector<PosLos> poslos;
  poslos.reserve(zen_grid.size() * azi_grid.size());

  std::vector<AntennaGeometrySample> samples;
  samples.reserve(zen_grid.size() * azi_grid.size());

  using Conversion::atan2d;

  for (Numeric izen : zen_grid) {
    Size degenerate_sample_index = std::numeric_limits<Size>::max();

    for (Numeric iazi : azi_grid) {
      const Vector3 local   = antenna_frame_los({izen, iazi});
      const bool degenerate = antenna_azimuth_is_degenerate(local);

      Size sample_index = degenerate_sample_index;
      bool is_representative = false;
      if (sample_index == std::numeric_limits<Size>::max()) {
        const Vector3 enu = normalized(local[0] * basis.v + local[1] * basis.h +
                                       local[2] * basis.k);
        sample_index      = poslos.size();
        poslos.push_back({.pos = pos, .los = enu2los(enu)});
        is_representative = true;

        if (degenerate) degenerate_sample_index = sample_index;
      }

      auto& sample = samples.emplace_back(AntennaGeometrySample{
          .sample_index      = sample_index,
          .has_response      = local[2] > 0.0,
          .is_representative = is_representative});

      if (sample.has_response) {
        sample.ant_zen = atan2d(local[0], local[2]);
        sample.ant_azi = atan2d(local[1], local[2]);
      }
    }
  }

  auto poslos_grid = std::make_shared<PosLosVector>(poslos.size());
  for (Size isample = 0; isample < poslos.size(); ++isample) {
    (*poslos_grid)[isample] = poslos[isample];
  }

  return {std::move(poslos_grid), std::move(samples)};
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

[[nodiscard]] Numeric gaussian_airy_response(Numeric ant_zen,
                                             Numeric ant_azi,
                                             Numeric inv_std_sq) {
  return std::exp(-0.5 * (ant_zen * ant_zen + ant_azi * ant_azi) * inv_std_sq);
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

  auto geometry =
      make_antenna_geometry_layout(zen_grid, azi_grid, pos, bore_los);
  SparseStokvecMatrix weight_matrix(geometry.poslos_grid->size(),
                                    channel_weights.size());

  Size igrid = 0;

  for (Size izen = 0; izen < zen_grid.size(); ++izen) {
    for (Size iazi = 0; iazi < azi_grid.size(); ++iazi) {
      const auto& antenna_weight = data[izen, iazi];
      if (not antenna_weight.is_zero()) {
        const Size sample_index = geometry.samples[igrid].sample_index;
        for (Size ifreq = 0; ifreq < channel_weights.size(); ++ifreq) {
          if (channel_weights[ifreq] == 0.0) continue;
          weight_matrix[sample_index, ifreq] +=
              channel_weights[ifreq] * antenna_weight;
        }
      }

      ++igrid;
    }
  }

  return {freq_grid, std::move(geometry.poslos_grid), std::move(weight_matrix)};
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

  auto geometry =
      make_antenna_geometry_layout(zen_grid, azi_grid, pos, bore_los);
  SparseStokvecMatrix weight_matrix(geometry.poslos_grid->size(),
                                    channel_weights.size());

  if (not weight.is_zero()) {
    Vector inv_std_sq(channel_weights.size(), 0.0);
    Vector normalization(channel_weights.size(), 0.0);
    Vector retained_scale(channel_weights.size(), 0.0);

    for (Size ifreq = 0; ifreq < channel_weights.size(); ++ifreq) {
      if (channel_weights[ifreq] == 0.0) continue;

      const Numeric airy_std =
          gaussian_airy_std((*freq_grid)[ifreq], aperture_diameter);
      inv_std_sq[ifreq] = 1.0 / (airy_std * airy_std);

      Numeric retained_mass = 0.0;
      for (const auto& sample : geometry.samples) {
        if (not sample.has_response) continue;
        normalization[ifreq] += gaussian_airy_response(
            sample.ant_zen, sample.ant_azi, inv_std_sq[ifreq]);
      }

      if (normalization[ifreq] == 0.0) continue;

      for (const auto& sample : geometry.samples) {
        if (not sample.has_response) continue;
        if (not sample.is_representative) continue;

        const Numeric single_weight =
            channel_weights[ifreq] * gaussian_airy_response(
                                        sample.ant_zen,
                                        sample.ant_azi,
                                        inv_std_sq[ifreq]) /
            normalization[ifreq];

        retained_mass += single_weight;
      }

      if (retained_mass > 0.0) {
        retained_scale[ifreq] = channel_weights[ifreq] / retained_mass;
      }
    }

    // Keep sparse insertions row-major so they append instead of shifting.
    for (const auto& sample : geometry.samples) {
      if (not sample.has_response) continue;
      if (not sample.is_representative) continue;

      for (Size ifreq = 0; ifreq < channel_weights.size(); ++ifreq) {
        if (normalization[ifreq] == 0.0) continue;
        if (retained_scale[ifreq] == 0.0) continue;

        const Numeric single_weight =
            channel_weights[ifreq] * gaussian_airy_response(
                                        sample.ant_zen,
                                        sample.ant_azi,
                                        inv_std_sq[ifreq]) /
            normalization[ifreq];
        if (single_weight == 0.0) continue;

        weight_matrix[sample.sample_index, ifreq] +=
            (retained_scale[ifreq] * single_weight) * weight;
      }
    }
  }

  return {freq_grid, std::move(geometry.poslos_grid), std::move(weight_matrix)};
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

void xml_io_stream<sensor::GriddedAntennaPattern>::write(
    std::ostream& os,
    const sensor::GriddedAntennaPattern& n,
    bofstream* pbofs,
    std::string_view name) {
  XMLTag tag(xml_io_stream_name_v<sensor::GriddedAntennaPattern>, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, n.data, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<sensor::GriddedAntennaPattern>::read(
    std::istream& is, sensor::GriddedAntennaPattern& n, bifstream* pbifs) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(xml_io_stream_name_v<sensor::GriddedAntennaPattern>);

  xml_read_from_stream(is, n.data, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(xml_io_stream_name_v<sensor::GriddedAntennaPattern>);
}

void xml_io_stream<sensor::GaussianAiryAntenna>::write(
    std::ostream& os,
    const sensor::GaussianAiryAntenna& n,
    bofstream* pbofs,
    std::string_view name) {
  XMLTag tag(xml_io_stream_name_v<sensor::GaussianAiryAntenna>, "name", name);
  tag.write_to_stream(os);
  xml_write_to_stream(os, n.aperture_diameter, pbofs);
  xml_write_to_stream(os, n.zen_grid, pbofs);
  xml_write_to_stream(os, n.azi_grid, pbofs);
  xml_write_to_stream(os, n.weight, pbofs);
  tag.write_to_end_stream(os);
}

void xml_io_stream<sensor::GaussianAiryAntenna>::read(
    std::istream& is, sensor::GaussianAiryAntenna& n, bifstream* pbifs) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(xml_io_stream_name_v<sensor::GaussianAiryAntenna>);

  xml_read_from_stream(is, n.aperture_diameter, pbifs);
  xml_read_from_stream(is, n.zen_grid, pbifs);
  xml_read_from_stream(is, n.azi_grid, pbifs);
  xml_read_from_stream(is, n.weight, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(xml_io_stream_name_v<sensor::GaussianAiryAntenna>);
}

void xml_io_stream<sensor::PencilBeamAntenna>::write(
    std::ostream& os,
    const sensor::PencilBeamAntenna& n,
    bofstream* pbofs,
    std::string_view name) {
  XMLTag tag(xml_io_stream_name_v<sensor::PencilBeamAntenna>, "name", name);
  tag.write_to_stream(os);
  xml_write_to_stream(os, n.weight, pbofs);
  tag.write_to_end_stream(os);
}

void xml_io_stream<sensor::PencilBeamAntenna>::read(
    std::istream& is, sensor::PencilBeamAntenna& n, bifstream* pbifs) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(xml_io_stream_name_v<sensor::PencilBeamAntenna>);
  xml_read_from_stream(is, n.weight, pbifs);
  tag.read_from_stream(is);
  tag.check_end_name(xml_io_stream_name_v<sensor::PencilBeamAntenna>);
}

void xml_io_stream<sensor::AntennaPattern>::write(std::ostream& os,
                                                  const sensor::AntennaPattern&,
                                                  bofstream*,
                                                  std::string_view name) {
  XMLTag tag(xml_io_stream_name_v<sensor::AntennaPattern>, "name", name);
  tag.write_to_stream(os);
  tag.write_to_end_stream(os);
}

void xml_io_stream<sensor::AntennaPattern>::read(std::istream& is,
                                                 sensor::AntennaPattern&,
                                                 bifstream*) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(xml_io_stream_name_v<sensor::AntennaPattern>);
  tag.read_from_stream(is);
  tag.check_end_name(xml_io_stream_name_v<sensor::AntennaPattern>);
}
