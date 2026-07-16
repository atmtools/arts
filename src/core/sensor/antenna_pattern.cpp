#include "antenna_pattern.h"

#include <arts_constants.h>
#include <arts_conversions.h>
#include <geodetic.h>

#include <cmath>
#include <memory>
#include <ranges>
#include <vector>

namespace sensor {

namespace {
constexpr Numeric gaussian_airy_hwhm_factor = Constant::bessel_j_n1_k1_zero / Constant::pi;

struct AntennaGeometrySample {
  PosLos  poslos;
  Numeric weight;
  Vector2 local_los;
};

template <typename ZenGridT, typename AziGridT>
[[nodiscard]] std::vector<AntennaGeometrySample> make_antenna_geometry_layout(const Matrix&   antenna_data,
                                                                              const ZenGridT& zen_grid,
                                                                              const AziGridT& azi_grid,
                                                                              const Vector3&  pos,
                                                                              const Vector2&  bore_los,
                                                                              const Vector2&  ell) {
  const auto [pos_ecef, los_ecef] = geodetic_los2ecef(pos, bore_los, ell);

  auto rot = [](Vector3 v, Vector3 k, Numeric aa) {
    const Numeric cos_azi = Conversion::cosd(aa);
    const Numeric sin_azi = Conversion::sind(aa);
    return v * cos_azi + cross(k, v) * sin_azi + k * dot(k, v) * (1.0 - cos_azi);
  };

  const std::vector<Vector3> zen_los{std::from_range, zen_grid | stdv::transform([&](Numeric zen) {
                                                        Vector2 bore_los2    = bore_los;
                                                        bore_los2[0]        += zen;
                                                        const auto [_, los]  = geodetic_los2ecef(pos, bore_los2, ell);
                                                        return los;
                                                      })};

  std::vector<AntennaGeometrySample> out{};
  out.reserve(zen_grid.size() * azi_grid.size());
  Numeric renorm{};
  for (Size izen = 0; izen < zen_grid.size(); ++izen) {
    if (zen_grid[izen] == 0.0) {
      const auto su  = sum(antenna_data[izen]);
      const auto mea = su / static_cast<Numeric>(azi_grid.size());
      out.emplace_back(
          AntennaGeometrySample{.poslos = {.pos = pos, .los = bore_los}, .weight = mea, .local_los = {0.0, 0.0}});
      renorm = su - mea;
      continue;
    }

    for (Size iazi = 0; iazi < azi_grid.size(); ++iazi) {
      auto enu      = rot(zen_los[izen], los_ecef, azi_grid[iazi]);
      auto [_, los] = ecef2geodetic_los(pos_ecef, enu, ell);
      out.emplace_back(AntennaGeometrySample{.poslos    = {.pos = pos, .los = los},
                                             .weight    = antenna_data[izen, iazi],
                                             .local_los = {zen_grid[izen], azi_grid[iazi]}});
    }
  }

  if (renorm != 0.0) {
    const Numeric v = 1.0 / (1.0 - renorm);
    for (auto& w : out) w.weight *= v;
  }

  return out;
}

[[nodiscard]] AntennaPatternField make_gaussian_field(ZenGrid zen_grid, AziGrid azi_grid, Numeric std) {
  ARTS_USER_ERROR_IF(std <= 0.0, "Gaussian antenna std must be positive")

  AntennaPatternField out{
      .data_name  = "gaussian"s,
      .data       = Matrix(zen_grid.size(), azi_grid.size()),
      .grid_names = {"zenith"s, "azimuth"s},
      .grids      = {std::move(zen_grid), std::move(azi_grid)},
  };

  for (Size izen = 0; izen < out.grid<0>().size(); ++izen) {
    out[izen] = std::exp(-0.5 * Math::pow2(out.grid<0>()[izen] / std));
  }

  if (auto su = sum(out.data); su > 0.0) out.data /= su;

  return out;
}

[[nodiscard]] Numeric gaussian_airy_std(Numeric frequency, Numeric aperture_diameter) {
  using Conversion::rad2deg, Conversion::hwhm2std, Conversion::freq2wavelen;
  const Numeric hwhm_rad = gaussian_airy_hwhm_factor * freq2wavelen(frequency) / aperture_diameter;

  return rad2deg(hwhm2std(hwhm_rad));
}

[[nodiscard]] Numeric gaussian_airy_response(Numeric offset_zenith, Numeric airy_std) {
  return std::exp(-0.5 * Math::pow2(offset_zenith / airy_std));
}

std::shared_ptr<const PosLosVector> make_single_poslos_grid(const Vector3& pos, const Vector2& los) {
  return std::make_shared<PosLosVector>(1, PosLos{.pos = pos, .los = los});
}
}  // namespace

PencilBeamAntenna::PencilBeamAntenna(Stokvec weight) : weight(weight) {}

Obsel PencilBeamAntenna::operator()(const Channel& channel,
                                    const Vector3& pos,
                                    const Vector2& bore_los,
                                    const Vector2&) const {
  const auto&         channel_weights = channel.weights();
  const auto          freq_grid       = std::make_shared<const AscendingGrid>(channel.freq_grid());
  SparseStokvecMatrix weight_matrix(1, channel_weights.size());

  if (not weight.is_zero()) {
    for (Size ifreq = 0; ifreq < channel_weights.size(); ++ifreq) {
      if (channel_weights[ifreq] == 0.0) continue;
      weight_matrix[0, ifreq] = channel_weights[ifreq] * weight;
    }
  }

  return {freq_grid, make_single_poslos_grid(pos, bore_los), std::move(weight_matrix)};
}

std::shared_ptr<const AntennaPattern> PencilBeamAntenna::clone() const {
  return std::make_shared<PencilBeamAntenna>(*this);
}

Obsel GriddedAntennaPattern::operator()(const Channel& channel,
                                        const Vector3& pos,
                                        const Vector2& bore_los,
                                        const Vector2& ell) const {
  ARTS_USER_ERROR_IF(not data.ok(), "GriddedAntennaPattern data shape does not match its grids")

  const auto& zen_grid        = data.grid<0>();
  const auto& azi_grid        = data.grid<1>();
  const auto& channel_weights = channel.weights();
  const auto  freq_grid       = std::make_shared<const AscendingGrid>(channel.freq_grid());

  const auto antenna_geom = make_antenna_geometry_layout(data.data, zen_grid, azi_grid, pos, bore_los, ell);

  SparseStokvecMatrix weight_matrix(antenna_geom.size(), channel_weights.size());

  for (Size i = 0; i < antenna_geom.size(); ++i) {
    if (antenna_geom[i].weight == 0.0) continue;

    for (Size j = 0; j < channel.size(); j++) {
      if (channel_weights[j] == 0.0) continue;

      weight_matrix[i, j] = channel_weights[j] * antenna_geom[i].weight * weight;
    }
  }

  return {freq_grid,
          std::make_shared<PosLosVector>(std::from_range,
                                         antenna_geom | stdv::transform([](const auto& s) { return s.poslos; })),
          std::move(weight_matrix)};
}

std::shared_ptr<const AntennaPattern> GriddedAntennaPattern::clone() const {
  return std::make_shared<GriddedAntennaPattern>(*this);
}

namespace {
AziGrid make_azi_grid(Size size) {
  ARTS_USER_ERROR_IF(size < 1, "Gaussian antenna azimuth size must be at least 1")

  return Vector{nlinspace(0.0, 360.0, size + 1)[Range(0, size)]};
}
}  // namespace

GaussianAntenna::GaussianAntenna(ZenGrid zen_grid, Numeric std, Size azi_grid_size, Stokvec weight_) {
  data   = make_gaussian_field(std::move(zen_grid), make_azi_grid(azi_grid_size), std);
  weight = weight_;
}

std::shared_ptr<const AntennaPattern> GaussianAntenna::clone() const {
  return std::make_shared<GaussianAntenna>(*this);
}

GaussianAiryAntenna::GaussianAiryAntenna(ZenGrid zen_grid,
                                         Numeric aperture_diameter,
                                         Size    azi_grid_size,
                                         Stokvec weight)
    : zen_grid(std::move(zen_grid)),
      azi_grid(make_azi_grid(azi_grid_size)),
      aperture_diameter(aperture_diameter),
      weight(weight) {}

Numeric GaussianAiryAntenna::std(Numeric frequency) const { return gaussian_airy_std(frequency, aperture_diameter); }

Obsel GaussianAiryAntenna::operator()(const Channel& channel,
                                      const Vector3& pos,
                                      const Vector2& bore_los,
                                      const Vector2& ell) const {
  ARTS_USER_ERROR_IF(aperture_diameter <= 0.0, "Gaussian Airy antenna aperture_diameter must be positive")

  const auto& channel_weights = channel.weights();
  const auto  freq_grid       = std::make_shared<const AscendingGrid>(channel.freq_grid());

  ARTS_USER_ERROR_IF(
      freq_grid->front() <= 0.0,
      "Gaussian Airy antenna requires strictly positive channel frequencies because sensor builder only provides the channel frequency grid")

  Matrix data(zen_grid.size(), azi_grid.size(), 1.0);

  const auto antenna_geom = make_antenna_geometry_layout(data, zen_grid, azi_grid, pos, bore_los, ell);

  Matrix airy_ws(channel.freq_grid().size(), antenna_geom.size(), 0.0);
  for (Size ifreq = 0; ifreq < channel.freq_grid().size(); ++ifreq) {
    const Numeric airy_std = std(channel.freq_grid()[ifreq]);
    for (Size ipos = 0; ipos < antenna_geom.size(); ++ipos) {
      const Numeric resp   = gaussian_airy_response(antenna_geom[ipos].local_los[0], airy_std) * channel_weights[ifreq];
      airy_ws[ifreq, ipos] = resp;
    }
  }
  airy_ws /= sum(airy_ws);

  SparseStokvecMatrix weight_matrix(antenna_geom.size(), channel_weights.size());

  for (Size i = 0; i < antenna_geom.size(); ++i) {
    for (Size j = 0; j < channel.size(); j++) {
      if (channel_weights[j] == 0.0) continue;

      weight_matrix[i, j] = airy_ws[j, i] * weight;
    }
  }

  return {freq_grid,
          std::make_shared<PosLosVector>(std::from_range,
                                         antenna_geom | stdv::transform([](const auto& s) { return s.poslos; })),
          std::move(weight_matrix)};
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

void xml_io_stream<sensor::GriddedAntennaPattern>::write(std::ostream&                        os,
                                                         const sensor::GriddedAntennaPattern& n,
                                                         bofstream*                           pbofs,
                                                         std::string_view                     name) {
  XMLTag tag(xml_io_stream_name_v<sensor::GriddedAntennaPattern>, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, n.data, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<sensor::GriddedAntennaPattern>::read(std::istream&                  is,
                                                        sensor::GriddedAntennaPattern& n,
                                                        bifstream*                     pbifs) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(xml_io_stream_name_v<sensor::GriddedAntennaPattern>);

  xml_read_from_stream(is, n.data, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(xml_io_stream_name_v<sensor::GriddedAntennaPattern>);
}

void xml_io_stream<sensor::GaussianAiryAntenna>::write(std::ostream&                      os,
                                                       const sensor::GaussianAiryAntenna& n,
                                                       bofstream*                         pbofs,
                                                       std::string_view                   name) {
  XMLTag tag(xml_io_stream_name_v<sensor::GaussianAiryAntenna>, "name", name);
  tag.write_to_stream(os);
  xml_write_to_stream(os, n.aperture_diameter, pbofs);
  xml_write_to_stream(os, n.zen_grid, pbofs);
  xml_write_to_stream(os, n.azi_grid, pbofs);
  xml_write_to_stream(os, n.weight, pbofs);
  tag.write_to_end_stream(os);
}

void xml_io_stream<sensor::GaussianAiryAntenna>::read(std::istream&                is,
                                                      sensor::GaussianAiryAntenna& n,
                                                      bifstream*                   pbifs) {
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

void xml_io_stream<sensor::PencilBeamAntenna>::write(std::ostream&                    os,
                                                     const sensor::PencilBeamAntenna& n,
                                                     bofstream*                       pbofs,
                                                     std::string_view                 name) {
  XMLTag tag(xml_io_stream_name_v<sensor::PencilBeamAntenna>, "name", name);
  tag.write_to_stream(os);
  xml_write_to_stream(os, n.weight, pbofs);
  tag.write_to_end_stream(os);
}

void xml_io_stream<sensor::PencilBeamAntenna>::read(std::istream& is, sensor::PencilBeamAntenna& n, bifstream* pbifs) {
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

void xml_io_stream<sensor::AntennaPattern>::read(std::istream& is, sensor::AntennaPattern&, bifstream*) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(xml_io_stream_name_v<sensor::AntennaPattern>);
  tag.read_from_stream(is);
  tag.check_end_name(xml_io_stream_name_v<sensor::AntennaPattern>);
}
