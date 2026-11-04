#include "antenna_pattern.h"

#include <arts_constants.h>
#include <arts_conversions.h>
#include <geodetic.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <ranges>
#include <vector>

namespace sensor {

namespace {
//! The half-power half width of the Airy pattern of a circular aperture is this times wavelength over diameter [rad]
constexpr Numeric gaussian_airy_hwhm_factor = Constant::airy_half_power_argument / Constant::pi;

//! The unit vector of a line of sight in the local east-north-up frame
Vector3 enu_direction(const Vector2& los) {
  using Conversion::cosd, Conversion::sind;
  return {sind(los[0]) * sind(los[1]), sind(los[0]) * cosd(los[1]), cosd(los[0])};
}

//! The polarization basis of the radiation received along a line of sight: that of its propagation direction
rtepack::polarization_basis received_basis(const Vector2& los) {
  using Conversion::cosd, Conversion::sind;
  const Numeric propagation_azimuth = los[1] + 180.0;
  return rtepack::meridional_basis(enu_direction(los) * -1.0,
                                   Vector3{cosd(propagation_azimuth), -sind(propagation_azimuth), 0.0});
}

//! v turned by angle [deg] about the unit axis k, right-handed
Vector3 turn(const Vector3& v, const Vector3& k, Numeric angle) {
  const Numeric c = Conversion::cosd(angle);
  const Numeric s = Conversion::sind(angle);
  return v * c + cross(k, v) * s + k * dot(k, v) * (1.0 - c);
}

/** The Stokes change of basis from a basis carried along with the antenna to that of the line of sight it sampled
 *
 *  The two propagation directions agree to the rounding of the geodetic conversions, except within the pole
 *  threshold of those conversions, where the stored azimuth of the line of sight is set to 0 and moves its direction
 *  by up to that threshold.  The radiative transfer of the stored line of sight uses its basis, so the carried basis
 *  is projected on it.
 */
Muelmat basis_change(const rtepack::polarization_basis& carried, const rtepack::polarization_basis& sampled) {
  const auto [c2, s2] = rtepack::detail::rotation_angle(carried, sampled);
  const Numeric norm  = std::hypot(c2, s2);
  return rtepack::stokes_rotation(c2 / norm, s2 / norm);
}

/** The change from the boresight's polarization basis to that of a line of sight
 *
 *  The basis is carried along the great circle from the boresight to the direction of the line of sight, the
 *  shortest rotation between them, so it does not spin about the boresight.  Directions are unit vectors in the
 *  local east-north-up frame.
 */
Muelmat carried_basis_change(const Vector3&                     bore,
                             const rtepack::polarization_basis& bore_basis,
                             const Vector3&                     direction,
                             const Vector2&                     los) {
  const Vector3 normal = cross(bore, direction);
  const Numeric sine   = hypot(normal);
  const Numeric angle  = Conversion::atan2d(sine, dot(bore, direction));
  const auto    carry  = [&](const Vector3& v) { return sine > 0.0 ? turn(v, normal / sine, angle) : v; };
  return basis_change({carry(bore_basis.first), carry(bore_basis.second)}, received_basis(los));
}

struct AntennaGeometrySample {
  PosLos poslos;
  //! The pattern's gain at the sample, its azimuthal mean at zero zenith offset
  Numeric gain;
  //! The solid angle [sr] of the cell of the sphere that the sample stands for
  Numeric solid_angle;
  //! [zenith, azimuth] offset [deg] of the sample from the boresight
  Vector2 local_los;
  //! The change from the boresight's polarization basis, carried with the antenna, to the sample's
  Muelmat basis_from_boresight;
};

//! The zenith offsets [deg] bounding the ring of each zenith sample: midpoints, clamped to [0, 180]
std::vector<Vector2> ring_bounds(const ZenGrid& zen) {
  const Size           n = zen.size();
  std::vector<Vector2> out(n, Vector2{zen[0], zen[0]});
  if (n < 2) return out;

  for (Size i = 0; i < n; i++) {
    out[i][0] = i == 0 ? std::max(0.0, zen[0] - 0.5 * (zen[1] - zen[0])) : 0.5 * (zen[i - 1] + zen[i]);
    out[i][1] = i + 1 == n ? std::min(180.0, zen[i] + 0.5 * (zen[i] - zen[i - 1])) : 0.5 * (zen[i] + zen[i + 1]);
  }

  return out;
}

//! The fraction of the full circle that each azimuth sample stands for: midpoints, periodic over 360 degrees
std::vector<Numeric> azimuth_fractions(const AziGrid& azi) {
  const Size           n = azi.size();
  std::vector<Numeric> out(n, 1.0);
  if (n < 2) return out;

  for (Size j = 0; j < n; j++) {
    const Numeric prev = j == 0 ? azi[n - 1] - 360.0 : azi[j - 1];
    const Numeric next = j + 1 == n ? azi[0] + 360.0 : azi[j + 1];
    out[j]             = 0.5 * (next - prev) / 360.0;
  }

  return out;
}

[[nodiscard]] std::vector<AntennaGeometrySample> make_antenna_geometry_layout(const Matrix&  gain,
                                                                              const ZenGrid& zen_grid,
                                                                              const AziGrid& azi_grid,
                                                                              const Vector3& pos,
                                                                              const Vector2& bore_los,
                                                                              const Vector2& ell) {
  using Conversion::cosd, Conversion::sind;

  if (zen_grid.empty() or azi_grid.empty()) arts::user_error("An antenna pattern needs zenith and azimuth offsets");

  const auto [pos_ecef, los_ecef] = geodetic_los2ecef(pos, bore_los, ell);
  const auto bounds               = ring_bounds(zen_grid);
  const auto azi_fraction         = azimuth_fractions(azi_grid);

  // In the local east-north-up frame: the boresight, the axis that turns it towards larger zenith angles, and the
  // polarization basis of the radiation it receives
  const Vector3 bore       = enu_direction(bore_los);
  const Vector3 tilt_axis  = {-cosd(bore_los[1]), sind(bore_los[1]), 0.0};
  const auto    bore_basis = received_basis(bore_los);

  std::vector<AntennaGeometrySample> out{};
  out.reserve(zen_grid.size() * azi_grid.size());
  for (Size izen = 0; izen < zen_grid.size(); ++izen) {
    const Numeric low  = bounds[izen][0];
    const Numeric high = bounds[izen][1];

    // A single zenith offset has no width; its samples then share one ring of arbitrary solid angle
    const Numeric ring = high > low ? Constant::two_pi * (cosd(low) - cosd(high)) : 1.0;

    if (zen_grid[izen] == 0.0) {
      Numeric mean_gain = 0.0;
      for (Size iazi = 0; iazi < azi_grid.size(); ++iazi) mean_gain += gain[izen, iazi] * azi_fraction[iazi];

      out.push_back({.poslos               = {.pos = pos, .los = bore_los},
                     .gain                 = mean_gain,
                     .solid_angle          = ring,
                     .local_los            = {0.0, 0.0},
                     .basis_from_boresight = rtepack::stokes_rotation(1.0, 0.0)});
      continue;
    }

    Vector2 tilted_los         = bore_los;
    tilted_los[0]             += zen_grid[izen];
    const Vector3 tilted_ecef  = geodetic_los2ecef(pos, tilted_los, ell).second;

    for (Size iazi = 0; iazi < azi_grid.size(); ++iazi) {
      const Vector3 sample_ecef = turn(tilted_ecef, los_ecef, azi_grid[iazi]);
      const Vector2 los         = ecef2geodetic_los(pos_ecef, sample_ecef, ell).second;

      // The antenna's field is carried along the great circle from the boresight to the sample; the tilt and turn
      // that place the sample would also spin it about the boresight
      const Vector3 sample = turn(turn(bore, tilt_axis, zen_grid[izen]), bore, azi_grid[iazi]);

      out.push_back({.poslos               = {.pos = pos, .los = los},
                     .gain                 = gain[izen, iazi],
                     .solid_angle          = ring * azi_fraction[iazi],
                     .local_los            = {zen_grid[izen], azi_grid[iazi]},
                     .basis_from_boresight = carried_basis_change(bore, bore_basis, sample, los)});
    }
  }

  return out;
}

//! The channel's integration weights, normalized to sum to one, and their frequencies
struct ChannelQuadrature {
  std::shared_ptr<const AscendingGrid> freq_grid;
  Vector                               weights;
};

[[nodiscard]] ChannelQuadrature normalized_quadrature(const Channel& channel) {
  auto          q     = channel.quadrature();
  const Numeric total = sum(q.data);
  if (total == 0.0)
    arts::user_error("The channel's response integrates to zero, so the channel receives nothing:\n{:B,}", q.data);
  q.data /= total;
  return {.freq_grid = std::make_shared<const AscendingGrid>(std::move(q.grid<0>())), .weights = std::move(q.data)};
}

[[nodiscard]] AntennaPatternField make_gaussian_field(ZenGrid zen_grid, AziGrid azi_grid, Numeric std) {
  if (std <= 0.0) arts::user_error("Gaussian antenna std must be positive");

  AntennaPatternField out{
      .data_name  = "gaussian"s,
      .data       = Matrix(zen_grid.size(), azi_grid.size()),
      .grid_names = {"zenith"s, "azimuth"s},
      .grids      = {std::move(zen_grid), std::move(azi_grid)},
  };

  for (Size izen = 0; izen < out.grid<0>().size(); ++izen) {
    out[izen] = std::exp(-0.5 * Math::pow2(out.grid<0>()[izen] / std));
  }

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

std::shared_ptr<const PosLosVector> make_poslos_grid(const std::vector<AntennaGeometrySample>& samples) {
  return std::make_shared<PosLosVector>(std::from_range,
                                        samples | stdv::transform([](const auto& s) { return s.poslos; }));
}

//! The obsels of an antenna that is a single beam: its one pixel
ArrayOfSensorObsel one_pixel(Obsel&& obsel) {
  ArrayOfSensorObsel out;
  out.push_back(std::move(obsel));
  return out;
}
}  // namespace

Size PixelGrid::size() const { return rows.size() * cols.size(); }

PixelGrid PencilBeamAntenna::pixels() const { return {}; }

PixelGrid GriddedAntennaPattern::pixels() const { return {}; }

PixelGrid GaussianAntenna::pixels() const { return {}; }

PixelGrid GaussianAiryAntenna::pixels() const { return {}; }

ArrayOfSensorObsel PencilBeamAntenna::operator()(const Channel& channel,
                                                 const Stokvec& response,
                                                 const Vector3& pos,
                                                 const Vector2& bore_los,
                                                 const Vector2&) const {
  const auto [freq_grid, channel_weights] = normalized_quadrature(channel);

  SparseStokvecMatrix weight_matrix(1, channel_weights.size());
  for (Size ifreq = 0; ifreq < channel_weights.size(); ++ifreq) {
    if (channel_weights[ifreq] == 0.0) continue;
    weight_matrix[0, ifreq] = channel_weights[ifreq] * response;
  }

  return one_pixel({freq_grid, make_single_poslos_grid(pos, bore_los), std::move(weight_matrix)});
}

ArrayOfSensorObsel GriddedAntennaPattern::operator()(const Channel& channel,
                                                     const Stokvec& response,
                                                     const Vector3& pos,
                                                     const Vector2& bore_los,
                                                     const Vector2& ell) const {
  if (not data.ok()) arts::user_error("GriddedAntennaPattern data shape does not match its grids");

  const auto [freq_grid, channel_weights] = normalized_quadrature(channel);

  const auto antenna_geom = make_antenna_geometry_layout(data.data, data.grid<0>(), data.grid<1>(), pos, bore_los, ell);

  Numeric beam_total = 0.0;
  for (const auto& sample : antenna_geom) beam_total += sample.gain * sample.solid_angle;
  if (beam_total == 0.0) arts::user_error("The antenna pattern has no gain on its grid");

  SparseStokvecMatrix weight_matrix(antenna_geom.size(), channel_weights.size());
  for (Size i = 0; i < antenna_geom.size(); ++i) {
    const Numeric beam = antenna_geom[i].gain * antenna_geom[i].solid_angle / beam_total;
    if (beam == 0.0) continue;

    const Stokvec sample_response = antenna_geom[i].basis_from_boresight * response;
    for (Size j = 0; j < channel_weights.size(); j++) {
      if (channel_weights[j] == 0.0) continue;
      weight_matrix[i, j] = channel_weights[j] * beam * sample_response;
    }
  }

  return one_pixel({freq_grid, make_poslos_grid(antenna_geom), std::move(weight_matrix)});
}

namespace {
AziGrid make_azi_grid(Size size) {
  if (size < 1) arts::user_error("Gaussian antenna azimuth size must be at least 1");

  return Vector{nlinspace(0.0, 360.0, size + 1)[Range(0, size)]};
}
}  // namespace

GaussianAntenna::GaussianAntenna(ZenGrid zen_grid, Numeric std, Size azi_grid_size)
    : zen_grid(std::move(zen_grid)), azi_grid(make_azi_grid(azi_grid_size)), std(std) {}

ArrayOfSensorObsel GaussianAntenna::operator()(const Channel& channel,
                                               const Stokvec& response,
                                               const Vector3& pos,
                                               const Vector2& bore_los,
                                               const Vector2& ell) const {
  const GriddedAntennaPattern pattern{.data = make_gaussian_field(zen_grid, azi_grid, this->std)};
  return pattern(channel, response, pos, bore_los, ell);
}

GaussianAiryAntenna::GaussianAiryAntenna(ZenGrid zen_grid, Numeric aperture_diameter, Size azi_grid_size)
    : zen_grid(std::move(zen_grid)), azi_grid(make_azi_grid(azi_grid_size)), aperture_diameter(aperture_diameter) {}

Numeric GaussianAiryAntenna::std(Numeric frequency) const { return gaussian_airy_std(frequency, aperture_diameter); }

ArrayOfSensorObsel GaussianAiryAntenna::operator()(const Channel& channel,
                                                   const Stokvec& response,
                                                   const Vector3& pos,
                                                   const Vector2& bore_los,
                                                   const Vector2& ell) const {
  if (aperture_diameter <= 0.0) arts::user_error("Gaussian Airy antenna aperture_diameter must be positive");

  const auto [freq_grid, channel_weights] = normalized_quadrature(channel);

  if (freq_grid->front() <= 0.0)
    arts::user_error(
        "Gaussian Airy antenna requires strictly positive channel frequencies because sensor builder only provides the channel frequency grid");

  const auto antenna_geom = make_antenna_geometry_layout(
      Matrix(zen_grid.size(), azi_grid.size(), 1.0), zen_grid, azi_grid, pos, bore_los, ell);

  // The beam at each frequency, normalized over the lines of sight
  Matrix beam(channel_weights.size(), antenna_geom.size(), 0.0);
  for (Size ifreq = 0; ifreq < channel_weights.size(); ++ifreq) {
    if (channel_weights[ifreq] == 0.0) continue;

    const Numeric airy_std = std((*freq_grid)[ifreq]);
    for (Size ipos = 0; ipos < antenna_geom.size(); ++ipos) {
      beam[ifreq, ipos] =
          gaussian_airy_response(antenna_geom[ipos].local_los[0], airy_std) * antenna_geom[ipos].solid_angle;
    }

    beam[ifreq] /= sum(beam[ifreq]);
  }

  SparseStokvecMatrix weight_matrix(antenna_geom.size(), channel_weights.size());
  for (Size i = 0; i < antenna_geom.size(); ++i) {
    const Stokvec sample_response = antenna_geom[i].basis_from_boresight * response;
    for (Size j = 0; j < channel_weights.size(); j++) {
      if (channel_weights[j] == 0.0 or beam[j, i] == 0.0) continue;
      weight_matrix[i, j] = channel_weights[j] * beam[j, i] * sample_response;
    }
  }

  return one_pixel({freq_grid, make_poslos_grid(antenna_geom), std::move(weight_matrix)});
}

CameraAntenna CameraAntenna::thin_lens(Size        n_rows,
                                       Size        n_cols,
                                       Numeric     ccd_height,
                                       Numeric     ccd_width,
                                       Numeric     focal_length,
                                       Numeric     focus_distance,
                                       BeamPattern beam) {
  if (n_rows < 1 or n_cols < 1)
    arts::user_error("A camera needs at least one row and one column of pixels, got {} by {}", n_rows, n_cols);
  if (ccd_height <= 0.0 or ccd_width <= 0.0)
    arts::user_error("A camera needs a CCD of positive size, got {} by {} m", ccd_height, ccd_width);
  if (focal_length <= 0.0) arts::user_error("A camera needs a positive focal length, got {} m", focal_length);
  if (focus_distance <= focal_length)
    arts::user_error("A camera must be focused beyond its focal length, got focus distance {} m and focal length {} m",
                     focus_distance,
                     focal_length);

  // The thin-lens equation 1 / f = 1 / d_o + 1 / d_i gives the distance from the lens to the CCD
  const Numeric image_distance = focal_length * focus_distance / (focus_distance - focal_length);

  // The angle of each pixel centre seen from the lens, along one side of the CCD
  const auto offsets = [image_distance](Size n, Numeric size) {
    const Numeric pitch = size / static_cast<Numeric>(n);
    Vector        out(n);
    for (Size i = 0; i < n; ++i) {
      const Numeric centre = (static_cast<Numeric>(i) - 0.5 * static_cast<Numeric>(n - 1)) * pitch;
      out[i]               = Conversion::atan2d(centre, image_distance);
    }
    return AscendingGrid{std::move(out)};
  };

  return {.grid = {.rows = offsets(n_rows, ccd_height), .cols = offsets(n_cols, ccd_width)}, .beam = std::move(beam)};
}

PixelGrid CameraAntenna::pixels() const { return grid; }

ArrayOfSensorObsel CameraAntenna::operator()(const Channel& channel,
                                             const Stokvec& response,
                                             const Vector3& pos,
                                             const Vector2& bore_los,
                                             const Vector2& ell) const {
  using Conversion::cosd, Conversion::sind, Conversion::tand;

  if (grid.rows.empty() or grid.cols.empty())
    arts::user_error("A camera needs at least one row and one column of pixels, got {} by {}",
                     grid.rows.size(),
                     grid.cols.size());
  for (const AscendingGrid* offsets : {&grid.rows, &grid.cols}) {
    if (offsets->front() <= -90.0 or offsets->back() >= 90.0)
      arts::user_error("A camera's pixel offsets must be within 90 degrees of the optical axis, got {:B,}", *offsets);
  }

  // In the local east-north-up frame: the optical axis, the directions of increasing zenith and azimuth there,
  // and the polarization basis of the radiation it receives
  const Numeric za         = bore_los[0];
  const Numeric aa         = bore_los[1];
  const Vector3 bore       = enu_direction(bore_los);
  const Vector3 e_zen      = {cosd(za) * sind(aa), cosd(za) * cosd(aa), -sind(za)};
  const Vector3 e_azi      = {cosd(aa), -sind(aa), 0.0};
  const auto    bore_basis = received_basis(bore_los);

  ArrayOfSensorObsel out;
  out.reserve(grid.size());
  for (const Numeric row : grid.rows) {
    for (const Numeric col : grid.cols) {
      // CCD coordinates: the lens inverts the image, so a positive column offset looks towards decreasing azimuth
      Vector3 direction = bore + e_zen * tand(row) - e_azi * tand(col);
      direction         = direction / hypot(direction);

      const Vector2 los{Conversion::acosd(std::clamp(direction[2], -1.0, 1.0)),
                        Conversion::atan2d(direction[0], direction[1])};

      const Stokvec pixel_response = carried_basis_change(bore, bore_basis, direction, los) * response;

      auto obsels = std::visit([&](const auto& b) { return b(channel, pixel_response, pos, los, ell); }, beam);
      out.push_back(std::move(obsels.front()));
    }
  }

  return out;
}

PixelGrid antenna_pixels(const AntennaPattern& antenna) {
  return std::visit([](const auto& pattern) { return pattern.pixels(); }, antenna);
}

ArrayOfSensorObsel antenna_obsels(const AntennaPattern& antenna,
                                  const Channel&        channel,
                                  const Stokvec&        response,
                                  const Vector3&        pos,
                                  const Vector2&        bore_los,
                                  const Vector2&        ell) {
  return std::visit([&](const auto& pattern) { return pattern(channel, response, pos, bore_los, ell); }, antenna);
}

static_assert(AntennaPatternSelection<GriddedAntennaPattern>);
static_assert(AntennaPatternSelection<PencilBeamAntenna>);
static_assert(AntennaPatternSelection<GaussianAntenna>);
static_assert(AntennaPatternSelection<GaussianAiryAntenna>);
static_assert(AntennaPatternSelection<CameraAntenna>);
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

void xml_io_stream<sensor::GaussianAntenna>::write(std::ostream&                  os,
                                                   const sensor::GaussianAntenna& n,
                                                   bofstream*                     pbofs,
                                                   std::string_view               name) {
  XMLTag tag(xml_io_stream_name_v<sensor::GaussianAntenna>, "name", name);
  tag.write_to_stream(os);
  xml_write_to_stream(os, n.std, pbofs);
  xml_write_to_stream(os, n.zen_grid, pbofs);
  xml_write_to_stream(os, n.azi_grid, pbofs);
  tag.write_to_end_stream(os);
}

void xml_io_stream<sensor::GaussianAntenna>::read(std::istream& is, sensor::GaussianAntenna& n, bifstream* pbifs) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(xml_io_stream_name_v<sensor::GaussianAntenna>);

  xml_read_from_stream(is, n.std, pbifs);
  xml_read_from_stream(is, n.zen_grid, pbifs);
  xml_read_from_stream(is, n.azi_grid, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(xml_io_stream_name_v<sensor::GaussianAntenna>);
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

  tag.read_from_stream(is);
  tag.check_end_name(xml_io_stream_name_v<sensor::GaussianAiryAntenna>);
}

void xml_io_stream<sensor::PencilBeamAntenna>::write(std::ostream& os,
                                                     const sensor::PencilBeamAntenna&,
                                                     bofstream*,
                                                     std::string_view name) {
  XMLTag tag(xml_io_stream_name_v<sensor::PencilBeamAntenna>, "name", name);
  tag.write_to_stream(os);
  tag.write_to_end_stream(os);
}

void xml_io_stream<sensor::PencilBeamAntenna>::read(std::istream& is, sensor::PencilBeamAntenna&, bifstream*) {
  XMLTag tag{};
  tag.read_from_stream(is);
  tag.check_name(xml_io_stream_name_v<sensor::PencilBeamAntenna>);
  tag.read_from_stream(is);
  tag.check_end_name(xml_io_stream_name_v<sensor::PencilBeamAntenna>);
}
