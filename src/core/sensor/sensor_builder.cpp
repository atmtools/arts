#include "sensor_builder.h"

#include <debug.h>
#include <frequency_bandpass_filters.h>

#include <cmath>
#include <format>
#include <memory>
#include <optional>

namespace sensor {
namespace {
SensorMetaInfo make_meta_info(Size nchannels, Size geometry_index, std::optional<Time> time = std::nullopt) {
  SortedGriddedField1 gf;
  gf.data_name     = std::format("sensor-builder-{}", geometry_index);
  gf.gridname<0>() = "channel";
  Vector channel_axis(nchannels);
  for (Size i = 0; i < nchannels; ++i) { channel_axis[i] = static_cast<Numeric>(i); }
  gf.grid<0>() = AscendingGrid{std::move(channel_axis)};
  gf.data.resize(nchannels);
  gf.data = 0.0;

  return SensorMetaInfo{.data = std::move(gf), .time = std::move(time)};
}
}  // namespace

Builder Builder::from(const Spectrometer& spectrometer, AntennaPattern antenna) {
  return {.channels = spectrometer.channels, .antenna = std::move(antenna), .preserve_common_frequency_grid = true};
}

Builder Builder::from(const Spectrometer& spectrometer, const FrequencyRange& backend, AntennaPattern antenna) {
  return {.channels                       = make_bandpass_channels(backend, spectrometer),
          .antenna                        = std::move(antenna),
          .preserve_common_frequency_grid = true};
}

namespace {
struct BoreGeometry {
  Vector3 pos;
  Vector2 los;
  //! Signed scan angle [deg], when the geometry is a scan
  std::optional<Numeric> scan_angle;
  std::optional<Time>    time;
};

std::pair<ArrayOfSensorObsel, ArrayOfSensorMetaInfo> build(const Builder&                   builder,
                                                           const std::vector<BoreGeometry>& geometries,
                                                           const Vector2&                   ell) {
  const auto& channels = builder.channels;

  if (channels.empty()) arts::user_error("Builder requires at least one channel");
  if (geometries.empty()) arts::user_error("Builder requires at least one geometry");

  ArrayOfSensorObsel out;
  out.reserve(geometries.size() * channels.size());

  ArrayOfSensorMetaInfo meta;
  meta.reserve(geometries.size());

  std::vector<std::shared_ptr<const AscendingGrid>> freq_grids(channels.size());

  for (Size igeom = 0; igeom < geometries.size(); ++igeom) {
    const auto&                         geometry = geometries[igeom];
    std::shared_ptr<const PosLosVector> poslos_grid;

    for (Size ichan = 0; ichan < channels.size(); ++ichan) {
      const auto& polarization = channels[ichan].polarization;
      if (polarization.scan_rotation != 0.0 and not geometry.scan_angle)
        arts::user_error(
            R"(Channel {} has a polarization that turns with the scan, scan_rotation = {},
but the side of the scan of a bare line of sight is not known.  Build the sensor from a scan geometry instead.
)",
            ichan,
            polarization.scan_rotation);

      const Stokvec response = polarization.at(geometry.scan_angle.value_or(0.0));

      auto obsel = antenna_obsel(builder.antenna, channels[ichan], response, geometry.pos, geometry.los, ell);

      // The same channel has the same frequencies at every geometry
      if (not freq_grids[ichan]) {
        freq_grids[ichan] = obsel.f_grid_ptr();
      } else {
        obsel.set_f_grid_ptr(freq_grids[ichan]);
      }

      // The antenna samples the same lines of sight for every channel
      if (not poslos_grid) {
        poslos_grid = obsel.poslos_grid_ptr();
      } else {
        obsel.set_poslos_grid_ptr(poslos_grid);
      }

      out.emplace_back(std::move(obsel));
    }

    meta.push_back(make_meta_info(channels.size(), igeom, geometry.time));
  }

  if (builder.preserve_common_frequency_grid) collect_frequency_grids(out);

  return {std::move(out), std::move(meta)};
}
}  // namespace

std::pair<ArrayOfSensorObsel, ArrayOfSensorMetaInfo> Builder::operator()(std::span<const Vector3> pos,
                                                                         std::span<const Vector2> los,
                                                                         const Vector2&           ell) const {
  return operator()(pos, los, {}, ell);
}

std::pair<ArrayOfSensorObsel, ArrayOfSensorMetaInfo> Builder::operator()(std::span<const Vector3> pos,
                                                                         std::span<const Vector2> los,
                                                                         std::span<const Time>    time,
                                                                         const Vector2&           ell) const {
  if (channels.empty()) arts::user_error("Builder requires at least one channel");
  if (pos.empty()) arts::user_error("Builder requires at least one sensor position");
  if (los.empty()) arts::user_error("Builder requires at least one bore LOS");
  if (pos.size() != los.size())
    arts::user_error("Builder requires matching position and LOS counts. Got {} positions and {} LOS values.",
                     pos.size(),
                     los.size());
  if (not time.empty() and time.size() != pos.size())
    arts::user_error("Builder requires either no times or one time per geometry. Got {} times and {} geometries.",
                     time.size(),
                     pos.size());

  std::vector<BoreGeometry> geometries;
  geometries.reserve(pos.size());
  for (Size i = 0; i < pos.size(); ++i) {
    geometries.push_back({.pos        = pos[i],
                          .los        = los[i],
                          .scan_angle = std::nullopt,
                          .time       = time.empty() ? std::nullopt : std::optional<Time>{time[i]}});
  }

  return build(*this, geometries, ell);
}

std::pair<ArrayOfSensorObsel, ArrayOfSensorMetaInfo> Builder::cross_track(std::span<const Vector3> pos,
                                                                          std::span<const Numeric> heading,
                                                                          const Vector&            scan_angles,
                                                                          std::span<const Time>    time,
                                                                          const Vector2&           ell) const {
  if (pos.empty()) arts::user_error("Builder requires at least one platform position");
  if (scan_angles.empty()) arts::user_error("Builder requires at least one scan angle");
  if (pos.size() != heading.size())
    arts::user_error(
        "Builder requires one heading per position. Got {} positions and {} headings.", pos.size(), heading.size());
  if (not time.empty() and time.size() != pos.size())
    arts::user_error("Builder requires either no times or one time per position. Got {} times and {} positions.",
                     time.size(),
                     pos.size());

  std::vector<BoreGeometry> geometries;
  geometries.reserve(pos.size() * scan_angles.size());
  for (Size i = 0; i < pos.size(); ++i) {
    for (const Numeric s : scan_angles) {
      const Numeric azimuth = std::remainder(heading[i] + (s < 0.0 ? -90.0 : 90.0), 360.0);
      geometries.push_back({.pos        = pos[i],
                            .los        = {180.0 - std::abs(s), azimuth},
                            .scan_angle = s,
                            .time       = time.empty() ? std::nullopt : std::optional<Time>{time[i]}});
    }
  }

  return build(*this, geometries, ell);
}
}  // namespace sensor
