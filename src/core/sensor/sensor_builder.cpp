#include "sensor_builder.h"

#include <debug.h>
#include <frequency_bandpass_filters.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <optional>

namespace sensor {
namespace {
//! The indices 0, 1, ..., n - 1
AscendingGrid index_axis(Size n) {
  Vector out(n);
  for (Size i = 0; i < n; ++i) out[i] = static_cast<Numeric>(i);
  return AscendingGrid{std::move(out)};
}

//! The structure of the elements: placement by channel, by pixel rows and columns when the antenna has more than one
SensorMetaInfo make_meta_info(const String& name, const PixelGrid& pixels, Size nplacements, Size nchannels) {
  if (pixels.size() == 1) {
    SortedGriddedField2 gf;
    gf.data_name     = name;
    gf.gridname<0>() = "placement";
    gf.gridname<1>() = "channel";
    gf.grid<0>()     = index_axis(nplacements);
    gf.grid<1>()     = index_axis(nchannels);
    gf.data.resize(nplacements, nchannels);
    gf.data = 0.0;

    return gf;
  }

  SortedGriddedField4 gf;
  gf.data_name     = name;
  gf.gridname<0>() = "placement";
  gf.gridname<1>() = "row_offset_deg";
  gf.gridname<2>() = "col_offset_deg";
  gf.gridname<3>() = "channel";
  gf.grid<0>()     = index_axis(nplacements);
  gf.grid<1>()     = pixels.rows;
  gf.grid<2>()     = pixels.cols;
  gf.grid<3>()     = index_axis(nchannels);
  gf.data.resize(nplacements, pixels.rows.size(), pixels.cols.size(), nchannels);
  gf.data = 0.0;

  return gf;
}
}  // namespace

Builder Builder::from(const Spectrometer& spectrometer, AntennaPattern antenna) {
  return {.channels = spectrometer.channels, .antenna = std::move(antenna)};
}

Builder Builder::from(const Spectrometer& spectrometer, const FrequencyRange& backend, AntennaPattern antenna) {
  return {.channels = make_bandpass_channels(backend, spectrometer), .antenna = std::move(antenna)};
}

namespace {
struct BoreGeometry {
  Vector3 pos;
  Vector2 los;
  //! Signed scan angle [deg], when the geometry is a scan
  std::optional<Numeric> scan_angle;
};

//! The name of the sensor: the given name, or one after its number of channels if the given name is empty
String sensor_name(const Builder& builder, const String& name) {
  return name.empty() ? std::format("sensor-with-{}-channels", builder.channels.size()) : name;
}

Sensor build(const Builder&                   builder,
             const std::vector<BoreGeometry>& placements,
             const Vector2&                   ell,
             const String&                    name) {
  const auto& channels = builder.channels;

  if (channels.empty()) arts::user_error(R"(Sensor "{}": the builder requires at least one channel)", name);
  if (placements.empty()) arts::user_error(R"(Sensor "{}": the builder requires at least one geometry)", name);

  const PixelGrid pixels  = antenna_pixels(builder.antenna);
  const Size      npixels = pixels.size();

  // The sensor of each channel at each placement, one element per pixel
  std::vector<std::vector<Sensor>> placed(placements.size());
  for (Size iplace = 0; iplace < placements.size(); ++iplace) {
    const auto& placement = placements[iplace];

    if (npixels > 1 and placement.scan_angle)
      arts::user_error(
          R"(Sensor "{}": a scan of an antenna with more than one pixel is not defined, got {} rows by {} columns)",
          name,
          pixels.rows.size(),
          pixels.cols.size());

    placed[iplace].reserve(channels.size());
    for (Size ichan = 0; ichan < channels.size(); ++ichan) {
      const auto& polarization = channels[ichan].polarization;
      if (polarization.scan_rotation != 0.0 and not placement.scan_angle)
        arts::user_error(
            R"(Sensor "{}": channel {} has a polarization that turns with the scan, scan_rotation = {},
but the side of the scan of a bare line of sight is not known.  Build the sensor from a scan geometry instead.
)",
            name,
            ichan,
            polarization.scan_rotation);

      const Stokvec response = polarization.at(placement.scan_angle.value_or(0.0));

      auto& sensor = placed[iplace].emplace_back([&] {
        try {
          return antenna_sensor(builder.antenna, channels[ichan], response, placement.pos, placement.los, ell);
        } catch (const std::exception& e) {
          throw std::runtime_error(std::format(R"(Sensor "{}", channel {}, placement {}:
{})",
                                               name,
                                               ichan,
                                               iplace,
                                               e.what()));
        }
      }());
      if (sensor.elements.size() != npixels)
        arts::user_error(
            R"(Sensor "{}": the antenna gave {} elements for its {} pixels)", name, sensor.elements.size(), npixels);
    }
  }

  // The frequencies of all channels, and where each channel's frequencies are among them
  std::vector<Numeric> all_freqs;
  for (const auto& sensor : placed.front()) all_freqs.append_range(sensor.freq_grid);
  stdr::sort(all_freqs);
  all_freqs.erase(std::unique(all_freqs.begin(), all_freqs.end()), all_freqs.end());
  AscendingGrid freq_grid{std::move(all_freqs)};

  std::vector<std::vector<Size>> freq_maps(channels.size());
  for (Size ichan = 0; ichan < channels.size(); ++ichan) {
    for (const Numeric f : placed.front()[ichan].freq_grid) {
      freq_maps[ichan].push_back(static_cast<Size>(stdr::distance(freq_grid.begin(), stdr::lower_bound(freq_grid, f))));
    }
  }

  // The lines of sight of each placement, which the antenna samples the same for every channel
  std::vector<Size> geometry_offsets(placements.size());
  Size              ngeom = 0;
  for (Size iplace = 0; iplace < placements.size(); ++iplace) {
    const auto& geometry = placed[iplace].front().geometry;
    for (Size ichan = 1; ichan < channels.size(); ++ichan) {
      if (not stdr::equal(placed[iplace][ichan].geometry, geometry))
        arts::user_error(
            R"(Sensor "{}": the antenna samples other lines of sight for channel {} than for channel 0 at placement {})",
            name,
            ichan,
            iplace);
    }

    geometry_offsets[iplace]  = ngeom;
    ngeom                    += geometry.size();
  }

  Sensor out{.freq_grid = std::move(freq_grid),
             .geometry  = PosLosVector(ngeom),
             .elements  = {},
             .meta      = make_meta_info(name, pixels, placements.size(), channels.size())};

  // By placement, then by pixel and then by channel, as the meta
  out.elements.reserve(placements.size() * npixels * channels.size());
  for (Size iplace = 0; iplace < placements.size(); ++iplace) {
    const auto& geometry = placed[iplace].front().geometry;
    for (Size i = 0; i < geometry.size(); ++i) out.geometry[geometry_offsets[iplace] + i] = geometry[i];

    for (Size ipix = 0; ipix < npixels; ++ipix) {
      for (Size ichan = 0; ichan < channels.size(); ++ichan) {
        out.elements.push_back(embed(placed[iplace][ichan].elements[ipix],
                                     ngeom,
                                     out.freq_grid.size(),
                                     geometry_offsets[iplace],
                                     freq_maps[ichan]));
      }
    }
  }

  out.check();
  return out;
}
}  // namespace

Sensor Builder::operator()(std::span<const Vector3> pos,
                           std::span<const Vector2> los,
                           const Vector2&           ell,
                           const String&            name) const {
  const String sensor = sensor_name(*this, name);

  if (pos.empty()) arts::user_error(R"(Sensor "{}": the builder requires at least one sensor position)", sensor);
  if (los.empty()) arts::user_error(R"(Sensor "{}": the builder requires at least one bore LOS)", sensor);
  if (pos.size() != los.size())
    arts::user_error(
        R"(Sensor "{}": the builder requires matching position and LOS counts. Got {} positions and {} LOS values.)",
        sensor,
        pos.size(),
        los.size());

  std::vector<BoreGeometry> placements;
  placements.reserve(pos.size());
  for (Size i = 0; i < pos.size(); ++i) {
    placements.push_back({.pos = pos[i], .los = los[i], .scan_angle = std::nullopt});
  }

  return build(*this, placements, ell, sensor);
}

Sensor Builder::cross_track(std::span<const Vector3> pos,
                            std::span<const Numeric> heading,
                            const Vector&            scan_angles,
                            const Vector2&           ell,
                            const String&            name) const {
  const String sensor = sensor_name(*this, name);

  if (pos.empty()) arts::user_error(R"(Sensor "{}": the builder requires at least one platform position)", sensor);
  if (scan_angles.empty()) arts::user_error(R"(Sensor "{}": the builder requires at least one scan angle)", sensor);
  if (pos.size() != heading.size())
    arts::user_error(R"(Sensor "{}": the builder requires one heading per position. Got {} positions and {} headings.)",
                     sensor,
                     pos.size(),
                     heading.size());

  std::vector<BoreGeometry> placements;
  placements.reserve(pos.size() * scan_angles.size());
  for (Size i = 0; i < pos.size(); ++i) {
    for (const Numeric s : scan_angles) {
      const Numeric azimuth = std::remainder(heading[i] + (s < 0.0 ? -90.0 : 90.0), 360.0);
      placements.push_back({.pos = pos[i], .los = {180.0 - std::abs(s), azimuth}, .scan_angle = s});
    }
  }

  return build(*this, placements, ell, sensor);
}
}  // namespace sensor
