#include "sensor_builder.h"

#include <debug.h>
#include <frequency_bandpass_filters.h>

#include <format>
#include <memory>

namespace sensor {
namespace {
std::shared_ptr<const AntennaPattern> clone_antenna(
    const std::shared_ptr<const AntennaPattern>& antenna) {
  return antenna ? antenna->clone() : nullptr;
}

SensorMetaInfo make_meta_info(Size nchannels, Size geometry_index) {
  SortedGriddedField1 gf;
  gf.data_name     = std::format("sensor-builder-{}", geometry_index);
  gf.gridname<0>() = "channel";
  Vector channel_axis(nchannels);
  for (Size i = 0; i < nchannels; ++i) {
    channel_axis[i] = static_cast<Numeric>(i);
  }
  gf.grid<0>() = AscendingGrid{std::move(channel_axis)};
  gf.data.resize(nchannels);
  gf.data = 0.0;

  return SensorMetaInfo{std::move(gf)};
}
}  // namespace

Builder::Builder() : antenna(PencilBeamAntenna{}.clone()) {}

Builder::Builder(std::vector<Channel> channels,
                 std::shared_ptr<const AntennaPattern> antenna)
    : channels(std::move(channels)), antenna(std::move(antenna)) {}

Builder::Builder(const Spectrometer& spectrometer,
                 std::shared_ptr<const AntennaPattern> antenna)
    : Builder(std::vector<Channel>{spectrometer.channels}, std::move(antenna)) {
  preserve_common_frequency_grid = spectrometer.is_synced();
}

Builder::Builder(const Spectrometer& spectrometer,
                 const FrequencyRange& backend,
                 std::shared_ptr<const AntennaPattern> antenna)
    : Builder(make_bandpass_channels(backend, spectrometer),
              std::move(antenna)) {
  preserve_common_frequency_grid = spectrometer.is_synced();
}

Builder& Builder::operator=(const Builder& other) {
  if (this != &other) {
    channels                       = other.channels;
    antenna                        = clone_antenna(other.antenna);
    preserve_common_frequency_grid = other.preserve_common_frequency_grid;
  }

  return *this;
}

std::pair<ArrayOfSensorObsel, ArrayOfSensorMetaInfo> Builder::operator()(
    std::span<const Vector3> pos,
    std::span<const Vector2> los,
    const Vector2& ell) const {
  ARTS_USER_ERROR_IF(channels.empty(), "Builder requires at least one channel")
  ARTS_USER_ERROR_IF(not antenna, "Builder requires an antenna pattern")
  ARTS_USER_ERROR_IF(pos.empty(),
                     "Builder requires at least one sensor position")
  ARTS_USER_ERROR_IF(los.empty(), "Builder requires at least one bore LOS")
  ARTS_USER_ERROR_IF(
      pos.size() != los.size(),
      "Builder requires matching position and LOS counts. Got {} positions and {} LOS values.",
      pos.size(),
      los.size())

  std::vector<std::shared_ptr<const AscendingGrid>> freq_grids;
  freq_grids.reserve(channels.size());
  for (const auto& channel : channels) {
    freq_grids.push_back(
        std::make_shared<const AscendingGrid>(channel.freq_grid()));
  }

  ArrayOfSensorObsel out;
  out.reserve(pos.size() * channels.size());

  ArrayOfSensorMetaInfo meta;
  meta.reserve(pos.size());

  const auto append_geometry = [&](const Vector3& sensor_pos,
                                   const Vector2& bore_los,
                                   const Vector2& ell,
                                   Size geometry_index) {
    std::shared_ptr<const PosLosVector> poslos_grid;

    for (Size ichan = 0; ichan < channels.size(); ++ichan) {
      auto obsel =
          antenna->operator()(channels[ichan], sensor_pos, bore_los, ell);
      obsel.set_f_grid_ptr(freq_grids[ichan]);

      if (not poslos_grid) {
        poslos_grid = obsel.poslos_grid_ptr();
      } else {
        obsel.set_poslos_grid_ptr(poslos_grid);
      }

      out.emplace_back(std::move(obsel));
    }

    meta.push_back(make_meta_info(channels.size(), geometry_index));
  };

  for (Size i = 0; i < pos.size(); ++i) append_geometry(pos[i], los[i], ell, i);

  if (preserve_common_frequency_grid) collect_frequency_grids(out);

  return {std::move(out), std::move(meta)};
}

static_assert(BuilderSelection<Builder>);
}  // namespace sensor
