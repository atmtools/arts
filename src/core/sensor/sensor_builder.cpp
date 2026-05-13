#include "sensor_builder.h"

#include <debug.h>

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

SensorBuilder::SensorBuilder() : antenna(PencilBeamAntenna{}.clone()) {}

SensorBuilder::SensorBuilder(std::vector<Channel> channels,
                             const AntennaPattern& antenna)
    : channels(std::move(channels)), antenna(antenna.clone()) {}

SensorBuilder::SensorBuilder(const SensorBuilder& other)
    : channels(other.channels), antenna(clone_antenna(other.antenna)) {}

SensorBuilder& SensorBuilder::operator=(const SensorBuilder& other) {
  if (this != &other) {
    channels = other.channels;
    antenna  = clone_antenna(other.antenna);
  }

  return *this;
}

const AntennaPattern& SensorBuilder::get_antenna() const {
  ARTS_USER_ERROR_IF(not antenna, "SensorBuilder requires an antenna pattern")
  return *antenna;
}

void SensorBuilder::set_antenna(const AntennaPattern& pattern) {
  antenna = pattern.clone();
}

std::pair<ArrayOfSensorObsel, ArrayOfSensorMetaInfo> SensorBuilder::operator()(
    std::span<const Vector3> pos, std::span<const Vector2> los) const {
  ARTS_USER_ERROR_IF(channels.empty(),
                     "SensorBuilder requires at least one channel")
  ARTS_USER_ERROR_IF(not antenna, "SensorBuilder requires an antenna pattern")
  ARTS_USER_ERROR_IF(pos.empty(),
                     "SensorBuilder requires at least one sensor position")
  ARTS_USER_ERROR_IF(los.empty(),
                     "SensorBuilder requires at least one bore LOS")
  ARTS_USER_ERROR_IF(
      pos.size() != los.size(),
      "SensorBuilder requires matching position and LOS counts. Got {} positions and {} LOS values.",
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
                                   Size geometry_index) {
    std::shared_ptr<const PosLosVector> poslos_grid;

    for (Size ichan = 0; ichan < channels.size(); ++ichan) {
      auto obsel = get_antenna()(channels[ichan], sensor_pos, bore_los);
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

  for (Size i = 0; i < pos.size(); ++i) append_geometry(pos[i], los[i], i);

  return {std::move(out), std::move(meta)};
}

static_assert(SensorBuilderSelection<SensorBuilder>);
}  // namespace sensor
