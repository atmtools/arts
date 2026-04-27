#include "sensor_builder.h"

#include <debug.h>

#include <format>
#include <memory>

namespace sensor {
namespace {
using AntennaSamples = std::vector<std::pair<Stokvec, Vector2>>;

std::shared_ptr<const PosLosVector> make_poslos_grid(
    const Vector3& pos, const AntennaSamples& antenna_samples) {
  auto out = std::make_shared<PosLosVector>(antenna_samples.size());

  for (Size i = 0; i < antenna_samples.size(); ++i) {
    (*out)[i] = {.pos = pos, .los = antenna_samples[i].second};
  }

  return out;
}

SparseStokvecMatrix make_weight_matrix(const AntennaSamples& antenna_samples,
                                       const Channel& channel) {
  const auto& channel_weights = channel.weights();

  SparseStokvecMatrix out(antenna_samples.size(), channel_weights.size());

  for (Size iposlos = 0; iposlos < antenna_samples.size(); ++iposlos) {
    const auto& [antenna_weight, ignored_los] = antenna_samples[iposlos];
    static_cast<void>(ignored_los);

    if (antenna_weight.is_zero()) continue;

    for (Size ifreq = 0; ifreq < channel_weights.size(); ++ifreq) {
      if (channel_weights[ifreq] == 0.0) continue;
      out[iposlos, ifreq] = channel_weights[ifreq] * antenna_weight;
    }
  }

  return out;
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

std::pair<ArrayOfSensorObsel, ArrayOfSensorMetaInfo> SensorBuilder::operator()(
    std::span<const Vector3> pos, std::span<const Vector2> los) const {
  ARTS_USER_ERROR_IF(channels.empty(),
                     "SensorBuilder requires at least one channel")
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

  std::vector<AntennaSamples> antenna_samples;
  antenna_samples.reserve(los.size());
  for (const auto& bore_los : los) antenna_samples.push_back(antenna(bore_los));

  std::vector<std::vector<SparseStokvecMatrix>> weight_cache(los.size());
  for (Size ilos = 0; ilos < los.size(); ++ilos) {
    auto& cached = weight_cache[ilos];
    cached.reserve(channels.size());
    for (const auto& channel : channels) {
      cached.push_back(make_weight_matrix(antenna_samples[ilos], channel));
    }
  }

  ArrayOfSensorObsel out;
  out.reserve(pos.size() * channels.size());

  ArrayOfSensorMetaInfo meta;
  meta.reserve(pos.size());

  const auto append_geometry =
      [&](const Vector3& sensor_pos, Size ilos, Size geometry_index) {
        auto poslos_grid = make_poslos_grid(sensor_pos, antenna_samples[ilos]);

        for (Size ichan = 0; ichan < channels.size(); ++ichan) {
          out.emplace_back(
              freq_grids[ichan], poslos_grid, weight_cache[ilos][ichan]);
        }

        meta.push_back(make_meta_info(channels.size(), geometry_index));
      };

  for (Size i = 0; i < pos.size(); ++i) append_geometry(pos[i], i, i);

  return {std::move(out), std::move(meta)};
}

static_assert(SensorBuilderSelection<SensorBuilder>);
}  // namespace sensor
