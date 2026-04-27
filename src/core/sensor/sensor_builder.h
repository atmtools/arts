#pragma once

#include <concepts>
#include <span>
#include <utility>
#include <vector>

#include "antenna_pattern.h"
#include "frequency_channel_selection.h"
#include "obsel.h"
#include "sensor_meta_info.h"

namespace sensor {
//! Combines channels with an antenna pattern and bore geometries into obsels.
struct SensorBuilder;

//! Concept for selecting sensor builders.
template <typename T>
concept SensorBuilderSelection = std::derived_from<T, SensorBuilder>;

struct SensorBuilder {
  std::vector<Channel> channels;
  AntennaPattern antenna;

  [[nodiscard]] std::pair<ArrayOfSensorObsel, ArrayOfSensorMetaInfo> operator()(
      std::span<const Vector3> pos, std::span<const Vector2> los) const;
};
}  // namespace sensor
