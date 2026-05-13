#pragma once

#include <concepts>
#include <memory>
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
  std::shared_ptr<const AntennaPattern> antenna;

  SensorBuilder();
  SensorBuilder(std::vector<Channel> channels, const AntennaPattern& antenna);
  SensorBuilder(const SensorBuilder& other);
  SensorBuilder(SensorBuilder&&) noexcept = default;
  SensorBuilder& operator=(const SensorBuilder& other);
  SensorBuilder& operator=(SensorBuilder&&) noexcept = default;

  [[nodiscard]] const AntennaPattern& get_antenna() const;
  void set_antenna(const AntennaPattern& pattern);

  [[nodiscard]] std::pair<ArrayOfSensorObsel, ArrayOfSensorMetaInfo> operator()(
      std::span<const Vector3> pos, std::span<const Vector2> los) const;
};
}  // namespace sensor
