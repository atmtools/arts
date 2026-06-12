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
struct Builder;
struct FrequencyRange;

//! Concept for selecting sensor builders.
template <typename T>
concept BuilderSelection = std::derived_from<T, Builder>;

struct Builder {
  std::vector<Channel> channels;
  std::shared_ptr<const AntennaPattern> antenna;
  bool preserve_common_frequency_grid{false};

  Builder();
  Builder(std::vector<Channel> channels,
          std::shared_ptr<const AntennaPattern> antenna);
  Builder(const Spectrometer& spectrometer,
          std::shared_ptr<const AntennaPattern> antenna);
  Builder(const Spectrometer& spectrometer,
          const FrequencyRange& backend,
          std::shared_ptr<const AntennaPattern> antenna);
  Builder(const Builder& other) = default;
  Builder(Builder&&) noexcept   = default;
  Builder& operator=(const Builder& other);
  Builder& operator=(Builder&&) noexcept = default;

  [[nodiscard]] std::pair<ArrayOfSensorObsel, ArrayOfSensorMetaInfo> operator()(
      std::span<const Vector3> pos, std::span<const Vector2> los) const;
};
}  // namespace sensor
