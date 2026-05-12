#pragma once

#include <matpack.h>
#include <rtepack.h>

#include <concepts>

#include "frequency_channel_selection.h"
#include "obsel.h"

namespace sensor {
//! A 2D angular antenna pattern on local zenith and azimuth offsets.
struct AntennaPattern;

//! A 1x1 antenna pattern that samples only the bore line of sight.
struct PencilBeamAntenna;

//! A 2D Gaussian antenna pattern on local zenith and azimuth offsets.
struct GaussianAntenna;

//! 2D gridded field of antenna weights on local zenith and azimuth offsets.
using AntennaPatternField = matpack::gridded_data_t<Stokvec, ZenGrid, AziGrid>;

struct AntennaPattern {
  AntennaPatternField data;  // center at [0, 0]

  //! Builds one sensor obsel for a channel at a sensor position and bore LOS.
  [[nodiscard]] Obsel operator()(const Channel& channel,
                                 const Vector3& pos,
                                 const Vector2& bore_los) const;
};

struct PencilBeamAntenna final : AntennaPattern {
  PencilBeamAntenna(Stokvec weight = {1.0, 0.0, 0.0, 0.0});
};

struct GaussianAntenna final : AntennaPattern {
  GaussianAntenna(ZenGrid zen_grid,
                  AziGrid azi_grid,
                  Numeric zenith_std,
                  Numeric azimuth_std,
                  Stokvec weight = {1.0, 0.0, 0.0, 0.0});
};

//! Concept for types that can build a sensor obsel from channel geometry.
template <typename T>
concept AntennaPatternSelection = requires(const T& antenna,
                                           const Channel& channel,
                                           const Vector3& pos,
                                           const Vector2& bore_los) {
  { antenna(channel, pos, bore_los) } -> std::same_as<Obsel>;
};
}  // namespace sensor

// AntennaPattern format tags and XML I/O

template <>
struct format_tag_aggregate<sensor::AntennaPattern> {
  constexpr static bool value = true;
};

template <>
struct xml_io_stream_name<sensor::AntennaPattern> {
  static constexpr std::string_view name = "SensorAntennaPattern";
};

template <>
struct xml_io_stream_aggregate<sensor::AntennaPattern> {
  static constexpr bool value = true;
};

// PencilBeamAntenna format tags and XML I/O

template <>
struct format_tag_aggregate<sensor::PencilBeamAntenna> {
  constexpr static bool value = true;
};

template <>
struct xml_io_stream_name<sensor::PencilBeamAntenna> {
  static constexpr std::string_view name = "SensorPencilBeamAntenna";
};

template <>
struct xml_io_stream_aggregate<sensor::PencilBeamAntenna> {
  static constexpr bool value = true;
};

// GaussianAntenna format tags and XML I/O

template <>
struct format_tag_aggregate<sensor::GaussianAntenna> {
  constexpr static bool value = true;
};

template <>
struct xml_io_stream_name<sensor::GaussianAntenna> {
  static constexpr std::string_view name = "SensorGaussianAntenna";
};

template <>
struct xml_io_stream_aggregate<sensor::GaussianAntenna> {
  static constexpr bool value = true;
};
