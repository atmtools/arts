#pragma once

#include <matpack.h>
#include <rtepack.h>

namespace sensor {
//! A 2D angular antenna pattern on local zenith and azimuth offsets.
struct AntennaPattern;

//! A 1x1 antenna pattern that samples only the bore line of sight.
struct PencilBeamAntenna;

//! A 2D Gaussian antenna pattern on local zenith and azimuth offsets.
struct GaussianAntenna;

//! 2D gridded field of antenna weights on local zenith and azimuth offsets.
using AntennaPatternField = matpack::gridded_data_t<Stokvec, ZenGrid, AziGrid>;

//! Concept for selecting antenna patterns.
template <typename T>
concept AntennaPatternSelection = std::derived_from<T, AntennaPattern>;

struct AntennaPattern {
  AntennaPatternField data;  // center at [0, 0]

  //! Maps the local antenna pattern to global LOS values around a new bore LOS.
  //!
  //! Local zenith 0 points along the bore. Local azimuth 0 points toward
  //! increasing global zenith, and local azimuth 90 points toward increasing
  //! global azimuth.
  [[nodiscard]] std::vector<std::pair<Stokvec, Vector2>> operator()(
      Vector2 bore_los) const;
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