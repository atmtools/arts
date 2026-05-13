#pragma once

#include <matpack.h>
#include <rtepack.h>

#include <concepts>
#include <memory>

#include "frequency_channel_selection.h"
#include "obsel.h"

namespace sensor {
//! A 2D angular antenna pattern on local zenith and azimuth offsets.
struct AntennaPattern;

//! A 2D gridded antenna pattern on local zenith and azimuth offsets.
struct GriddedAntennaPattern;

//! A 1x1 antenna pattern that samples only the bore line of sight.
struct PencilBeamAntenna;

//! A 2D Gaussian antenna pattern on local zenith and azimuth offsets.
struct GaussianAntenna;

//! A Gaussianized Airy antenna pattern with frequency-dependent width.
struct GaussianAiryAntenna;

//! 2D gridded field of antenna weights on local zenith and azimuth offsets.
using AntennaPatternField = matpack::gridded_data_t<Stokvec, ZenGrid, AziGrid>;

struct AntennaPattern {
  virtual ~AntennaPattern()                        = default;
  AntennaPattern()                                 = default;
  AntennaPattern(const AntennaPattern&)            = default;
  AntennaPattern& operator=(const AntennaPattern&) = default;
  AntennaPattern(AntennaPattern&&)                 = default;
  AntennaPattern& operator=(AntennaPattern&&)      = default;

  //! Builds one sensor obsel for a channel at a sensor position and bore LOS.
  [[nodiscard]] virtual Obsel operator()(const Channel& channel,
                                         const Vector3& pos,
                                         const Vector2& bore_los) const = 0;

  //! Creates an owning copy preserving the dynamic antenna type.
  [[nodiscard]] virtual std::shared_ptr<const AntennaPattern> clone() const = 0;
};

struct GriddedAntennaPattern : AntennaPattern {
  AntennaPatternField data;  // center at [0, 0]

  [[nodiscard]] Obsel operator()(const Channel& channel,
                                 const Vector3& pos,
                                 const Vector2& bore_los) const override;

  [[nodiscard]] std::shared_ptr<const AntennaPattern> clone() const override;
};

struct PencilBeamAntenna final : AntennaPattern {
  Stokvec weight{1.0, 0.0, 0.0, 0.0};

  PencilBeamAntenna(Stokvec weight = {1.0, 0.0, 0.0, 0.0});

  [[nodiscard]] Obsel operator()(const Channel& channel,
                                 const Vector3& pos,
                                 const Vector2& bore_los) const override;

  [[nodiscard]] std::shared_ptr<const AntennaPattern> clone() const override;
};

struct GaussianAntenna final : GriddedAntennaPattern {
  GaussianAntenna()                                  = default;
  GaussianAntenna(const GaussianAntenna&)            = default;
  GaussianAntenna& operator=(const GaussianAntenna&) = default;
  GaussianAntenna(GaussianAntenna&&)                 = default;
  GaussianAntenna& operator=(GaussianAntenna&&)      = default;

  GaussianAntenna(ZenGrid zen_grid,
                  AziGrid azi_grid,
                  Numeric zenith_std,
                  Numeric azimuth_std,
                  Stokvec weight = {1.0, 0.0, 0.0, 0.0});

  GaussianAntenna(ZenGrid zen_grid,
                  AziGrid azi_grid,
                  Numeric std,
                  Stokvec weight = {1.0, 0.0, 0.0, 0.0});

  [[nodiscard]] std::shared_ptr<const AntennaPattern> clone() const override;
};

struct GaussianAiryAntenna final : AntennaPattern {
  ZenGrid zen_grid;
  AziGrid azi_grid;
  Numeric aperture_diameter;
  Stokvec weight;

  GaussianAiryAntenna()                                      = default;
  GaussianAiryAntenna(const GaussianAiryAntenna&)            = default;
  GaussianAiryAntenna& operator=(const GaussianAiryAntenna&) = default;
  GaussianAiryAntenna(GaussianAiryAntenna&&)                 = default;
  GaussianAiryAntenna& operator=(GaussianAiryAntenna&&)      = default;

  GaussianAiryAntenna(ZenGrid zen_grid,
                      AziGrid azi_grid,
                      Numeric aperture_diameter,
                      Stokvec weight = {1.0, 0.0, 0.0, 0.0});

  [[nodiscard]] Obsel operator()(const Channel& channel,
                                 const Vector3& pos,
                                 const Vector2& bore_los) const override;

  [[nodiscard]] std::shared_ptr<const AntennaPattern> clone() const override;
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
  constexpr static bool value = false;
};

template <>
struct xml_io_stream_name<sensor::AntennaPattern> {
  static constexpr std::string_view name = "SensorAntennaPattern";
};

template <>
struct xml_io_stream_aggregate<sensor::AntennaPattern> {
  static constexpr bool value = false;
};

// GriddedAntennaPattern format tags and XML I/O

template <>
struct format_tag_aggregate<sensor::GriddedAntennaPattern> {
  constexpr static bool value = true;
};

template <>
struct xml_io_stream_name<sensor::GriddedAntennaPattern> {
  static constexpr std::string_view name = "SensorGriddedAntennaPattern";
};

template <>
struct xml_io_stream_aggregate<sensor::GriddedAntennaPattern> {
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

// GaussianAiryAntenna format tags and XML I/O

template <>
struct format_tag_aggregate<sensor::GaussianAiryAntenna> {
  constexpr static bool value = true;
};

template <>
struct xml_io_stream_name<sensor::GaussianAiryAntenna> {
  static constexpr std::string_view name = "SensorGaussianAiryAntenna";
};

template <>
struct xml_io_stream_aggregate<sensor::GaussianAiryAntenna> {
  static constexpr bool value = true;
};
