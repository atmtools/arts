#pragma once

#include <format_tags.h>
#include <matpack.h>
#include <rtepack.h>
#include <xml.h>

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
                                         const Vector2& bore_los,
                                         const Vector2& ell) const = 0;

  //! Creates an owning copy preserving the dynamic antenna type.
  [[nodiscard]] virtual std::shared_ptr<const AntennaPattern> clone() const = 0;
};

struct PencilBeamAntenna final : AntennaPattern {
  Stokvec weight{1.0, 0.0, 0.0, 0.0};

  PencilBeamAntenna(Stokvec weight = {1.0, 0.0, 0.0, 0.0});

  [[nodiscard]] Obsel operator()(const Channel& channel,
                                 const Vector3& pos,
                                 const Vector2& bore_los,
                                 const Vector2& ell) const override;

  [[nodiscard]] std::shared_ptr<const AntennaPattern> clone() const override;
};

//! 2D gridded field of antenna weights on local zenith and azimuth offsets.
using AntennaPatternField = matpack::gridded_data_t<Numeric, ZenGrid, AziGrid>;

struct GriddedAntennaPattern : AntennaPattern {
  AntennaPatternField data;  // center at [0, 0]
  Stokvec weight{1.0, 0.0, 0.0, 0.0};

  [[nodiscard]] Obsel operator()(const Channel& channel,
                                 const Vector3& pos,
                                 const Vector2& bore_los,
                                 const Vector2& ell) const override;

  [[nodiscard]] std::shared_ptr<const AntennaPattern> clone() const override;
};

struct GaussianAntenna final : GriddedAntennaPattern {
  GaussianAntenna()                                  = default;
  GaussianAntenna(const GaussianAntenna&)            = default;
  GaussianAntenna& operator=(const GaussianAntenna&) = default;
  GaussianAntenna(GaussianAntenna&&)                 = default;
  GaussianAntenna& operator=(GaussianAntenna&&)      = default;

  GaussianAntenna(ZenGrid grid,
                  Numeric std,
                  Size azi_grid_size = 2,
                  Stokvec weight     = {1.0, 0.0, 0.0, 0.0});

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
                      Numeric aperture_diameter,
                      Size azi_grid_size = 2,
                      Stokvec weight     = {1.0, 0.0, 0.0, 0.0});

  [[nodiscard]] Obsel operator()(const Channel& channel,
                                 const Vector3& pos,
                                 const Vector2& bore_los,
                                 const Vector2& ell) const override;

  [[nodiscard]] std::shared_ptr<const AntennaPattern> clone() const override;
  [[nodiscard]] Numeric std(Numeric frequency) const;
};

//! Concept for types that can build a sensor obsel from channel geometry.
template <typename T>
concept AntennaPatternSelection = requires(const T& antenna,
                                           const Channel& channel,
                                           const Vector3& pos,
                                           const Vector2& bore_los,
                                           const Vector2& ell) {
  { antenna(channel, pos, bore_los, ell) } -> std::same_as<Obsel>;
};
}  // namespace sensor

// AntennaPattern format tags and XML I/O

template <>
struct std::formatter<sensor::AntennaPattern> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const sensor::AntennaPattern&,
                              FmtContext& ctx) const {
    return tags.format(ctx, "SensorAntennaPattern"sv);
  }
};

template <>
struct xml_io_stream_name<sensor::AntennaPattern> {
  static constexpr std::string_view name = "SensorAntennaPattern";
};

template <>
struct xml_io_stream<sensor::AntennaPattern> {
  static constexpr std::string_view type_name =
      xml_io_stream_name_v<sensor::AntennaPattern>;

  static void write(std::ostream& os,
                    const sensor::AntennaPattern& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   sensor::AntennaPattern& n,
                   bifstream* pbifs = nullptr);
};

// GriddedAntennaPattern format tags and XML I/O

template <>
struct std::formatter<sensor::GriddedAntennaPattern> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const sensor::GriddedAntennaPattern& v,
                              FmtContext& ctx) const {
    return tags.format(ctx, v.data);
  }
};

template <>
struct xml_io_stream_name<sensor::GriddedAntennaPattern> {
  static constexpr std::string_view name = "SensorGriddedAntennaPattern";
};

template <>
struct xml_io_stream<sensor::GriddedAntennaPattern> {
  static constexpr std::string_view type_name =
      xml_io_stream_name_v<sensor::GriddedAntennaPattern>;

  static void write(std::ostream& os,
                    const sensor::GriddedAntennaPattern& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   sensor::GriddedAntennaPattern& n,
                   bifstream* pbifs = nullptr);
};

// PencilBeamAntenna format tags and XML I/O

template <>
struct std::formatter<sensor::PencilBeamAntenna> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const sensor::PencilBeamAntenna& v,
                              FmtContext& ctx) const {
    return tags.format(ctx, v.weight);
  }
};

template <>
struct xml_io_stream_name<sensor::PencilBeamAntenna> {
  static constexpr std::string_view name = "SensorPencilBeamAntenna";
};

template <>
struct xml_io_stream<sensor::PencilBeamAntenna> {
  static constexpr std::string_view type_name =
      xml_io_stream_name_v<sensor::PencilBeamAntenna>;

  static void write(std::ostream& os,
                    const sensor::PencilBeamAntenna& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   sensor::PencilBeamAntenna& n,
                   bifstream* pbifs = nullptr);
};

// GaussianAntenna format tags and XML I/O

template <>
struct std::formatter<sensor::GaussianAntenna>
    : format_tag_inherit<sensor::GriddedAntennaPattern,
                         sensor::GaussianAntenna> {};

template <>
struct xml_io_stream_name<sensor::GaussianAntenna> {
  static constexpr std::string_view name = "SensorGaussianAntenna";
};

template <>
struct xml_io_stream<sensor::GaussianAntenna>
    : xml_io_stream_inherit<sensor::GriddedAntennaPattern,
                            sensor::GaussianAntenna> {};

// GaussianAiryAntenna format tags and XML I/O

template <>
struct std::formatter<sensor::GaussianAiryAntenna> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const sensor::GaussianAiryAntenna& v,
                              FmtContext& ctx) const {
    auto sep = tags.sep();
    return tags.format(ctx,
                       v.zen_grid,
                       sep,
                       v.azi_grid,
                       sep,
                       v.aperture_diameter,
                       sep,
                       v.weight);
  }
};

template <>
struct xml_io_stream_name<sensor::GaussianAiryAntenna> {
  static constexpr std::string_view name = "SensorGaussianAiryAntenna";
};

template <>
struct xml_io_stream<sensor::GaussianAiryAntenna> {
  static constexpr std::string_view type_name =
      xml_io_stream_name_v<sensor::GaussianAiryAntenna>;

  static void write(std::ostream& os,
                    const sensor::GaussianAiryAntenna& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   sensor::GaussianAiryAntenna& n,
                   bifstream* pbifs = nullptr);
};
