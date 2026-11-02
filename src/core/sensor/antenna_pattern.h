#pragma once

#include <format_tags.h>
#include <matpack.h>
#include <rtepack.h>
#include <xml.h>

#include <concepts>
#include <variant>

#include "frequency_channel_selection.h"
#include "obsel.h"

namespace sensor {
//! A 2D gridded antenna pattern on local zenith and azimuth offsets.
struct GriddedAntennaPattern;

//! A 1x1 antenna pattern that samples only the bore line of sight.
struct PencilBeamAntenna;

//! A 2D Gaussian antenna pattern on local zenith and azimuth offsets.
struct GaussianAntenna;

//! A Gaussianized Airy antenna pattern with frequency-dependent width.
struct GaussianAiryAntenna;

struct PencilBeamAntenna {
  [[nodiscard]] Obsel operator()(const Channel& channel,
                                 const Stokvec& response,
                                 const Vector3& pos,
                                 const Vector2& bore_los,
                                 const Vector2& ell) const;
};

//! 2D gridded field of antenna gain on local zenith and azimuth offsets.
using AntennaPatternField = matpack::gridded_data_t<Numeric, ZenGrid, AziGrid>;

/** An antenna gain sampled on rings of zenith offset around the boresight
 *
 *  Each sample stands for the cell of the sphere between the midpoints to
 *  its neighbouring zenith and azimuth offsets, the outermost ring extending
 *  as far beyond its offset as half the distance to its inner neighbour.
 *  The beam weight of a sample is its gain times the solid angle of its cell.
 *  The ring at zero zenith offset is one sample, the cap inside the first
 *  midpoint, with the azimuthal mean of its gains.
 */
struct GriddedAntennaPattern {
  AntennaPatternField data;  // center at [0, 0]

  [[nodiscard]] Obsel operator()(const Channel& channel,
                                 const Stokvec& response,
                                 const Vector3& pos,
                                 const Vector2& bore_los,
                                 const Vector2& ell) const;
};

//! A gain exp(-(offset / std)^2 / 2) of the angular offset [deg] from the boresight, sampled as GriddedAntennaPattern
struct GaussianAntenna {
  ZenGrid zen_grid;
  AziGrid azi_grid;
  Numeric std{};  // [deg]

  GaussianAntenna()                                  = default;
  GaussianAntenna(const GaussianAntenna&)            = default;
  GaussianAntenna& operator=(const GaussianAntenna&) = default;
  GaussianAntenna(GaussianAntenna&&)                 = default;
  GaussianAntenna& operator=(GaussianAntenna&&)      = default;

  GaussianAntenna(ZenGrid zen_grid, Numeric std, Size azi_grid_size = 2);

  [[nodiscard]] Obsel operator()(const Channel& channel,
                                 const Stokvec& response,
                                 const Vector3& pos,
                                 const Vector2& bore_los,
                                 const Vector2& ell) const;
};

//! A Gaussian gain whose half width scales with wavelength over aperture diameter, sampled as GriddedAntennaPattern
struct GaussianAiryAntenna {
  ZenGrid zen_grid;
  AziGrid azi_grid;
  Numeric aperture_diameter;

  GaussianAiryAntenna()                                      = default;
  GaussianAiryAntenna(const GaussianAiryAntenna&)            = default;
  GaussianAiryAntenna& operator=(const GaussianAiryAntenna&) = default;
  GaussianAiryAntenna(GaussianAiryAntenna&&)                 = default;
  GaussianAiryAntenna& operator=(GaussianAiryAntenna&&)      = default;

  GaussianAiryAntenna(ZenGrid zen_grid, Numeric aperture_diameter, Size azi_grid_size = 2);

  [[nodiscard]] Obsel operator()(const Channel& channel,
                                 const Stokvec& response,
                                 const Vector3& pos,
                                 const Vector2& bore_los,
                                 const Vector2& ell) const;

  [[nodiscard]] Numeric std(Numeric frequency) const;
};

/** Any antenna pattern
 *
 *  Each pattern builds one sensor obsel for a channel at a sensor position and
 *  bore LOS, given the channel's polarization response at the boresight (see
 *  Polarization::at).  The weights are the channel's integration weights
 *  normalized to sum to one, times the beam weight of each line of sight
 *  normalized to sum to one at each frequency, times the polarization response
 *  carried from the boresight to that line of sight.
 */
using AntennaPattern = std::variant<PencilBeamAntenna, GaussianAntenna, GaussianAiryAntenna, GriddedAntennaPattern>;

//! The obsel of a channel through an antenna pattern, see AntennaPattern
[[nodiscard]] Obsel antenna_obsel(const AntennaPattern& antenna,
                                  const Channel&        channel,
                                  const Stokvec&        response,
                                  const Vector3&        pos,
                                  const Vector2&        bore_los,
                                  const Vector2&        ell);

//! Concept for types that can build a sensor obsel from channel geometry.
template <typename T>
concept AntennaPatternSelection = requires(const T&       antenna,
                                           const Channel& channel,
                                           const Stokvec& response,
                                           const Vector3& pos,
                                           const Vector2& bore_los,
                                           const Vector2& ell) {
  { antenna(channel, response, pos, bore_los, ell) } -> std::same_as<Obsel>;
};
}  // namespace sensor

// GriddedAntennaPattern format tags and XML I/O

template <> struct std::formatter<sensor::GriddedAntennaPattern> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const sensor::GriddedAntennaPattern& v, FmtContext& ctx) const {
    return tags.format(ctx, v.data);
  }
};

template <> struct xml_io_stream_name<sensor::GriddedAntennaPattern> {
  static constexpr std::string_view name = "SensorGriddedAntennaPattern";
};

template <> struct xml_io_stream<sensor::GriddedAntennaPattern> {
  static constexpr std::string_view type_name = xml_io_stream_name_v<sensor::GriddedAntennaPattern>;

  static void write(std::ostream&                        os,
                    const sensor::GriddedAntennaPattern& n,
                    bofstream*                           pbofs = nullptr,
                    std::string_view                     name  = ""sv);

  static void read(std::istream& is, sensor::GriddedAntennaPattern& n, bifstream* pbifs = nullptr);
};

// PencilBeamAntenna format tags and XML I/O

template <> struct std::formatter<sensor::PencilBeamAntenna> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const sensor::PencilBeamAntenna&, FmtContext& ctx) const {
    return tags.format(ctx, "SensorPencilBeamAntenna"sv);
  }
};

template <> struct xml_io_stream_name<sensor::PencilBeamAntenna> {
  static constexpr std::string_view name = "SensorPencilBeamAntenna";
};

template <> struct xml_io_stream<sensor::PencilBeamAntenna> {
  static constexpr std::string_view type_name = xml_io_stream_name_v<sensor::PencilBeamAntenna>;

  static void write(std::ostream&                    os,
                    const sensor::PencilBeamAntenna& n,
                    bofstream*                       pbofs = nullptr,
                    std::string_view                 name  = ""sv);

  static void read(std::istream& is, sensor::PencilBeamAntenna& n, bifstream* pbifs = nullptr);
};

// GaussianAntenna format tags and XML I/O

template <> struct std::formatter<sensor::GaussianAntenna> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const sensor::GaussianAntenna& v, FmtContext& ctx) const {
    auto sep = tags.sep();
    return tags.format(ctx, v.zen_grid, sep, v.azi_grid, sep, v.std);
  }
};

template <> struct xml_io_stream_name<sensor::GaussianAntenna> {
  static constexpr std::string_view name = "SensorGaussianAntenna";
};

template <> struct xml_io_stream<sensor::GaussianAntenna> {
  static constexpr std::string_view type_name = xml_io_stream_name_v<sensor::GaussianAntenna>;

  static void write(std::ostream&                  os,
                    const sensor::GaussianAntenna& n,
                    bofstream*                     pbofs = nullptr,
                    std::string_view               name  = ""sv);

  static void read(std::istream& is, sensor::GaussianAntenna& n, bifstream* pbifs = nullptr);
};

// GaussianAiryAntenna format tags and XML I/O

template <> struct std::formatter<sensor::GaussianAiryAntenna> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const sensor::GaussianAiryAntenna& v, FmtContext& ctx) const {
    auto sep = tags.sep();
    return tags.format(ctx, v.zen_grid, sep, v.azi_grid, sep, v.aperture_diameter);
  }
};

template <> struct xml_io_stream_name<sensor::GaussianAiryAntenna> {
  static constexpr std::string_view name = "SensorGaussianAiryAntenna";
};

template <> struct xml_io_stream<sensor::GaussianAiryAntenna> {
  static constexpr std::string_view type_name = xml_io_stream_name_v<sensor::GaussianAiryAntenna>;

  static void write(std::ostream&                      os,
                    const sensor::GaussianAiryAntenna& n,
                    bofstream*                         pbofs = nullptr,
                    std::string_view                   name  = ""sv);

  static void read(std::istream& is, sensor::GaussianAiryAntenna& n, bifstream* pbifs = nullptr);
};
