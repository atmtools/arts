#pragma once

#include <format_tags.h>
#include <matpack.h>
#include <rtepack.h>
#include <xml.h>

#include <concepts>
#include <span>
#include <variant>

#include "frequency_channel_selection.h"
#include "sensor.h"

namespace sensor {
/** The layout of an antenna's pixels
 *
 *  An antenna makes one observation element per pixel for each channel.  Its
 *  pixels form a grid of rows and columns, given by the angular offsets [deg]
 *  of their boresights from the antenna's boresight, and its elements are
 *  ordered by row and then by column.  A kind with more than one pixel defines
 *  the directions of its rows and columns.  An antenna that is a single beam
 *  has one pixel at no offset.
 */
struct PixelGrid {
  AscendingGrid rows{Vector{0.0}};
  AscendingGrid cols{Vector{0.0}};

  //! The number of pixels
  [[nodiscard]] Size size() const;
};

//! A 2D gridded antenna pattern on local zenith and azimuth offsets.
struct GriddedAntennaPattern;

//! A 1x1 antenna pattern that samples only the bore line of sight.
struct PencilBeamAntenna;

//! A 2D Gaussian antenna pattern on local zenith and azimuth offsets.
struct GaussianAntenna;

//! A Gaussianized Airy antenna pattern with frequency-dependent width.
struct GaussianAiryAntenna;

//! A 2D Gaussian antenna pattern on a rectangular grid of angular offsets.
struct GaussianRectangularAntenna;

//! A camera: a grid of pixels, each with the same beam.
struct CameraAntenna;

struct PencilBeamAntenna {
  //! One pixel at the boresight
  [[nodiscard]] PixelGrid pixels() const;

  [[nodiscard]] Sensor operator()(std::span<const Channel> channels,
                                  std::span<const Stokvec> responses,
                                  const Vector3&           pos,
                                  const Vector2&           bore_los,
                                  const Vector2&           ell) const;
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

  //! One pixel at the boresight
  [[nodiscard]] PixelGrid pixels() const;

  [[nodiscard]] Sensor operator()(std::span<const Channel> channels,
                                  std::span<const Stokvec> responses,
                                  const Vector3&           pos,
                                  const Vector2&           bore_los,
                                  const Vector2&           ell) const;
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

  //! One pixel at the boresight
  [[nodiscard]] PixelGrid pixels() const;

  [[nodiscard]] Sensor operator()(std::span<const Channel> channels,
                                  std::span<const Stokvec> responses,
                                  const Vector3&           pos,
                                  const Vector2&           bore_los,
                                  const Vector2&           ell) const;
};

/** A gain exp(-(offset / std)^2 / 2) of the angular offset [deg] from the boresight, sampled on a rectangular grid
 *
 *  The samples are at all pairs of the angular offsets zen_offsets [deg]
 *  towards increasing zenith angle and azi_offsets [deg] towards increasing
 *  azimuth at the boresight, by zenith offset and then azimuth offset.  A pair
 *  (z, a) looks along the boresight plus tan(z) and tan(a) times the unit
 *  vectors of increasing zenith and azimuth there.  Each sample stands for the
 *  cell between the midpoints to its neighbouring offsets, the outermost cells
 *  extending as far beyond their offsets as half the distance to their inner
 *  neighbours, and its beam weight is its gain times the solid angle of its
 *  cell.  Unlike the rings of GaussianAntenna, the grid can be fine in one
 *  direction only, such as across a limb.
 */
struct GaussianRectangularAntenna {
  AscendingGrid zen_offsets{};
  AscendingGrid azi_offsets{};
  Numeric       std{};  // [deg]

  //! One pixel at the boresight
  [[nodiscard]] PixelGrid pixels() const;

  [[nodiscard]] Sensor operator()(std::span<const Channel> channels,
                                  std::span<const Stokvec> responses,
                                  const Vector3&           pos,
                                  const Vector2&           bore_los,
                                  const Vector2&           ell) const;
};

/** A Gaussian gain with the half-power width of the Airy pattern of a circular aperture, sampled as GriddedAntennaPattern
 *
 *  The Airy pattern (2 J1(x) / x)^2, x = pi D sin(offset) / wavelength, is that of a uniformly illuminated
 *  aperture of diameter D [m].  Its half-power half width for small angles, Constant::airy_half_power_argument
 *  / pi times wavelength over D [rad], is the half width at half maximum of the Gaussian at each frequency.
 */
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

  //! One pixel at the boresight
  [[nodiscard]] PixelGrid pixels() const;

  [[nodiscard]] Sensor operator()(std::span<const Channel> channels,
                                  std::span<const Stokvec> responses,
                                  const Vector3&           pos,
                                  const Vector2&           bore_los,
                                  const Vector2&           ell) const;

  [[nodiscard]] Numeric std(Numeric frequency) const;
};

//! The kinds of antenna pattern that are a single beam
using BeamPattern = std::
    variant<PencilBeamAntenna, GaussianAntenna, GaussianAiryAntenna, GriddedAntennaPattern, GaussianRectangularAntenna>;

/** A camera: a grid of pixels in CCD coordinates, each with the same beam
 *
 *  The pixel grid holds the angular offsets [deg] of the pixels from the
 *  optical axis, the boresight, in CCD coordinates.  A pixel at row offset r
 *  and column offset c looks along g + tan(r) e_zen - tan(c) e_azi, where g is
 *  the boresight and e_zen and e_azi point towards increasing zenith and
 *  azimuth at the boresight; the lens inverts the image on the CCD.  Each
 *  pixel's polarization response is the boresight response carried along the
 *  great circle to the pixel, and its beam is centred on the pixel.
 */
struct CameraAntenna {
  PixelGrid   grid;
  BeamPattern beam{};

  /** A thin-lens camera
   *
   *  n_rows by n_cols pixels on a CCD of ccd_height by ccd_width [m], behind a
   *  lens of focal_length [m] focused at focus_distance [m].  The CCD sits at
   *  the image distance of the thin-lens equation, and each pixel's offset is
   *  the angle of its centre seen from the lens.
   */
  [[nodiscard]] static CameraAntenna thin_lens(Size        n_rows,
                                               Size        n_cols,
                                               Numeric     ccd_height,
                                               Numeric     ccd_width,
                                               Numeric     focal_length,
                                               Numeric     focus_distance,
                                               BeamPattern beam = {});

  [[nodiscard]] PixelGrid pixels() const;

  [[nodiscard]] Sensor operator()(std::span<const Channel> channels,
                                  std::span<const Stokvec> responses,
                                  const Vector3&           pos,
                                  const Vector2&           bore_los,
                                  const Vector2&           ell) const;
};

/** Any antenna pattern
 *
 *  An antenna pattern describes the layout of its pixels, see PixelGrid, and
 *  the beam of each pixel.  For channels at a sensor position and bore LOS,
 *  given each channel's polarization response at the boresight (see
 *  Polarization::at), it builds a sensor with one element per pixel and
 *  channel, by pixel and then channel, on the frequencies of all channels and
 *  every line of sight its beams sample.  The lines of sight do not depend on
 *  the channel, so the antenna is laid out once for all channels.  The weights
 *  of an element are the channel's integration weights normalized to sum to
 *  one, times the beam weight of each line of sight normalized to sum to one at
 *  each frequency, times the polarization response carried from the boresight
 *  to that line of sight.  The kinds in BeamPattern are a single beam, one
 *  pixel at the boresight; a CameraAntenna has a grid of pixels.
 */
using AntennaPattern = std::variant<PencilBeamAntenna,
                                    GaussianAntenna,
                                    GaussianAiryAntenna,
                                    GriddedAntennaPattern,
                                    CameraAntenna,
                                    GaussianRectangularAntenna>;

//! The layout of the pixels of an antenna pattern, see PixelGrid
[[nodiscard]] PixelGrid antenna_pixels(const AntennaPattern& antenna);

//! The sensor of channels through an antenna pattern, one element per pixel and channel, see AntennaPattern
[[nodiscard]] Sensor antenna_sensor(const AntennaPattern&    antenna,
                                    std::span<const Channel> channels,
                                    std::span<const Stokvec> responses,
                                    const Vector3&           pos,
                                    const Vector2&           bore_los,
                                    const Vector2&           ell);

//! The sensor of one channel through an antenna pattern, one element per pixel by row and then column
[[nodiscard]] Sensor antenna_sensor(const AntennaPattern& antenna,
                                    const Channel&        channel,
                                    const Stokvec&        response,
                                    const Vector3&        pos,
                                    const Vector2&        bore_los,
                                    const Vector2&        ell);

//! Concept for antenna patterns: the layout of their pixels, and a sensor with one element per pixel and channel
template <typename T>
concept AntennaPatternSelection = requires(const T&                 antenna,
                                           std::span<const Channel> channels,
                                           std::span<const Stokvec> responses,
                                           const Vector3&           pos,
                                           const Vector2&           bore_los,
                                           const Vector2&           ell) {
  { antenna.pixels() } -> std::same_as<PixelGrid>;
  { antenna(channels, responses, pos, bore_los, ell) } -> std::same_as<Sensor>;
};
}  // namespace sensor

template <> struct format_tag_aggregate<sensor::PixelGrid> {
  constexpr static bool value = true;
};

template <> struct xml_io_stream_name<sensor::PixelGrid> {
  static constexpr std::string_view name = "SensorPixelGrid";
};

template <> struct xml_io_stream_aggregate<sensor::PixelGrid> {
  static constexpr bool value = true;
};

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

// GaussianRectangularAntenna format tags and XML I/O

template <> struct std::formatter<sensor::GaussianRectangularAntenna> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const sensor::GaussianRectangularAntenna& v, FmtContext& ctx) const {
    auto sep = tags.sep();
    return tags.format(ctx, v.zen_offsets, sep, v.azi_offsets, sep, v.std);
  }
};

template <> struct xml_io_stream_name<sensor::GaussianRectangularAntenna> {
  static constexpr std::string_view name = "SensorGaussianRectangularAntenna";
};

template <> struct xml_io_stream<sensor::GaussianRectangularAntenna> {
  static constexpr std::string_view type_name = xml_io_stream_name_v<sensor::GaussianRectangularAntenna>;

  static void write(std::ostream&                             os,
                    const sensor::GaussianRectangularAntenna& n,
                    bofstream*                                pbofs = nullptr,
                    std::string_view                          name  = ""sv);

  static void read(std::istream& is, sensor::GaussianRectangularAntenna& n, bifstream* pbifs = nullptr);
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

// CameraAntenna format tags and XML I/O

template <> struct format_tag_aggregate<sensor::CameraAntenna> {
  constexpr static bool value = true;
};

template <> struct xml_io_stream_name<sensor::CameraAntenna> {
  static constexpr std::string_view name = "SensorCameraAntenna";
};

template <> struct xml_io_stream_aggregate<sensor::CameraAntenna> {
  static constexpr bool value = true;
};
