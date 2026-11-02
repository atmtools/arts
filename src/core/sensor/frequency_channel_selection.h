#pragma once

#include <format_tags.h>
#include <matpack.h>
#include <xml_io_stream.h>

#include "sensor_polarization.h"

namespace sensor {
//! Free-form channel struct.  Others inherit from this.
struct Channel;

//! A channel that responds evenly between a lower and upper frequency.
struct BoxChannel;

//! A channel with a single frequency and all weight on that frequency.
struct DiracChannel;

//! A channel with Gaussian response centered at the center frequency with given standard deviation.
struct GaussianChannel;

//! Concept for selecting frequencies for a channel
template <typename T>
concept FrequencyChannelSelection = std::derived_from<T, Channel>;

/** A receiver channel: its frequency response and its polarization response
 *
 *  The response is a sum of pieces.  Each piece is linear between its
 *  frequencies and zero outside them, so the frequencies of a piece are both
 *  where the response is defined and where the radiance is sampled: their
 *  density is chosen with the channel.  A piece with a single frequency is a
 *  delta function of that weight.  Pieces may overlap, and a response with a
 *  gap, such as the two sidebands of a heterodyne channel, is two pieces.
 *
 *  The measurement is the integral of the response times the radiance over
 *  the integral of the response.  quadrature() gives the integration weights
 *  of that integral on the frequencies of the pieces: the trapezoid rule
 *  within each piece, delta weights as they are, summed where pieces share a
 *  frequency.
 */
struct Channel {
  std::vector<SortedGriddedField1> response;

  //! What the channel receives of the Stokes vector at the boresight
  Polarization polarization{};

  //! The integration weights of the response on the frequencies of its pieces, not normalized
  [[nodiscard]] SortedGriddedField1 quadrature() const;
};

struct BoxChannel final : Channel {
  BoxChannel(Numeric lower, Numeric upper, Size N);  // [lower, upper]
  BoxChannel(Numeric hw, Size N);                    // [-hw, hw]
  BoxChannel(AscendingGrid f);                       // f
  BoxChannel() = default;                            // no response
};

struct DiracChannel final : Channel {
  DiracChannel(Numeric f);  // f
  DiracChannel();           // f = 0
};

struct GaussianChannel final : Channel {
  GaussianChannel(AscendingGrid f, Numeric f0, Numeric std);  // std around f0
  GaussianChannel(Numeric f0, Numeric std, Size N, Size M);   // f0 +- M*std
  GaussianChannel(AscendingGrid f, Numeric std);              // f, f0 = 0
  GaussianChannel(Numeric std, Size N, Size M);               // +-M*std, f0 = 0
  GaussianChannel() = default;                                // no response
};

/** Channels that share a backend
 *
 *  A sensor built from a spectrometer puts all its channels on one common
 *  frequency grid.
 */
struct Spectrometer {
  std::vector<Channel> channels;

  /** Make a simple spectrometer

  Every channel is base_channel shifted by one of freq_offsets.  So if
  base_channel responds on [-1, 1] and freq_offsets is [2, 3, 4], the channels
  respond on [1, 3], [2, 4], and [3, 5].

  @param[in] base_channel The base channel
  @param[in] freq_offsets The frequency offsets for each channel
*/
  Spectrometer(const Channel& base_channel, const AscendingGrid& freq_offsets);

  //! Dirac channels at freq_offsets
  Spectrometer(const AscendingGrid& freq_offsets);

  Spectrometer(std::vector<Channel> channels) : channels(std::move(channels)) {}

  operator const std::vector<Channel>&() const { return channels; }
};
}  // namespace sensor

// Channel format tags and XML I/O

template <> struct format_tag_aggregate<sensor::Channel> {
  constexpr static bool value = true;
};

template <> struct xml_io_stream_name<sensor::Channel> {
  static constexpr std::string_view name = "SensorChannel";
};

template <> struct xml_io_stream_aggregate<sensor::Channel> {
  static constexpr bool value = true;
};

// BoxChannel format tags and XML I/O

template <> struct std::formatter<sensor::BoxChannel> : format_tag_inherit<sensor::Channel, sensor::BoxChannel> {};

template <> struct xml_io_stream<sensor::BoxChannel> : xml_io_stream_inherit<sensor::Channel, sensor::BoxChannel> {};

// DiracChannel format tags and XML I/O

template <> struct std::formatter<sensor::DiracChannel> : format_tag_inherit<sensor::Channel, sensor::DiracChannel> {};

template <> struct xml_io_stream<sensor::DiracChannel> : xml_io_stream_inherit<sensor::Channel, sensor::DiracChannel> {
};

// GaussianChannel format tags and XML I/O

template <> struct std::formatter<sensor::GaussianChannel>
    : format_tag_inherit<sensor::Channel, sensor::GaussianChannel> {};

template <> struct xml_io_stream<sensor::GaussianChannel>
    : xml_io_stream_inherit<sensor::Channel, sensor::GaussianChannel> {};

// Spectrometer format tags and XML I/O

template <> struct std::formatter<sensor::Spectrometer> : std::formatter<std::vector<sensor::Channel>> {};

template <> struct xml_io_stream<sensor::Spectrometer> {
  static void write(std::ostream&, const sensor::Spectrometer&, bofstream* = nullptr, std::string_view = ""sv);
  static void read(std::istream&, sensor::Spectrometer&, bifstream* = nullptr);
};
