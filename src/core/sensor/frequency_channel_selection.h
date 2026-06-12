#pragma once

#include <format_tags.h>
#include <matpack.h>
#include <xml_io_stream.h>

namespace sensor {
//! Free-form channel struct.  Others inherit from this.
struct Channel;

//! A channel that is even between a lower and upper frequency, with equal weights.
struct BoxChannel;

//! A channel with a single frequency and all weight on that frequency.
struct DiracChannel;

//! A channel with Gaussian weights centered at the center frequency with given standard deviation.
struct GaussianChannel;

//! Concept for selecting frequencies for a channel
template <typename T>
concept FrequencyChannelSelection = std::derived_from<T, Channel>;

struct Channel {
  SortedGriddedField1 channel;

  [[nodiscard]] const AscendingGrid& freq_grid() const;
  [[nodiscard]] const Vector& weights() const;
  [[nodiscard]] Size size() const { return channel.data.size(); }
};

struct BoxChannel final : Channel {
  BoxChannel(Numeric lower, Numeric upper, Size N);  // [lower, upper]
  BoxChannel(Numeric hw, Size N);                    // [-hw, hw]
  BoxChannel(AscendingGrid f);                       // f
  BoxChannel() = default;                            // [0, 0]
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
  GaussianChannel() = default;                                // f0 = 0, std = 0

  [[nodiscard]] AscendingGrid center_grid() const;
};

struct Spectrometer {
  std::vector<Channel> channels;

 private:
  void sync_frequency_grids();

 public:
  [[nodiscard]] bool is_synced() const;

  /** Make a simple spectrometer

  The spectrometer will have the same channel response for each position
  but will be shifted by the frequencies in freq_offsets.

  So if base_channel is f: [-1, 1], w: [0.5, 0.5], and freq_offset is [2, 3, 4],
  the resulting channels will be [[1, 2], [2, 3], [3, 4]] with the weights [0.5, 0.5] for each channel.

  But, since this is supposed to simulate a spectrometer, the true frequency grid of each channel
  will be [1, 2, 3, 4] for all channels but the weights will be [0.5, 0.5, 0, 0], [0, 0.5, 0.5, 0], and
  [0, 0, 0.5, 0.5] for each channel respectively.

  @param[in] base_channel The base channel
  @param[in] freq_offsets The frequency offsets for each channel
*/
  Spectrometer(const Channel& base_channel, const AscendingGrid& freq_offsets);

  Spectrometer(const AscendingGrid& freq_offsets);

  Spectrometer(std::vector<Channel> channels) : channels(std::move(channels)) {}

  operator const std::vector<Channel>&() const { return channels; }
};
}  // namespace sensor

// Channel format tags and XML I/O

template <>
struct format_tag_aggregate<sensor::Channel> {
  constexpr static bool value = true;
};

template <>
struct xml_io_stream_name<sensor::Channel> {
  static constexpr std::string_view name = "SensorChannel";
};

template <>
struct xml_io_stream_aggregate<sensor::Channel> {
  static constexpr bool value = true;
};

// BoxChannel format tags and XML I/O

template <>
struct std::formatter<sensor::BoxChannel>
    : format_tag_inherit<sensor::Channel, sensor::BoxChannel> {};

template <>
struct xml_io_stream<sensor::BoxChannel>
    : xml_io_stream_inherit<sensor::Channel, sensor::BoxChannel> {};

// DiracChannel format tags and XML I/O

template <>
struct std::formatter<sensor::DiracChannel>
    : format_tag_inherit<sensor::Channel, sensor::DiracChannel> {};

template <>
struct xml_io_stream<sensor::DiracChannel>
    : xml_io_stream_inherit<sensor::Channel, sensor::DiracChannel> {};

// GaussianChannel format tags and XML I/O

template <>
struct std::formatter<sensor::GaussianChannel>
    : format_tag_inherit<sensor::Channel, sensor::GaussianChannel> {};

template <>
struct xml_io_stream<sensor::GaussianChannel>
    : xml_io_stream_inherit<sensor::Channel, sensor::GaussianChannel> {};

// Spectrometer format tags and XML I/O

template <>
struct std::formatter<sensor::Spectrometer>
    : std::formatter<std::vector<sensor::Channel>> {};

template <>
struct xml_io_stream<sensor::Spectrometer> {
  static void write(std::ostream&,
                    const sensor::Spectrometer&,
                    bofstream*       = nullptr,
                    std::string_view = ""sv);
  static void read(std::istream&, sensor::Spectrometer&, bifstream* = nullptr);
};
