#pragma once

#include <matpack.h>

#include <limits>
#include <optional>

namespace sensor {
//! Free-form frequency range struct.  Others inherit from this.
struct FrequencyRange;

//! A frequency range for heterodyne channels.  True frequencies are computed.
struct HeterodyneFrequencyRange;

//! A single affine path through a staged heterodyne chain.
struct FrequencyResponsePath;

//! Concept for selecting frequency ranges for a set of channels
template <typename T>
concept FrequencyRangeSelection = std::derived_from<T, FrequencyRange>;

struct FrequencyResponsePath {
  static constexpr Numeric inf = std::numeric_limits<Numeric>::infinity();

  Numeric intercept{0.0};
  Numeric slope{1.0};
  Vector2 local_range{0.0, inf};
  std::vector<SortedGriddedField1> filters{};

  [[nodiscard]] Vector2 global_range() const;
  [[nodiscard]] Numeric map_to_global(Numeric local_frequency) const;
  [[nodiscard]] std::optional<Numeric> map_to_local(
      Numeric global_frequency) const;
  [[nodiscard]] Numeric local_weight(Numeric local_frequency) const;
  [[nodiscard]] Numeric global_weight(Numeric global_frequency) const;
};

struct FrequencyRange {
  static constexpr Numeric inf = FrequencyResponsePath::inf;

  std::vector<FrequencyResponsePath> response_paths{{}};

  std::vector<Vector2> global_ranges{{0, inf}};
  std::vector<Vector2> local_ranges{{0, inf}};

  [[nodiscard]] Size size() const;
  [[nodiscard]] const FrequencyResponsePath& path(Size index) const;
  [[nodiscard]] const std::vector<FrequencyResponsePath>& paths() const;

  void sync_ranges();
};

struct HeterodyneFrequencyRange final : FrequencyRange {
  HeterodyneFrequencyRange() = default;
  HeterodyneFrequencyRange(Numeric clock_frequency,
                           const Vector2& bandpass_range);
  HeterodyneFrequencyRange(const std::span<const Numeric>& clock_frequencies,
                           const std::span<const Vector2>& bandpass_ranges);

  void apply_lowpass(Numeric upper_frequency);
  void apply_highpass(Numeric lower_frequency);
  void apply_bandpass(const Vector2& bandpass_range);
  void apply_bandpass(const SortedGriddedField1& bandpass_filter);
  void apply_mixer(Numeric clock_frequency);

  [[nodiscard]] Vector local_response(ConstVectorView local_frequency_grid,
                                      Size path_index = 0) const;
  [[nodiscard]] Vector global_response(ConstVectorView global_frequency_grid,
                                       Size path_index = 0) const;
};
}  // namespace sensor

// FrequencyResponsePath format tags and XML I/O

template <>
struct format_tag_aggregate<sensor::FrequencyResponsePath> {
  constexpr static bool value = true;
};

template <>
struct xml_io_stream_name<sensor::FrequencyResponsePath> {
  static constexpr std::string_view name = "SensorFrequencyResponsePath";
};

template <>
struct xml_io_stream_aggregate<sensor::FrequencyResponsePath> {
  static constexpr bool value = true;
};

// FrequencyRange format tags and XML I/O

template <>
struct format_tag_aggregate<sensor::FrequencyRange> {
  constexpr static bool value = true;
};

template <>
struct xml_io_stream_name<sensor::FrequencyRange> {
  static constexpr std::string_view name = "SensorFrequencyRange";
};

template <>
struct xml_io_stream_aggregate<sensor::FrequencyRange> {
  static constexpr bool value = true;
};

// HeterodyneFrequencyRange format tags and XML I/O

template <>
struct std::formatter<sensor::HeterodyneFrequencyRange>
    : format_tag_inherit<sensor::FrequencyRange,
                         sensor::HeterodyneFrequencyRange> {};

template <>
struct xml_io_stream<sensor::HeterodyneFrequencyRange>
    : xml_io_stream_inherit<sensor::FrequencyRange,
                            sensor::HeterodyneFrequencyRange> {};
