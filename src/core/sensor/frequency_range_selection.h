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

  std::vector<Vector2> global_ranges{{0, inf}};
  std::vector<Vector2> local_ranges{{0, inf}};

  [[nodiscard]] Size size() const;
  [[nodiscard]] const FrequencyResponsePath& path(Size index) const;
  [[nodiscard]] const std::vector<FrequencyResponsePath>& paths() const;

 protected:
  std::vector<FrequencyResponsePath> response_paths{{}};

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

template <>
struct std::formatter<sensor::FrequencyRange> {
  format_tags tags{};

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const sensor::FrequencyRange& v,
                              FmtContext& ctx) const {
    return tags.format(
        ctx, "GLOBALS: "sv, v.global_ranges, "; LOCALS : "sv, v.local_ranges);
  }
};

template <>
struct std::formatter<sensor::HeterodyneFrequencyRange> final
    : std::formatter<sensor::FrequencyRange> {};