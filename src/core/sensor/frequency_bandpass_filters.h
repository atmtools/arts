#pragma once

#include <matpack.h>

#include <span>

#include "frequency_channel_selection.h"
#include "frequency_range_selection.h"

namespace sensor {
//! Real frequency bandpass filter.  Others inherit from this.
struct BandpassFilter;

//! Sets the bandpass filter from weights on derived frequeny ranges.
struct FrequencyRangeBandpassFilter;

//! Concept that creates a valid frequency bandpass filter for a given set of channels and frequency ranges.
template <typename T>
concept FrequencyBandpassFilter = std::derived_from<T, BandpassFilter>;

struct BandpassFilter {
  std::vector<SortedGriddedField1> filters;

  [[nodiscard]] Numeric operator()(Numeric f) const;
  [[nodiscard]] Vector operator()(ConstVectorView f) const;
};

struct FrequencyRangeBandpassFilter final : BandpassFilter {
  FrequencyRangeBandpassFilter(const FrequencyRange& range,
                               const std::span<const Channel>& channels);
  FrequencyRangeBandpassFilter(const FrequencyRange& range,
                               const Spectrometer& spectrometer);
};
}  // namespace sensor

// BandpassFilter format tags and XML I/O

template <>
struct format_tag_aggregate<sensor::BandpassFilter> {
  constexpr static bool value = true;
};

template <>
struct xml_io_stream_name<sensor::BandpassFilter> {
  static constexpr std::string_view name = "SensorBandpassFilter";
};

template <>
struct xml_io_stream_aggregate<sensor::BandpassFilter> {
  static constexpr bool value = true;
};

// FrequencyRangeBandpassFilter format tags and XML I/O

template <>
struct format_tag_aggregate<sensor::FrequencyRangeBandpassFilter> {
  constexpr static bool value = true;
};

template <>
struct xml_io_stream_name<sensor::FrequencyRangeBandpassFilter> {
  static constexpr std::string_view name = "SensorFrequencyRangeBandpassFilter";
};

template <>
struct xml_io_stream_aggregate<sensor::FrequencyRangeBandpassFilter> {
  static constexpr bool value = true;
};
