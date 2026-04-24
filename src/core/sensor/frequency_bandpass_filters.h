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
                               const std::vector<Channel>& channels);
};
}  // namespace sensor
