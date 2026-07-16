#pragma once

#include <matpack.h>

#include <span>

#include "frequency_channel_selection.h"
#include "frequency_range_selection.h"

namespace sensor {
[[nodiscard]] std::vector<Channel> make_bandpass_channels(const FrequencyRange&           range,
                                                          const std::span<const Channel>& channels);

[[nodiscard]] std::vector<Channel> make_bandpass_channels(const FrequencyRange& range,
                                                          const Spectrometer&   spectrometer);
}  // namespace sensor
