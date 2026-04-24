#pragma once

#include <matpack.h>

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
  [[nodiscard]] bool is_always_relative() const;
};

struct BoxChannel final : Channel {
  BoxChannel(Numeric lower, Numeric upper, Size N);  // [lower, upper]
  BoxChannel(Numeric hw, Size N);                    // [-hw, hw]
  BoxChannel(AscendingGrid f);                       // f
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
};
}  // namespace sensor
