#include "frequency_channel_selection.h"

#include <matpack.h>

#include <boost/math/distributions/normal.hpp>
#include <ranges>

namespace sensor {
const AscendingGrid& Channel::freq_grid() const { return channel.grid<0>(); }

const Vector& Channel::weights() const { return channel.data; }

bool Channel::is_always_relative() const { return freq_grid().front() <= 0; }

DiracChannel::DiracChannel(Numeric f)
    : Channel{.channel = {
                  .data_name  = "dirac"s,
                  .data       = Vector(1, 1.0),
                  .grid_names = std::array<String, 1>{"frequency"s},
                  .grids = std::array<AscendingGrid, 1>{AscendingGrid{f}}}} {}

DiracChannel::DiracChannel() : DiracChannel(0) {}

BoxChannel::BoxChannel(AscendingGrid f)
    : Channel{
          .channel = {
              .data_name = "box"s,
              .data = Vector(f.size(), 1.0 / static_cast<Numeric>(f.size())),
              .grid_names = std::array<String, 1>{"frequency"s},
              .grids      = std::array<AscendingGrid, 1>{std::move(f)}}} {}

BoxChannel::BoxChannel(Numeric lower, Numeric upper, Size N)
    : BoxChannel{nlinspace(lower, upper, N)} {}

BoxChannel::BoxChannel(Numeric hw, Size N) : BoxChannel(-hw, hw, N) {}

namespace {
SortedGriddedField1 gauss(AscendingGrid&& f, Numeric f0, Numeric std) {
  using gauss = boost::math::normal_distribution<Numeric>;

  Vector data{std::from_range,
              f | stdv::transform([dist = gauss(f0, std)](Numeric fi) {
                return pdf(dist, fi);
              })};

  return {.data_name  = "gaussian"s,
          .data       = std::move(data),
          .grid_names = std::array<String, 1>{"frequency"s},
          .grids      = std::array<AscendingGrid, 1>{std::move(f)}};
}
}  // namespace

GaussianChannel::GaussianChannel(AscendingGrid f, Numeric f0, Numeric std)
    : Channel{.channel = gauss(std::move(f), f0, std)} {}

GaussianChannel::GaussianChannel(Numeric f0, Numeric std, Size N, Size M)
    : GaussianChannel(
          nlinspace(
              -static_cast<Numeric>(M) * std, static_cast<Numeric>(M) * std, N),
          f0,
          std) {}

GaussianChannel::GaussianChannel(AscendingGrid f, Numeric std)
    : GaussianChannel(std::move(f), 0, std) {}

GaussianChannel::GaussianChannel(Numeric std, Size N, Size M)
    : GaussianChannel(
          nlinspace(
              -static_cast<Numeric>(M) * std, static_cast<Numeric>(M) * std, N),
          std) {}

static_assert(FrequencyChannelSelection<Channel>);
static_assert(FrequencyChannelSelection<BoxChannel>);
static_assert(FrequencyChannelSelection<DiracChannel>);
static_assert(FrequencyChannelSelection<GaussianChannel>);
}  // namespace sensor
