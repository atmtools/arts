#include "frequency_channel_selection.h"

#include <math_funcs.h>
#include <matpack.h>

#include <algorithm>
#include <boost/math/distributions/normal.hpp>
#include <limits>
#include <ranges>
#include <vector>

namespace sensor {
namespace {
constexpr Numeric frequency_eps = 64 * std::numeric_limits<Numeric>::epsilon();

bool is_close(Numeric a, Numeric b) {
  return std::abs(a - b) <= frequency_eps * std::max({Numeric{1}, std::abs(a), std::abs(b)});
}

SortedGriddedField1 make_piece(String name, AscendingGrid f, Vector response) {
  return {
      .data_name = std::move(name), .data = std::move(response), .grid_names = {"frequency"s}, .grids = {std::move(f)}};
}
}  // namespace

SortedGriddedField1 Channel::quadrature() const {
  std::vector<std::pair<Numeric, Numeric>> samples;

  Vector weights;
  for (const auto& piece : response) {
    if (not piece.ok()) arts::user_error("A channel response piece must have as many values as frequencies");

    const auto& f = piece.grid<0>();
    if (f.size() == 1) {
      samples.emplace_back(f[0], piece.data[0]);
      continue;
    }

    calculate_int_weights_arbitrary_grid(weights, f);
    for (Size i = 0; i < f.size(); i++) samples.emplace_back(f[i], weights[i] * piece.data[i]);
  }

  stdr::sort(samples, {}, &std::pair<Numeric, Numeric>::first);

  std::vector<Numeric> f;
  std::vector<Numeric> w;
  for (const auto& [fi, wi] : samples) {
    if (not f.empty() and is_close(f.back(), fi)) {
      w.back() += wi;
    } else {
      f.push_back(fi);
      w.push_back(wi);
    }
  }

  return make_piece("quadrature", AscendingGrid{std::move(f)}, Vector{std::move(w)});
}

DiracChannel::DiracChannel(Numeric f) : Channel{.response = {make_piece("dirac", AscendingGrid{f}, Vector{1.0})}} {}

DiracChannel::DiracChannel() : DiracChannel(0) {}

BoxChannel::BoxChannel(AscendingGrid f) : Channel{} {
  Vector one(f.size(), 1.0);
  response.push_back(make_piece("box", std::move(f), std::move(one)));
}

BoxChannel::BoxChannel(Numeric lower, Numeric upper, Size N) : BoxChannel{nlinspace(lower, upper, N)} {}

BoxChannel::BoxChannel(Numeric hw, Size N) : BoxChannel(-hw, hw, N) {}

namespace {
SortedGriddedField1 gauss(AscendingGrid&& f, Numeric f0, Numeric std) {
  using gauss = boost::math::normal_distribution<Numeric>;

  Vector data{std::from_range, f | stdv::transform([dist = gauss(f0, std)](Numeric fi) { return pdf(dist, fi); })};

  return make_piece("gaussian", std::move(f), std::move(data));
}
}  // namespace

GaussianChannel::GaussianChannel(AscendingGrid f, Numeric f0, Numeric std)
    : Channel{.response = {gauss(std::move(f), f0, std)}} {}

GaussianChannel::GaussianChannel(Numeric f0, Numeric std, Size N, Size M)
    : GaussianChannel(nlinspace(-static_cast<Numeric>(M) * std, static_cast<Numeric>(M) * std, N), f0, std) {}

GaussianChannel::GaussianChannel(AscendingGrid f, Numeric std) : GaussianChannel(std::move(f), 0, std) {}

GaussianChannel::GaussianChannel(Numeric std, Size N, Size M)
    : GaussianChannel(nlinspace(-static_cast<Numeric>(M) * std, static_cast<Numeric>(M) * std, N), std) {}

static_assert(FrequencyChannelSelection<Channel>);
static_assert(FrequencyChannelSelection<BoxChannel>);
static_assert(FrequencyChannelSelection<DiracChannel>);
static_assert(FrequencyChannelSelection<GaussianChannel>);

Spectrometer::Spectrometer(const Channel& base_channel, const AscendingGrid& freq_offsets)
    : channels(freq_offsets.size(), base_channel) {
  for (Size i = 0; i < channels.size(); ++i) {
    for (auto& piece : channels[i].response) {
      Vector f         = piece.grid<0>();
      f               += freq_offsets[i];
      piece.grid<0>()  = AscendingGrid{std::move(f)};
    }
  }
}

Spectrometer::Spectrometer(const AscendingGrid& freq_offsets)
    : channels(std::from_range, freq_offsets | stdv::transform([](auto f) -> Channel { return DiracChannel{f}; })) {}
}  // namespace sensor

void xml_io_stream<sensor::Spectrometer>::write(std::ostream&               os,
                                                const sensor::Spectrometer& spec,
                                                bofstream*                  pbofs,
                                                std::string_view            name) {
  xml_io_stream<std::vector<sensor::Channel>>::write(os, spec.channels, pbofs, name);
}

void xml_io_stream<sensor::Spectrometer>::read(std::istream& is, sensor::Spectrometer& spec, bifstream* pbifs) {
  xml_io_stream<std::vector<sensor::Channel>>::read(is, spec.channels, pbifs);
}
