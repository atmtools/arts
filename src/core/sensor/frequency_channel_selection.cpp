#include "frequency_channel_selection.h"

#include <matpack.h>

#include <algorithm>
#include <boost/math/distributions/normal.hpp>
#include <numeric>
#include <ranges>
#include <vector>

#include "compare.h"

namespace sensor {
const AscendingGrid& Channel::freq_grid() const { return channel.grid<0>(); }

const Vector& Channel::weights() const { return channel.data; }

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

void Spectrometer::sync_frequency_grids() {
  std::vector<Numeric> f;
  f.reserve(std::transform_reduce(channels.begin(),
                                  channels.end(),
                                  Size{},
                                  std::plus<>{},
                                  [](const Channel& c) { return c.size(); }));

  for (const auto& channel : channels) {
    const auto& freqs = channel.freq_grid();
    f.insert(f.end(), freqs.begin(), freqs.end());
  }

  stdr::sort(f);
  f.erase(std::unique(f.begin(), f.end()), f.end());

  const auto fs = AscendingGrid(std::move(f));
  Vector ws(fs.size(), 0);

  for (Channel& channel : channels) {
    const auto& ch_fs = channel.freq_grid();
    const auto& ch_ws = channel.weights();
    const Size n      = ch_fs.size();

    const auto fs_ptr_first = stdr::lower_bound(fs, ch_fs.front());
    const Size OFFSET       = std::distance(fs.begin(), fs_ptr_first);
    const auto ws_ptr_first = ws.begin() + OFFSET;

    auto fs_ptr = fs_ptr_first;
    auto ws_ptr = ws_ptr_first;
    for (Size i = 0; i < n; i++) {
      while (*fs_ptr < ch_fs[i]) {
        ++fs_ptr;
        ++ws_ptr;
      }
      *ws_ptr = ch_ws[i];
    }

    channel.channel.grid<0>() = fs;
    channel.channel.data      = ws;

    fs_ptr = fs_ptr_first;
    ws_ptr = ws_ptr_first;
    for (Size i = 0; i < n; i++) {
      while (*fs_ptr < ch_fs[i]) {
        ++fs_ptr;
        ++ws_ptr;
      }
      *ws_ptr = 0.0;
    }
  }
}

bool Spectrometer::is_synced() const {
  return channels.empty() or stdr::all_of(channels,
                                          Cmp::eq(channels.front().freq_grid()),
                                          &Channel::freq_grid);
}

Spectrometer::Spectrometer(const Channel& base_channel,
                           const AscendingGrid& freq_offsets)
    : channels(freq_offsets.size(), base_channel) {
  for (Size i = 0; i < channels.size(); ++i) {
    Vector f                       = channels[i].channel.grid<0>();
    f                             += freq_offsets[i];
    channels[i].channel.grid<0>()  = f;
  }

  sync_frequency_grids();
}

Spectrometer::Spectrometer(const AscendingGrid& freq_offsets)
    : channels(std::from_range,
               freq_offsets |
                   stdv::transform([](auto f) -> DiracChannel { return f; })) {
  sync_frequency_grids();
}
}  // namespace sensor

void xml_io_stream<sensor::Spectrometer>::write(
    std::ostream& os,
    const sensor::Spectrometer& spec,
    bofstream* pbofs,
    std::string_view name) {
  xml_io_stream<std::vector<sensor::Channel>>::write(
      os, spec.channels, pbofs, name);
}

void xml_io_stream<sensor::Spectrometer>::read(std::istream& is,
                                               sensor::Spectrometer& spec,
                                               bifstream* pbifs) {
  xml_io_stream<std::vector<sensor::Channel>>::read(is, spec.channels, pbifs);
}
