#include "frequency_bandpass_filters.h"

#include <lagrange_interp.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <limits>
#include <string_view>
#include <variant>

namespace sensor {
namespace {
constexpr Numeric response_eps = 64 * std::numeric_limits<Numeric>::epsilon();

Numeric scale(Numeric x) { return std::max<Numeric>(1.0, std::abs(x)); }

bool is_close(Numeric a, Numeric b) { return std::abs(a - b) <= response_eps * std::max(scale(a), scale(b)); }

Numeric sample_filter(const SortedGriddedField1& filter, Numeric f) {
  const auto& grid = filter.grid<0>();

  if (grid.empty() or f < grid.front() or f > grid.back()) return 0.0;

  auto it = stdr::lower_bound(grid, f);

  if (it == grid.begin()) return filter.data.front();
  if (it == grid.end()) return filter.data.back();

  const Size upper = static_cast<Size>(std::distance(grid.begin(), it));
  if (is_close(*it, f)) return filter.data[upper];

  const Size    lower = upper - 1;
  const Numeric x0    = grid[lower];
  const Numeric x1    = grid[upper];
  const Numeric y0    = filter.data[lower];
  const Numeric y1    = filter.data[upper];
  const Numeric t     = (f - x0) / (x1 - x0);

  return y0 + t * (y1 - y0);
}

void add_support_points(std::vector<Numeric>& points, const AscendingGrid& grid, Numeric low, Numeric high) {
  if (grid.empty() or low > high) return;

  points.push_back(low);
  if (not is_close(low, high)) points.push_back(high);

  auto lower = stdr::lower_bound(grid, low);
  auto upper = stdr::upper_bound(grid, high);
  for (auto it = lower; it != upper; ++it) points.push_back(*it);
}

void sort_unique(std::vector<Numeric>& points) {
  stdr::sort(points);
  points.erase(std::unique(points.begin(), points.end(), &is_close), points.end());
}

SortedGriddedField1 make_filter(const std::vector<std::pair<Numeric, Numeric>>& samples, std::string_view name) {
  std::vector<Numeric> grid(samples.size());
  Vector               data(samples.size());

  for (Size i = 0; i < samples.size(); i++) {
    grid[i] = samples[i].first;
    data[i] = samples[i].second;
  }

  return {.data_name  = String{name},
          .data       = std::move(data),
          .grid_names = std::array<String, 1>{"frequency"s},
          .grids      = std::array<AscendingGrid, 1>{AscendingGrid{grid}}};
}

void deduplicate_zero_frequency_images(std::vector<std::pair<Numeric, Numeric>>& samples) {
  std::vector<std::pair<Numeric, Numeric>> unique_samples;
  unique_samples.reserve(samples.size());

  for (const auto& sample : samples) {
    const bool duplicate =
        stdr::any_of(unique_samples, [&](const auto& kept) { return is_close(kept.first, sample.first); });

    if (not duplicate) unique_samples.push_back(sample);
  }

  samples = std::move(unique_samples);
}

using InterpolationPoint = std::variant<lagrange_interp::lag_t<0, lagrange_interp::grid_identity>,
                                        lagrange_interp::lag_t<1, lagrange_interp::grid_identity>>;

struct FoldedLocalPoint {
  InterpolationPoint                       interpolation;
  std::vector<std::pair<Numeric, Numeric>> folded_samples;
};

InterpolationPoint make_interpolation_point(const AscendingGrid& grid, Numeric f) {
  return lagrange_interp::variant_lag<lagrange_interp::grid_identity>(grid, f);
}

Numeric interpolate_channel_weight(const Vector& weights, const InterpolationPoint& point) {
  return std::visit([&weights](const auto& lag) { return lagrange_interp::interp(weights, lag); }, point);
}

std::vector<Numeric> collect_local_points(const FrequencyRange& range, const AscendingGrid& channel_grid) {
  std::vector<Numeric> local_points;

  for (const auto& path : range.paths()) {
    if (channel_grid.empty()) continue;

    const Numeric local_low  = std::max(path.local_range[0], channel_grid.front());
    const Numeric local_high = std::min(path.local_range[1], channel_grid.back());

    if (local_low > local_high and not is_close(local_low, local_high)) { continue; }

    add_support_points(local_points, channel_grid, local_low, local_high);
    for (const auto& filter : path.filters) {
      add_support_points(local_points, filter.grid<0>(), local_low, local_high);
    }
  }

  sort_unique(local_points);
  return local_points;
}

std::vector<std::pair<Numeric, Numeric>> fold_unit_samples(const FrequencyRange& range, Numeric local_frequency) {
  std::vector<std::pair<Numeric, Numeric>> folded_samples;
  folded_samples.reserve(range.size());

  for (const auto& path : range.paths()) {
    const Numeric path_weight = path.local_weight(local_frequency);
    if (path_weight == 0.0) continue;

    folded_samples.emplace_back(path.map_to_global(local_frequency), path_weight);
  }

  if (folded_samples.empty()) return folded_samples;

  const Numeric fold_count = static_cast<Numeric>(folded_samples.size());
  for (auto& sample : folded_samples) { sample.second /= fold_count; }

  if (is_close(local_frequency, 0.0)) { deduplicate_zero_frequency_images(folded_samples); }

  return folded_samples;
}

void combine_samples(std::vector<std::pair<Numeric, Numeric>>& samples) {
  std::sort(samples.begin(), samples.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

  std::vector<std::pair<Numeric, Numeric>> combined;
  combined.reserve(samples.size());
  for (const auto& sample : samples) {
    if (combined.empty() or not is_close(combined.back().first, sample.first)) {
      combined.push_back(sample);
    } else {
      combined.back().second += sample.second;
    }
  }

  samples = std::move(combined);
}

std::vector<Channel> build_channels(const FrequencyRange& range, const std::span<const Channel>& channels) {
  std::vector<Channel> out;
  out.reserve(channels.size());

  for (Size ichan = 0; ichan < channels.size(); ichan++) {
    const auto&                              channel      = channels[ichan];
    const auto&                              channel_grid = channel.freq_grid();
    std::vector<std::pair<Numeric, Numeric>> samples;
    std::vector<Numeric>                     local_points;

    for (const auto& path : range.paths()) {
      if (channel_grid.empty()) continue;

      const Numeric local_low  = std::max(path.local_range[0], channel_grid.front());
      const Numeric local_high = std::min(path.local_range[1], channel_grid.back());

      if (local_low > local_high and not is_close(local_low, local_high)) continue;

      add_support_points(local_points, channel_grid, local_low, local_high);
      for (const auto& filter : path.filters) {
        add_support_points(local_points, filter.grid<0>(), local_low, local_high);
      }
    }

    sort_unique(local_points);

    for (Numeric local_frequency : local_points) {
      const Numeric channel_weight = sample_filter(channel.channel, local_frequency);
      if (channel_weight == 0.0) continue;

      auto folded_samples = fold_unit_samples(range, local_frequency);
      if (folded_samples.empty()) continue;

      for (const auto& sample : folded_samples) { samples.emplace_back(sample.first, sample.second * channel_weight); }
    }

    combine_samples(samples);
    out.push_back(Channel{.channel = make_filter(samples, std::format("channel-response-{}", ichan))});
  }

  return out;
}

std::vector<Channel> build_synced_channels(const FrequencyRange& range, const Spectrometer& spectrometer) {
  const auto& channels = spectrometer.channels;
  if (channels.empty()) return {};

  if (not spectrometer.is_synced()) {
    return build_channels(range, std::span<const Channel>{channels.data(), channels.size()});
  }

  const auto& channel_grid = channels.front().freq_grid();
  const auto  local_points = collect_local_points(range, channel_grid);

  std::vector<FoldedLocalPoint> folded_points;
  folded_points.reserve(local_points.size());
  for (Numeric local_frequency : local_points) {
    auto folded_samples = fold_unit_samples(range, local_frequency);
    if (folded_samples.empty()) continue;

    folded_points.push_back({make_interpolation_point(channel_grid, local_frequency), std::move(folded_samples)});
  }

  std::vector<Channel> out;
  out.reserve(channels.size());

  for (Size ichan = 0; ichan < channels.size(); ++ichan) {
    const auto&                              weights = channels[ichan].weights();
    std::vector<std::pair<Numeric, Numeric>> samples;

    for (const auto& point : folded_points) {
      const Numeric channel_weight = interpolate_channel_weight(weights, point.interpolation);
      if (channel_weight == 0.0) continue;

      for (const auto& sample : point.folded_samples) {
        samples.emplace_back(sample.first, sample.second * channel_weight);
      }
    }

    combine_samples(samples);
    out.push_back(Channel{.channel = make_filter(samples, std::format("channel-response-{}", ichan))});
  }

  return out;
}
}  // namespace

std::vector<Channel> make_bandpass_channels(const FrequencyRange& range, const std::span<const Channel>& channels) {
  return build_channels(range, channels);
}

std::vector<Channel> make_bandpass_channels(const FrequencyRange& range, const Spectrometer& spectrometer) {
  return build_synced_channels(range, spectrometer);
}
}  // namespace sensor
