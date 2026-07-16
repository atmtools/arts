#include "frequency_range_selection.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <span>

#include "matpack_mdspan_cdata_t.h"

namespace sensor {
Size FrequencyRange::size() const { return response_paths.size(); }

const FrequencyResponsePath& FrequencyRange::path(Size index) const { return response_paths.at(index); }

const std::vector<FrequencyResponsePath>& FrequencyRange::paths() const { return response_paths; }

void FrequencyRange::sync_ranges() {
  global_ranges.clear();
  local_ranges.clear();

  global_ranges.reserve(response_paths.size());
  local_ranges.reserve(response_paths.size());

  for (const auto& path : response_paths) {
    global_ranges.push_back(path.global_range());
    local_ranges.push_back(path.local_range);
  }
}

namespace {
constexpr Numeric response_eps = 64 * std::numeric_limits<Numeric>::epsilon();

Numeric scale(Numeric x) { return std::max<Numeric>(1.0, std::abs(x)); }

bool is_close(Numeric a, Numeric b) {
  if (std::isinf(a) or std::isinf(b)) return a == b;
  return std::abs(a - b) <= response_eps * std::max(scale(a), scale(b));
}

bool is_in_closed_interval(Numeric x, Numeric low, Numeric high) {
  return x >= low - response_eps * scale(low) and x <= high + response_eps * scale(high);
}

void assert_frange(const Vector2& sideband) {
  if (sideband[0] > sideband[1] or is_close(sideband[0], sideband[1]) or sideband[0] < 0) {
    throw std::invalid_argument(std::format(
        "Frequency range must be unique, sorted, and non-negative.  Got: [{}, {}].", sideband[0], sideband[1]));
  }
}

void assert_bandpass_ranges(const std::span<const Vector2>& bandpasses) {
  if (bandpasses.empty()) { throw std::invalid_argument("At least one sideband must be provided."); }

  for (auto sideband : bandpasses) assert_frange(sideband);
}

void assert_lo_nonnegative(const Numeric& LO) {
  if (LO < 0) { throw std::invalid_argument(std::format("LO frequency must be non-negative.  Got: {}.", LO)); }
}

void assert_clocks(const std::span<const Numeric>& clocks) {
  if (clocks.empty()) { throw std::invalid_argument("At least one LO frequency must be provided."); }

  for (auto clock : clocks) assert_lo_nonnegative(clock);
}

void assert_weighted_bandpass(const SortedGriddedField1& bandpass_filter) {
  const auto& grid = bandpass_filter.grid<0>();

  if (grid.empty()) { throw std::invalid_argument("Bandpass filter grid must not be empty."); }

  if (grid.size() != bandpass_filter.data.size()) {
    throw std::invalid_argument(
        std::format("Bandpass filter grid and data must have the same size.  Got {} grid points and {} weights.",
                    grid.size(),
                    bandpass_filter.data.size()));
  }

  for (Size i = 1; i < grid.size(); i++) {
    if (grid[i - 1] > grid[i] or is_close(grid[i - 1], grid[i])) {
      throw std::invalid_argument("Bandpass filter grid must be unique and sorted.");
    }
  }
}

Numeric sample_filter(const SortedGriddedField1& filter, Numeric f) {
  const auto& grid = filter.grid<0>();

  if (grid.empty() or f < grid.front() or f > grid.back()) return 0.0;

  auto it = std::lower_bound(grid.begin(), grid.end(), f);

  if (it == grid.begin()) return filter.data.front();
  if (it == grid.end()) return filter.data.back();

  const Size upper = static_cast<Size>(std::distance(grid.begin(), it));
  if (is_close(*it, f)) return filter.data[upper];

  const Size    lower = upper - 1;
  const Numeric x0    = grid[lower];
  const Numeric x1    = grid[upper];
  const Numeric y0    = filter.data[lower];
  const Numeric y1    = filter.data[upper];

  const Numeric t = (f - x0) / (x1 - x0);
  return y0 + t * (y1 - y0);
}

std::vector<SortedGriddedField1> transformed_filters(const std::vector<SortedGriddedField1>& filters,
                                                     Numeric                                 LO,
                                                     bool                                    upper) {
  std::vector<SortedGriddedField1> out;
  out.reserve(filters.size());

  for (const auto& filter : filters) {
    const auto&          grid = filter.grid<0>();
    std::vector<Numeric> mapped(grid.size());
    Vector               weights(filter.data.size());

    if (upper) {
      for (Size i = 0; i < grid.size(); i++) {
        mapped[i]  = grid[i] - LO;
        weights[i] = filter.data[i];
      }
    } else {
      for (Size i = 0; i < grid.size(); i++) {
        const Size j = grid.size() - 1 - i;
        mapped[i]    = LO - grid[j];
        weights[i]   = filter.data[j];
      }
    }

    out.push_back({.data_name  = filter.data_name,
                   .data       = std::move(weights),
                   .grid_names = filter.grid_names,
                   .grids      = std::array<AscendingGrid, 1>{AscendingGrid{mapped}}});
  }

  return out;
}

void apply_interval_clip(std::vector<FrequencyResponsePath>& paths, Numeric low, Numeric high) {
  std::erase_if(paths, [low, high](FrequencyResponsePath& path) {
    path.local_range[0] = std::max(path.local_range[0], low);
    path.local_range[1] = std::min(path.local_range[1], high);
    return path.local_range[0] >= path.local_range[1] or is_close(path.local_range[0], path.local_range[1]);
  });
}
}  // namespace

Vector2 FrequencyResponsePath::global_range() const {
  return {map_to_global(local_range[0]), map_to_global(local_range[1])};
}

Numeric FrequencyResponsePath::map_to_global(Numeric local_frequency) const {
  return intercept + slope * local_frequency;
}

std::optional<Numeric> FrequencyResponsePath::map_to_local(Numeric global_frequency) const {
  const Numeric local_frequency = (global_frequency - intercept) / slope;

  if (not is_in_closed_interval(local_frequency, local_range[0], local_range[1])) { return std::nullopt; }

  return std::clamp(local_frequency, local_range[0], local_range[1]);
}

Numeric FrequencyResponsePath::local_weight(Numeric local_frequency) const {
  if (not is_in_closed_interval(local_frequency, local_range[0], local_range[1])) { return 0.0; }

  Numeric weight = 1.0;

  for (const auto& filter : filters) {
    weight *= sample_filter(filter, local_frequency);
    if (weight == 0.0) break;
  }

  return weight;
}

Numeric FrequencyResponsePath::global_weight(Numeric global_frequency) const {
  const auto local_frequency = map_to_local(global_frequency);
  return local_frequency.has_value() ? local_weight(*local_frequency) : 0.0;
}

HeterodyneFrequencyRange::HeterodyneFrequencyRange(const std::span<const Numeric>& clocks,
                                                   const std::span<const Vector2>& bandpasses)
    : FrequencyRange() {
  const Size N = clocks.size();

  if (N != bandpasses.size()) {
    throw std::invalid_argument(
        std::format("Number of clock frequencies and bandpasses must match.  Got: {} clock "
                    "frequencies and {} bandpasses.",
                    clocks.size(),
                    bandpasses.size()));
  }

  assert_bandpass_ranges(bandpasses);
  assert_clocks(clocks);

  for (Size i = 0; i < N; i++) {
    apply_bandpass(bandpasses[i]);
    apply_mixer(clocks[i]);
  }
}

HeterodyneFrequencyRange::HeterodyneFrequencyRange(Numeric clock_frequency, const Vector2& sideband)
    : HeterodyneFrequencyRange(std::span{&clock_frequency, 1}, std::span{&sideband, 1}) {}

void HeterodyneFrequencyRange::apply_lowpass(Numeric upper_frequency) {
  if (upper_frequency < 0) {
    throw std::invalid_argument(std::format("Lowpass cutoff must be non-negative.  Got: {}.", upper_frequency));
  }

  apply_interval_clip(response_paths, 0.0, upper_frequency);
  sync_ranges();
}

void HeterodyneFrequencyRange::apply_highpass(Numeric lower_frequency) {
  if (lower_frequency < 0) {
    throw std::invalid_argument(std::format("Highpass cutoff must be non-negative.  Got: {}.", lower_frequency));
  }

  apply_interval_clip(response_paths, lower_frequency, inf);
  sync_ranges();
}

void HeterodyneFrequencyRange::apply_bandpass(const Vector2& bandpass_range) {
  assert_frange(bandpass_range);
  apply_interval_clip(response_paths, bandpass_range[0], bandpass_range[1]);
  sync_ranges();
}

void HeterodyneFrequencyRange::apply_bandpass(const SortedGriddedField1& bandpass_filter) {
  assert_weighted_bandpass(bandpass_filter);

  apply_interval_clip(response_paths, bandpass_filter.grid<0>().front(), bandpass_filter.grid<0>().back());

  for (auto& path : response_paths) path.filters.push_back(bandpass_filter);

  sync_ranges();
}

void HeterodyneFrequencyRange::apply_mixer(Numeric clock_frequency) {
  assert_lo_nonnegative(clock_frequency);

  std::vector<FrequencyResponsePath> mixed_paths;
  mixed_paths.reserve(2 * response_paths.size());

  for (const auto& path : response_paths) {
    const Numeric upper_low  = std::max(path.local_range[0], clock_frequency);
    const Numeric upper_high = path.local_range[1];

    if (upper_low < upper_high and not is_close(upper_low, upper_high)) {
      auto upper         = path;
      upper.intercept   += upper.slope * clock_frequency;
      upper.local_range  = {upper_low - clock_frequency, upper_high - clock_frequency};
      upper.filters      = transformed_filters(path.filters, clock_frequency, true);
      mixed_paths.push_back(std::move(upper));
    }

    const Numeric lower_low  = path.local_range[0];
    const Numeric lower_high = std::min(path.local_range[1], clock_frequency);

    if (lower_low < lower_high and not is_close(lower_low, lower_high)) {
      auto lower         = path;
      lower.intercept   += lower.slope * clock_frequency;
      lower.slope       *= -1;
      lower.local_range  = {clock_frequency - lower_high, clock_frequency - lower_low};
      lower.filters      = transformed_filters(path.filters, clock_frequency, false);
      mixed_paths.push_back(std::move(lower));
    }
  }

  response_paths = std::move(mixed_paths);
  sync_ranges();
}

Vector HeterodyneFrequencyRange::local_response(ConstVectorView local_frequency_grid, Size path_index) const {
  const auto& response_path = path(path_index);
  Vector      response(local_frequency_grid.size());

  for (Size i = 0; i < response.size(); i++) { response[i] = response_path.local_weight(local_frequency_grid[i]); }

  return response;
}

Vector HeterodyneFrequencyRange::global_response(ConstVectorView global_frequency_grid, Size path_index) const {
  const auto& response_path = path(path_index);
  Vector      response(global_frequency_grid.size());

  for (Size i = 0; i < response.size(); i++) { response[i] = response_path.global_weight(global_frequency_grid[i]); }

  return response;
}

static_assert(FrequencyRangeSelection<FrequencyRange>);
static_assert(FrequencyRangeSelection<HeterodyneFrequencyRange>);
}  // namespace sensor
