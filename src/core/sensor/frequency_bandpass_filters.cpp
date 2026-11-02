#include "frequency_bandpass_filters.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <utility>

namespace sensor {
namespace {
constexpr Numeric response_eps = 64 * std::numeric_limits<Numeric>::epsilon();

Numeric scale(Numeric x) { return std::max<Numeric>(1.0, std::abs(x)); }

bool is_close(Numeric a, Numeric b) { return std::abs(a - b) <= response_eps * std::max(scale(a), scale(b)); }

//! The piece's response at f: linear between its frequencies, zero outside them
Numeric sample_piece(const SortedGriddedField1& piece, Numeric f) {
  const auto& grid = piece.grid<0>();

  if (grid.empty() or f < grid.front() or f > grid.back()) return 0.0;

  auto it = stdr::lower_bound(grid, f);

  if (it == grid.begin()) return piece.data.front();
  if (it == grid.end()) return piece.data.back();

  const Size upper = static_cast<Size>(std::distance(grid.begin(), it));
  if (is_close(*it, f)) return piece.data[upper];

  const Size    lower = upper - 1;
  const Numeric t     = (f - grid[lower]) / (grid[upper] - grid[lower]);

  return piece.data[lower] + t * (piece.data[upper] - piece.data[lower]);
}

//! The frequencies of grid strictly between low and high
void add_inner_points(std::vector<Numeric>& points, const AscendingGrid& grid, Numeric low, Numeric high) {
  for (auto it = stdr::upper_bound(grid, low); it != grid.end() and *it < high; ++it) points.push_back(*it);
}

/** The image of one response piece through one response path
 *
 *  A delta piece maps to a delta weighted by the path's response.  Any other
 *  piece is clipped to the path's local range, gets the frequencies of the
 *  path's filters as further points, is multiplied by the path's response, and
 *  is mapped to the global frequencies.  Images of the different paths are
 *  separate pieces, so the response of a heterodyne channel is the sum over
 *  its sidebands.
 */
std::optional<SortedGriddedField1> map_piece(const SortedGriddedField1& piece, const FrequencyResponsePath& path) {
  const auto& grid = piece.grid<0>();
  if (grid.empty()) return std::nullopt;

  if (grid.size() == 1) {
    const Numeric gain = path.local_weight(grid[0]);
    if (gain == 0.0) return std::nullopt;
    return SortedGriddedField1{.data_name  = piece.data_name,
                               .data       = Vector{piece.data[0] * gain},
                               .grid_names = piece.grid_names,
                               .grids      = {AscendingGrid{path.map_to_global(grid[0])}}};
  }

  const Numeric low  = std::max(grid.front(), path.local_range[0]);
  const Numeric high = std::min(grid.back(), path.local_range[1]);
  if (not(low < high) or is_close(low, high)) return std::nullopt;

  std::vector<Numeric> local{low, high};
  add_inner_points(local, grid, low, high);
  for (const auto& filter : path.filters) add_inner_points(local, filter.grid<0>(), low, high);
  stdr::sort(local);
  local.erase(std::unique(local.begin(), local.end(), &is_close), local.end());

  std::vector<std::pair<Numeric, Numeric>> samples;
  samples.reserve(local.size());
  for (const Numeric f : local) {
    samples.emplace_back(path.map_to_global(f), sample_piece(piece, f) * path.local_weight(f));
  }
  stdr::sort(samples, {}, &std::pair<Numeric, Numeric>::first);

  std::vector<Numeric> global(samples.size());
  Vector               response(samples.size());
  for (Size i = 0; i < samples.size(); i++) std::tie(global[i], response[i]) = samples[i];

  return SortedGriddedField1{.data_name  = piece.data_name,
                             .data       = std::move(response),
                             .grid_names = piece.grid_names,
                             .grids      = {AscendingGrid{global}}};
}
}  // namespace

std::vector<Channel> make_bandpass_channels(const FrequencyRange& range, const std::span<const Channel>& channels) {
  std::vector<Channel> out;
  out.reserve(channels.size());

  for (const auto& channel : channels) {
    Channel& mapped = out.emplace_back(Channel{.response = {}, .polarization = channel.polarization});

    for (const auto& piece : channel.response) {
      for (const auto& path : range.paths()) {
        if (auto image = map_piece(piece, path)) mapped.response.push_back(std::move(*image));
      }
    }
  }

  return out;
}

std::vector<Channel> make_bandpass_channels(const FrequencyRange& range, const Spectrometer& spectrometer) {
  return make_bandpass_channels(range, std::span<const Channel>{spectrometer.channels});
}
}  // namespace sensor
