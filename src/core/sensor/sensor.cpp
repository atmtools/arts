#include "sensor.h"

#include <debug.h>

#include <algorithm>
#include <array>
#include <format>
#include <numeric>
#include <ranges>
#include <utility>

bool SensorKey::operator==(const SensorKey& other) const {
  return other.sensor_elem == sensor_elem and other.type == type;
}

namespace sensor {
const String& Sensor::name() const {
  return std::visit([](const auto& gf) -> const String& { return gf.data_name; }, meta);
}

String& Sensor::name() {
  return std::visit([](auto& gf) -> String& { return gf.data_name; }, meta);
}

void Sensor::check() const {
  const std::array shape{static_cast<Index>(geometry.size()), static_cast<Index>(freq_grid.size())};

  for (Size i = 0; i < elements.size(); i++) {
    if (elements[i].shape() != shape)
      arts::user_error(
          R"(Sensor "{}": the weights of element {} must have one row per geometry and one column per frequency:
  weights shape: {:B,}
  geometries:    {}
  frequencies:   {}
)",
          name(),
          i,
          elements[i].shape(),
          geometry.size(),
          freq_grid.size());
  }

  if (count(meta) != elements.size())
    arts::user_error(
        R"(Sensor "{}": the meta describes {} elements but the sensor has {})", name(), count(meta), elements.size());
}

void Sensor::normalize(const Stokvec& pol) {
  for (auto& w : elements) sensor::normalize(w, pol);
}

Vector Sensor::extract(const StokvecMatrix& spectral_rad) const {
  if (spectral_rad.shape() != std::array{static_cast<Index>(geometry.size()), static_cast<Index>(freq_grid.size())})
    arts::user_error(R"(Sensor "{}": the spectral radiance must have one row per geometry and one column per frequency:
  spectral radiance shape: {:B,}
  geometries:              {}
  frequencies:             {}
)",
                     name(),
                     spectral_rad.shape(),
                     geometry.size(),
                     freq_grid.size());

  Vector out(elements.size(), 0.0);
  for (Size ie = 0; ie < elements.size(); ie++) {
    for (Size ip = 0; ip < geometry.size(); ip++) out[ie] += sumup(elements[ie], spectral_rad[ip], ip);
  }
  return out;
}

Size Sensor::flat_size(const SensorKeyType& key) const {
  switch (key) {
    using enum SensorKeyType;
    case freq: return freq_grid.size();
    case zen:  return geometry.size();
    case azi:  return geometry.size();
    case alt:  return geometry.size();
    case lat:  return geometry.size();
    case lon:  return geometry.size();
  }

  std::unreachable();
}

namespace {
//! The coordinate of a geometry that a sensor key selects
template <typename Self> decltype(auto) coordinate(Self&& poslos, const SensorKeyType& key) {
  switch (key) {
    using enum SensorKeyType;
    case zen:  return (std::forward<Self>(poslos).los[0]);
    case azi:  return (std::forward<Self>(poslos).los[1]);
    case alt:  return (std::forward<Self>(poslos).pos[0]);
    case lat:  return (std::forward<Self>(poslos).pos[1]);
    case lon:  return (std::forward<Self>(poslos).pos[2]);
    case freq: break;
  }

  std::unreachable();
}
}  // namespace

void Sensor::flat(VectorView x, const SensorKeyType& key) const {
  if (x.size() != flat_size(key))
    arts::user_error(R"(Sensor "{}": bad size. x.size(): {}, expected: {})", name(), x.size(), flat_size(key));

  if (key == SensorKeyType::freq) {
    x = freq_grid;
    return;
  }

  std::transform(
      geometry.begin(), geometry.end(), x.begin(), [&key](const PosLos& poslos) { return coordinate(poslos, key); });
}

Vector Sensor::flat(const SensorKeyType& key) const {
  Vector out(flat_size(key));
  flat(out, key);
  return out;
}

void Sensor::unflatten(const ConstVectorView& x, const SensorKeyType& key) {
  if (x.size() != flat_size(key))
    arts::user_error(R"(Sensor "{}": bad size. x.size(): {}, expected: {})", name(), x.size(), flat_size(key));

  if (key == SensorKeyType::freq) {
    freq_grid = AscendingGrid{Vector{x}};
    return;
  }

  for (Size i = 0; i < geometry.size(); i++) coordinate(geometry[i], key) = x[i];
}

GeometryIndex::GeometryIndex(const Sensor& sensor)
    : starts(sensor.geometry.size() + 1, 0), counts(sensor.elements.size(), 0) {
  const Size ngeom = sensor.geometry.size();
  const Size nfreq = sensor.freq_grid.size();

  //! Count the elements at each geometry: one for every run of an element's weights at one geometry
  for (Size ie = 0; ie < sensor.elements.size(); ie++) {
    const auto& w = sensor.elements[ie].vector();
    for (Size i = 0; i < w.size(); i++) {
      if (w[i].irow >= ngeom or w[i].icol >= nfreq)
        arts::user_error(R"(Sensor "{}": every weight must be at one of its {} geometries and {} frequencies)",
                         sensor.name(),
                         ngeom,
                         nfreq);

      if (i == 0 or w[i].irow != w[i - 1].irow) {
        starts[w[i].irow + 1]++;
        counts[ie]++;
      }
    }
  }

  std::partial_sum(starts.begin(), starts.end(), starts.begin());

  //! Place the runs, in increasing element order at each geometry
  entries.resize(starts.back());
  std::vector<Size> next(starts.begin(), starts.end() - 1);
  for (Size ie = 0; ie < sensor.elements.size(); ie++) {
    const std::span<const SparseStokvec> w{sensor.elements[ie].vector()};
    for (Size i = 0; i < w.size();) {
      Size n = 1;
      while (i + n < w.size() and w[i + n].irow == w[i].irow) n++;
      entries[next[w[i].irow]++] = {.element = ie, .weights = w.subspan(i, n)};
      i += n;
    }
  }
}

std::span<const GeometryWeights> GeometryIndex::operator[](Size ip) const {
  assert(ip + 1 < starts.size());
  return std::span{entries}.subspan(starts[ip], starts[ip + 1] - starts[ip]);
}

Size GeometryIndex::geometry_count(Size ie) const {
  assert(ie < counts.size());
  return counts[ie];
}
}  // namespace sensor

Size measurement_size(std::span<const Sensor> sensors) {
  Size n = 0;
  for (const auto& s : sensors) n += s.elements.size();
  return n;
}

Size measurement_start(std::span<const Sensor> sensors, Size isensor) {
  if (isensor > sensors.size())
    arts::user_error("Sensor index out of bounds: {} not in [0, {}]", isensor, sensors.size());

  return measurement_size(sensors.first(isensor));
}

std::pair<Size, Size> measurement_element(std::span<const Sensor> sensors, Size iy) {
  Size start = 0;
  for (Size i = 0; i < sensors.size(); i++) {
    const Size n = sensors[i].elements.size();
    if (iy < start + n) return {i, iy - start};
    start += n;
  }

  arts::user_error("Measurement element out of bounds: {} not in [0, {})", iy, start);
}

String cost_report(std::span<const Sensor> sensors) {
  constexpr std::string_view row = "{:>6}  {:<20}  {:>9}  {:<16}  {:>12}  {:>11}  {:>15}  {:>10}\n";

  String out = std::format(
      row, "Sensor", "Name", "Elements", "Meta shape", "Agenda calls", "Frequencies", "Spectral points", "Weights");

  Size total_elements = 0, total_calls = 0, total_points = 0, total_weights = 0;
  for (Size i = 0; i < sensors.size(); i++) {
    const Sensor& s = sensors[i];

    const Size calls   = s.geometry.size();
    const Size points  = calls * s.freq_grid.size();
    Size       weights = 0;
    for (const auto& w : s.elements) weights += w.size();

    const String shape = std::visit([](const auto& gf) { return std::format("{:B,}", gf.shape()); }, s.meta);

    out += std::format(row, i, s.name(), s.elements.size(), shape, calls, s.freq_grid.size(), points, weights);

    total_elements += s.elements.size();
    total_calls    += calls;
    total_points   += points;
    total_weights  += weights;
  }

  out += std::format(row, "Total", "", total_elements, "", total_calls, "", total_points, total_weights);

  // Geometries that more than one sensor simulates: each such sensor makes its own agenda call
  std::vector<std::pair<std::array<Numeric, 5>, Size>> geometries;
  for (Size i = 0; i < sensors.size(); i++) {
    for (const auto& g : sensors[i].geometry) {
      geometries.push_back({{g.pos[0], g.pos[1], g.pos[2], g.los[0], g.los[1]}, i});
    }
  }
  stdr::sort(geometries);
  geometries.erase(std::unique(geometries.begin(), geometries.end()), geometries.end());

  Size shared = 0, extra_calls = 0;
  for (auto it = geometries.begin(); it != geometries.end();) {
    const auto last = std::find_if(it, geometries.end(), [&it](const auto& g) { return g.first != it->first; });
    if (const auto n = static_cast<Size>(std::distance(it, last)); n > 1) {
      shared++;
      extra_calls += n - 1;
    }
    it = last;
  }

  if (shared != 0)
    out += std::format(
        "\nGeometries simulated by more than one sensor: {}, costing {} agenda calls more than simulating each once.\n",
        shared,
        extra_calls);

  return out;
}
