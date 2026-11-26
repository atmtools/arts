#include "subsurf_spectral_field.h"

#include <debug.h>
#include <field_concepts.h>
#include <lagrange_interp.h>

#include <cmath>
#include <optional>
#include <string_view>
#include <utility>

namespace Subsurface {
SpectralData::SpectralData()                                    = default;
SpectralData::SpectralData(const SpectralData &)                = default;
SpectralData::SpectralData(SpectralData &&) noexcept            = default;
SpectralData &SpectralData::operator=(const SpectralData &)     = default;
SpectralData &SpectralData::operator=(SpectralData &&) noexcept = default;

void SpectralData::adjust_interpolation_extrapolation() {
  using enum InterpolationExtrapolation;

  if (const auto *field = std::get_if<SpectralGeodeticField4>(&data)) {
    if (field->grid<0>().size() == 1) frq_upp = frq_low = Nearest;
    if (field->grid<1>().size() == 1) alt_upp = alt_low = Nearest;
    if (field->grid<2>().size() == 1) lat_upp = lat_low = Nearest;
    if (field->grid<3>().size() == 1) lon_upp = lon_low = Nearest;
  } else {
    frq_upp = frq_low = alt_upp = alt_low = lat_upp = lat_low = lon_upp = lon_low = Nearest;
  }
}

SpectralData::SpectralData(Numeric x) : data(x) { adjust_interpolation_extrapolation(); }

SpectralData &SpectralData::operator=(Numeric x) {
  data = x;
  adjust_interpolation_extrapolation();
  return *this;
}

SpectralData::SpectralData(SpectralGeodeticField4 x) : data(std::move(x)) { adjust_interpolation_extrapolation(); }

SpectralData &SpectralData::operator=(SpectralGeodeticField4 x) {
  data = std::move(x);
  adjust_interpolation_extrapolation();
  return *this;
}

SpectralData::SpectralData(SpectralFunctionalData x) : data(std::move(x)) { adjust_interpolation_extrapolation(); }

SpectralData &SpectralData::operator=(SpectralFunctionalData x) {
  data = std::move(x);
  adjust_interpolation_extrapolation();
  return *this;
}

String SpectralData::data_type() const {
  if (std::holds_alternative<SpectralGeodeticField4>(data)) return "SpectralGeodeticField4";
  if (std::holds_alternative<Numeric>(data)) return "Numeric";
  if (std::holds_alternative<SpectralFunctionalData>(data)) return "SpectralFunctionalData";

  arts::user_error("Cannot understand data type; is this a new type");
}

namespace {
using lagrange_interp::grid_identity;
using lagrange_interp::loncross;

//! Linear interpolation, or constant on a grid of one point
template <typename transform> using lags =
    std::variant<lagrange_interp::lag_t<0, transform>, lagrange_interp::lag_t<1, transform>>;

template <typename transform> lags<transform> lag(const auto &grid, Numeric x) {
  if (grid.empty()) arts::user_error("Cannot interpolate on an empty grid of the subsurface spectral data");
  if (grid.size() == 1) return grid.template lag<0, transform>(x);
  return grid.template lag<1, transform>(x);
}

const SpectralGeodeticField4 &checked(const SpectralGeodeticField4 &field) {
  if (not field.ok())
    arts::user_error("The subsurface spectral data shape {:B,} does not match its grid sizes {:B,}",
                     field.data.shape(),
                     field.shape());
  return field;
}

//! The coordinate to interpolate at after extrapolating outside the grid, or nothing where the data is zero
std::optional<Numeric> extrapolate(const auto                &grid,
                                   const Numeric              x,
                                   InterpolationExtrapolation low,
                                   InterpolationExtrapolation upp,
                                   std::string_view           coordinate) {
  using enum InterpolationExtrapolation;

  const bool below = x < grid.front();
  if (not below and x <= grid.back()) return x;

  switch (below ? low : upp) {
    case None:
      arts::user_error(
          "The {} {} is outside the grid [{}, {}] of the subsurface spectral data, which does not extrapolate there",
          coordinate,
          x,
          grid.front(),
          grid.back());
    case Zero:    return std::nullopt;
    case Nearest: return below ? grid.front() : grid.back();
    case Linear:  return x;
  }

  std::unreachable();
}

std::optional<lags<grid_identity>> frequency(const SpectralData           &d,
                                             const SpectralGeodeticField4 &field,
                                             const Numeric                 frq) {
  const auto f = extrapolate(field.grid<0>(), frq, d.frq_low, d.frq_upp, "frequency"sv);
  if (not f) return std::nullopt;
  return lag<grid_identity>(field.grid<0>(), *f);
}

struct spatial_lags {
  lags<grid_identity> alt;
  lags<grid_identity> lat;
  lags<loncross>      lon;
};

//! All coordinates are extrapolated first, so that None throws even where another coordinate gives zero
std::optional<spatial_lags> spatial(const SpectralData           &d,
                                    const SpectralGeodeticField4 &field,
                                    const Numeric                 alt,
                                    const Numeric                 lat,
                                    const Numeric                 lon) {
  const auto a = extrapolate(field.grid<1>(), alt, d.alt_low, d.alt_upp, "altitude"sv);
  const auto b = extrapolate(field.grid<2>(), lat, d.lat_low, d.lat_upp, "latitude"sv);
  const auto c = extrapolate(field.grid<3>(), lon, d.lon_low, d.lon_upp, "longitude"sv);
  if (not a or not b or not c) return std::nullopt;
  return spatial_lags{.alt = lag<grid_identity>(field.grid<1>(), *a),
                      .lat = lag<grid_identity>(field.grid<2>(), *b),
                      .lon = lag<loncross>(field.grid<3>(), *c)};
}

//! The spatially interpolated data at the frequency index k
Numeric spatial_value(const SpectralGeodeticField4 &field, const Index k, const spatial_lags &s, const bool log) {
  return std::visit(
      [data = field.data[k], log](const auto &a, const auto &b, const auto &c) {
        return log ? lagrange_interp::interp<lagrange_interp::log_transform>(data, a, b, c)
                   : lagrange_interp::interp(data, a, b, c);
      },
      s.alt,
      s.lat,
      s.lon);
}

//! The frequency interpolation of the spatially interpolated data, value(k) at frequency index k
Numeric frequency_value(const lags<grid_identity> &f, const auto &value, const bool log) {
  return std::visit(
      [&value, log](const auto &l) {
        Numeric out = 0.0;
        for (Size j = 0; j < l.size(); j++) {
          // The term is its own expression, so that it is not contracted with the sum
          const Numeric term  = l.data[j] * (log ? std::log(value(l.indx[j])) : value(l.indx[j]));
          out                += term;
        }
        return log ? std::exp(out) : out;
      },
      f);
}

//! The derivative of the log-interpolated value v = exp(sum w ln x) by a stored x is w v / x
void log_scale(std::vector<std::pair<Index, Numeric>> &weights, const ConstVectorView x) {
  Numeric ln_v = 0.0;
  for (const auto &[j, w] : weights) ln_v += w * std::log(x[j]);
  const Numeric v = std::exp(ln_v);
  for (auto &[j, w] : weights) w *= v / x[j];
}

//! The number of stored values per frequency
Index spatial_size(const SpectralGeodeticField4 &field) {
  return static_cast<Index>(field.grid<1>().size() * field.grid<2>().size() * field.grid<3>().size());
}
}  // namespace

Numeric SpectralData::at(const Numeric frq, const Numeric alt, const Numeric lat, const Numeric lon) const {
  return std::visit(
      [&]<typename T>(const T &x) -> Numeric {
        if constexpr (std::same_as<T, SpectralGeodeticField4>) {
          const SpectralGeodeticField4 &field = checked(x);
          const auto                    f     = frequency(*this, field, frq);
          const auto                    s     = spatial(*this, field, alt, lat, lon);
          if (not f or not s) return 0.0;
          return frequency_value(
              *f, [&](Index k) { return spatial_value(field, k, *s, log_interpolation); }, log_interpolation);
        } else if constexpr (std::same_as<T, Numeric>) {
          return x;
        } else {
          return x(frq, alt, lat, lon);
        }
      },
      data);
}

Vector SpectralData::at(const AscendingGrid &freq_grid, const Numeric alt, const Numeric lat, const Numeric lon) const {
  Vector out(freq_grid.size());

  std::visit(
      [&]<typename T>(const T &x) {
        if constexpr (std::same_as<T, SpectralGeodeticField4>) {
          const SpectralGeodeticField4 &field = checked(x);
          const auto                    s     = spatial(*this, field, alt, lat, lon);

          // The data at each of its frequencies, interpolated to the position once
          Vector values;
          if (s) {
            values.resize(field.grid<0>().size());
            for (Size k = 0; k < values.size(); k++) values[k] = spatial_value(field, k, *s, log_interpolation);
          }

          for (Size i = 0; i < freq_grid.size(); i++) {
            const auto f = frequency(*this, field, freq_grid[i]);
            out[i] = (f and s) ? frequency_value(*f, [&values](Index k) { return values[k]; }, log_interpolation) : 0.0;
          }
        } else if constexpr (std::same_as<T, Numeric>) {
          out = x;
        } else {
          for (Size i = 0; i < freq_grid.size(); i++) out[i] = x(freq_grid[i], alt, lat, lon);
        }
      },
      data);

  return out;
}

void SpectralData::flat_weight(std::vector<std::pair<Index, Numeric>> &weights,
                               const Numeric                           frq,
                               const Numeric                           alt,
                               const Numeric                           lat,
                               const Numeric                           lon) const {
  weights.clear();

  std::visit(
      [&]<typename T>(const T &x) {
        if constexpr (std::same_as<T, SpectralGeodeticField4>) {
          const SpectralGeodeticField4 &field = checked(x);
          const auto                    f     = frequency(*this, field, frq);
          const auto                    s     = spatial(*this, field, alt, lat, lon);
          if (not f or not s) return;

          const Index n    = spatial_size(field);
          const Index nlat = field.grid<2>().size();
          const Index nlon = field.grid<3>().size();

          std::visit(
              [&](const auto &l, const auto &a, const auto &b, const auto &c) {
                const auto w = lagrange_interp::interpweights(a, b, c);
                for (Size m = 0; m < l.size(); m++) {
                  if (l.data[m] == 0.0) continue;
                  for (Size i = 0; i < a.size(); i++) {
                    for (Size j = 0; j < b.size(); j++) {
                      for (Size k = 0; k < c.size(); k++) {
                        if (w[i, j, k] == 0.0) continue;
                        weights.emplace_back(l.indx[m] * n + (a.indx[i] * nlat + b.indx[j]) * nlon + c.indx[k],
                                             l.data[m] * w[i, j, k]);
                      }
                    }
                  }
                }
              },
              *f,
              s->alt,
              s->lat,
              s->lon);

          if (log_interpolation) log_scale(weights, field.data.view_as(field.data.size()));
        } else if constexpr (std::same_as<T, Numeric>) {
          weights.emplace_back(0, 1.0);
        }
      },
      data);
}

ConstVectorView SpectralData::flat_view() const {
  return std::visit(
      []<typename T>(const T &x) -> ConstVectorView {
        if constexpr (std::same_as<T, SpectralGeodeticField4>)
          return x.data.view_as(x.data.size());
        else if constexpr (std::same_as<T, Numeric>)
          return ConstVectorView{x};
        else
          return ConstVectorView{};
      },
      data);
}

VectorView SpectralData::flat_view() {
  return std::visit(
      []<typename T>(T &x) -> VectorView {
        if constexpr (std::same_as<T, SpectralGeodeticField4>)
          return x.data.view_as(x.data.size());
        else if constexpr (std::same_as<T, Numeric>)
          return VectorView{x};
        else
          return VectorView{};
      },
      data);
}

bool SpectralData::ok() const {
  if (const auto *v = std::get_if<SpectralGeodeticField4>(&data)) {
    return v->ok() and loncross::cycle(v->grid<3>().front()) == v->grid<3>().front() and
           loncross::cycle(v->grid<3>().back()) == v->grid<3>().back();
  }

  return true;
}

SpectralData &SpectralField::operator[](const SubsurfaceSpectralPropertyTag &key) { return props[key]; }

const SpectralData &SpectralField::operator[](const SubsurfaceSpectralPropertyTag &key) const {
  const auto ptr = props.find(key);
  if (ptr == props.end()) arts::user_error("The subsurface spectral field has no {}", key);
  return ptr->second;
}

bool SpectralField::contains(const SubsurfaceSpectralPropertyTag &key) const { return props.contains(key); }

std::vector<SubsurfaceSpectralPropertyTag> SpectralField::keys() const {
  std::vector<SubsurfaceSpectralPropertyTag> out;
  out.reserve(props.size());
  for (const auto &key : props | stdv::keys) out.push_back(key);
  return out;
}

Size SpectralField::size() const { return props.size(); }
}  // namespace Subsurface

static_assert(spectral_field<SubsurfaceSpectralField, SubsurfaceSpectralPropertyTag, Numeric, Numeric, Numeric>,
              "The subsurface spectral field must satisfy the spectral field concepts");
