#include "surf_spectral_field.h"

#include <debug.h>
#include <field_concepts.h>
#include <lagrange_interp.h>

#include <optional>
#include <string_view>
#include <utility>

namespace Surf {
SpectralData::SpectralData()                                    = default;
SpectralData::SpectralData(const SpectralData &)                = default;
SpectralData::SpectralData(SpectralData &&) noexcept            = default;
SpectralData &SpectralData::operator=(const SpectralData &)     = default;
SpectralData &SpectralData::operator=(SpectralData &&) noexcept = default;

void SpectralData::adjust_interpolation_extrapolation() {
  using enum InterpolationExtrapolation;

  if (const auto *field = std::get_if<SpectralGeodeticField3>(&data)) {
    if (field->grid<0>().size() == 1) frq_upp = frq_low = Nearest;
    if (field->grid<1>().size() == 1) lat_upp = lat_low = Nearest;
    if (field->grid<2>().size() == 1) lon_upp = lon_low = Nearest;
  } else {
    frq_upp = frq_low = lat_upp = lat_low = lon_upp = lon_low = Nearest;
  }
}

SpectralData::SpectralData(Numeric x) : data(x) { adjust_interpolation_extrapolation(); }

SpectralData &SpectralData::operator=(Numeric x) {
  data = x;
  adjust_interpolation_extrapolation();
  return *this;
}

SpectralData::SpectralData(SpectralGeodeticField3 x) : data(std::move(x)) { adjust_interpolation_extrapolation(); }

SpectralData &SpectralData::operator=(SpectralGeodeticField3 x) {
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
  if (std::holds_alternative<SpectralGeodeticField3>(data)) return "SpectralGeodeticField3";
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
  if (grid.empty()) arts::user_error("Cannot interpolate on an empty grid of the surface spectral data");
  if (grid.size() == 1) return grid.template lag<0, transform>(x);
  return grid.template lag<1, transform>(x);
}

const SpectralGeodeticField3 &checked(const SpectralGeodeticField3 &field) {
  if (not field.ok())
    arts::user_error(
        "The surface spectral data shape {:B,} does not match its grid sizes {:B,}", field.data.shape(), field.shape());
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
          "The {} {} is outside the grid [{}, {}] of the surface spectral data, which does not extrapolate there",
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
                                             const SpectralGeodeticField3 &field,
                                             const Numeric                 frq) {
  const auto f = extrapolate(field.grid<0>(), frq, d.frq_low, d.frq_upp, "frequency"sv);
  if (not f) return std::nullopt;
  return lag<grid_identity>(field.grid<0>(), *f);
}

struct spatial_lags {
  lags<grid_identity> lat;
  lags<loncross>      lon;
};

//! All coordinates are extrapolated first, so that None throws even where another coordinate gives zero
std::optional<spatial_lags> spatial(const SpectralData           &d,
                                    const SpectralGeodeticField3 &field,
                                    const Numeric                 lat,
                                    const Numeric                 lon) {
  const auto a = extrapolate(field.grid<1>(), lat, d.lat_low, d.lat_upp, "latitude"sv);
  const auto b = extrapolate(field.grid<2>(), lon, d.lon_low, d.lon_upp, "longitude"sv);
  if (not a or not b) return std::nullopt;
  return spatial_lags{.lat = lag<grid_identity>(field.grid<1>(), *a), .lon = lag<loncross>(field.grid<2>(), *b)};
}

//! The spatially interpolated data at the frequency index k
Numeric spatial_value(const SpectralGeodeticField3 &field, const Index k, const spatial_lags &s) {
  return std::visit(
      [data = field.data[k]](const auto &a, const auto &b) { return lagrange_interp::interp(data, a, b); },
      s.lat,
      s.lon);
}

//! The frequency interpolation of the spatially interpolated data, value(k) at frequency index k
Numeric frequency_value(const lags<grid_identity> &f, const auto &value) {
  return std::visit(
      [&value](const auto &l) {
        Numeric out = 0.0;
        for (Size j = 0; j < l.size(); j++) {
          // The term is its own expression, so that it is not contracted with the sum
          const Numeric term  = l.data[j] * value(l.indx[j]);
          out                += term;
        }
        return out;
      },
      f);
}

//! The number of stored values per frequency
Index spatial_size(const SpectralGeodeticField3 &field) {
  return static_cast<Index>(field.grid<1>().size() * field.grid<2>().size());
}
}  // namespace

Numeric SpectralData::at(const Numeric frq, const Numeric lat, const Numeric lon) const {
  return std::visit(
      [&]<typename T>(const T &x) -> Numeric {
        if constexpr (std::same_as<T, SpectralGeodeticField3>) {
          const SpectralGeodeticField3 &field = checked(x);
          const auto                    f     = frequency(*this, field, frq);
          const auto                    s     = spatial(*this, field, lat, lon);
          if (not f or not s) return 0.0;
          return frequency_value(*f, [&](Index k) { return spatial_value(field, k, *s); });
        } else if constexpr (std::same_as<T, Numeric>) {
          return x;
        } else {
          return x(frq, lat, lon);
        }
      },
      data);
}

Vector SpectralData::at(const AscendingGrid &freq_grid, const Numeric lat, const Numeric lon) const {
  Vector out(freq_grid.size());

  std::visit(
      [&]<typename T>(const T &x) {
        if constexpr (std::same_as<T, SpectralGeodeticField3>) {
          const SpectralGeodeticField3 &field = checked(x);
          const auto                    s     = spatial(*this, field, lat, lon);

          // The data at each of its frequencies, interpolated to the position once
          Vector values;
          if (s) {
            values.resize(field.grid<0>().size());
            for (Size k = 0; k < values.size(); k++) values[k] = spatial_value(field, k, *s);
          }

          for (Size i = 0; i < freq_grid.size(); i++) {
            const auto f = frequency(*this, field, freq_grid[i]);
            out[i]       = (f and s) ? frequency_value(*f, [&values](Index k) { return values[k]; }) : 0.0;
          }
        } else if constexpr (std::same_as<T, Numeric>) {
          out = x;
        } else {
          for (Size i = 0; i < freq_grid.size(); i++) out[i] = x(freq_grid[i], lat, lon);
        }
      },
      data);

  return out;
}

void SpectralData::flat_weight(std::vector<std::pair<Index, Numeric>> &weights,
                               const Numeric                           frq,
                               const Numeric                           lat,
                               const Numeric                           lon) const {
  weights.clear();

  std::visit(
      [&]<typename T>(const T &x) {
        if constexpr (std::same_as<T, SpectralGeodeticField3>) {
          const SpectralGeodeticField3 &field = checked(x);
          const auto                    f     = frequency(*this, field, frq);
          const auto                    s     = spatial(*this, field, lat, lon);
          if (not f or not s) return;

          const Index n    = spatial_size(field);
          const Index nlon = field.grid<2>().size();

          std::visit(
              [&](const auto &l, const auto &a, const auto &b) {
                const auto w = lagrange_interp::interpweights(a, b);
                for (Size m = 0; m < l.size(); m++) {
                  if (l.data[m] == 0.0) continue;
                  for (Size i = 0; i < a.size(); i++) {
                    for (Size j = 0; j < b.size(); j++) {
                      if (w[i, j] == 0.0) continue;
                      weights.emplace_back(l.indx[m] * n + a.indx[i] * nlon + b.indx[j], l.data[m] * w[i, j]);
                    }
                  }
                }
              },
              *f,
              s->lat,
              s->lon);
        } else if constexpr (std::same_as<T, Numeric>) {
          weights.emplace_back(0, 1.0);
        }
      },
      data);
}

ConstVectorView SpectralData::flat_view() const {
  return std::visit(
      []<typename T>(const T &x) -> ConstVectorView {
        if constexpr (std::same_as<T, SpectralGeodeticField3>)
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
        if constexpr (std::same_as<T, SpectralGeodeticField3>)
          return x.data.view_as(x.data.size());
        else if constexpr (std::same_as<T, Numeric>)
          return VectorView{x};
        else
          return VectorView{};
      },
      data);
}

bool SpectralData::ok() const {
  if (const auto *v = std::get_if<SpectralGeodeticField3>(&data)) {
    return v->ok() and loncross::cycle(v->grid<2>().front()) == v->grid<2>().front() and
           loncross::cycle(v->grid<2>().back()) == v->grid<2>().back();
  }

  return true;
}

SpectralData &SpectralField::operator[](const SurfaceSpectralPropertyTag &key) { return props[key]; }

const SpectralData &SpectralField::operator[](const SurfaceSpectralPropertyTag &key) const {
  const auto ptr = props.find(key);
  if (ptr == props.end()) arts::user_error("The surface spectral field has no {}", key);
  return ptr->second;
}

bool SpectralField::contains(const SurfaceSpectralPropertyTag &key) const { return props.contains(key); }

std::vector<SurfaceSpectralPropertyTag> SpectralField::keys() const {
  std::vector<SurfaceSpectralPropertyTag> out;
  out.reserve(props.size());
  for (const auto &key : props | stdv::keys) out.push_back(key);
  return out;
}

Size SpectralField::size() const { return props.size(); }
}  // namespace Surf

static_assert(spectral_field<SurfaceSpectralField, SurfaceSpectralPropertyTag, Numeric, Numeric>,
              "The surface spectral field must satisfy the spectral field concepts");
