#include "atm_field.h"

#include <compare.h>
#include <configtypes.h>
#include <debug.h>
#include <enumsFieldComponent.h>
#include <functional_atm_field_interp.h>
#include <hitran_species.h>
#include <isotopologues.h>
#include <matpack.h>
#include <physics_funcs.h>

#include <algorithm>
#include <limits>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

bool AtmPoint::zero_wind() const noexcept {
  return std::all_of(wind.begin(), wind.end(), Cmp::eq(0));
}

bool AtmPoint::zero_mag() const noexcept {
  return std::all_of(mag.begin(), mag.end(), Cmp::eq(0));
}

void Atm::Data::adjust_interpolation_extrapolation() {
  if (std::holds_alternative<SortedGriddedField3>(data)) {
    auto &field = std::get<SortedGriddedField3>(data);

    if (field.grid<0>().size() == 1) {
      alt_upp = InterpolationExtrapolation::Nearest;
      alt_low = InterpolationExtrapolation::Nearest;
    }

    if (field.grid<1>().size() == 1) {
      lat_upp = InterpolationExtrapolation::Nearest;
      lat_low = InterpolationExtrapolation::Nearest;
    }

    if (field.grid<2>().size() == 1) {
      lon_upp = InterpolationExtrapolation::Nearest;
      lon_low = InterpolationExtrapolation::Nearest;
    }
  } else {
    alt_upp = InterpolationExtrapolation::Nearest;
    alt_low = InterpolationExtrapolation::Nearest;
    lat_upp = InterpolationExtrapolation::Nearest;
    lat_low = InterpolationExtrapolation::Nearest;
    lon_upp = InterpolationExtrapolation::Nearest;
    lon_low = InterpolationExtrapolation::Nearest;
  }
}

AtmKey to_wind(const String &x) {
  switch (to<FieldComponent>(x)) {
    case FieldComponent::u: return AtmKey::wind_u;
    case FieldComponent::v: return AtmKey::wind_v;
    case FieldComponent::w: return AtmKey::wind_w;
  }
  std::unreachable();
}

AtmKey to_mag(const String &x) {
  switch (to<FieldComponent>(x)) {
    case FieldComponent::u: return AtmKey::mag_u;
    case FieldComponent::v: return AtmKey::mag_v;
    case FieldComponent::w: return AtmKey::mag_w;
  }
  std::unreachable();
}

Numeric AtmPoint::number_density() const {
  return ::number_density(pressure, temperature);
}

Numeric AtmPoint::number_density(const SpeciesEnum &spec) const {
  return specs.at(spec) * number_density();
}

Numeric AtmPoint::number_density(const SpeciesIsotope &spec) const {
  return isots.at(spec) * number_density(spec.spec);
}

Numeric AtmPoint::operator[](SpeciesEnum x) const try {
  return specs.at(x);
} catch (std::out_of_range &) {
  ARTS_USER_ERROR("Species VMR not found: \"{}\"", toString<1>(x))
}

Numeric AtmPoint::operator[](const SpeciesIsotope &x) const try {
  return isots.at(x);
} catch (std::out_of_range &) {
  ARTS_USER_ERROR("Isotopologue ratio not found: \"{}\"", x)
}

Numeric AtmPoint::operator[](const QuantumLevelIdentifier &x) const try {
  return nlte.at(x);
} catch (std::out_of_range &) {
  ARTS_USER_ERROR("QuantumLevelIdentifier not found: \"{}\"", x)
}

Numeric AtmPoint::operator[](const ScatteringSpeciesProperty &x) const try {
  return ssprops.at(x);
} catch (std::out_of_range &) {
  ARTS_USER_ERROR("ScatteringSpeciesProperty not found: \"", x, '"')
}

Numeric AtmPoint::operator[](AtmKey x) const {
  switch (x) {
    case AtmKey::t:      return temperature;
    case AtmKey::p:      return pressure;
    case AtmKey::wind_u: return wind[0];
    case AtmKey::wind_v: return wind[1];
    case AtmKey::wind_w: return wind[2];
    case AtmKey::mag_u:  return mag[0];
    case AtmKey::mag_v:  return mag[1];
    case AtmKey::mag_w:  return mag[2];
  }
  std::unreachable();
}

Numeric &AtmPoint::operator[](SpeciesEnum x) { return specs[x]; }

Numeric &AtmPoint::operator[](const SpeciesIsotope &x) { return isots[x]; }

Numeric &AtmPoint::operator[](const QuantumLevelIdentifier &x) {
  return nlte[x];
}

Numeric &AtmPoint::operator[](const ScatteringSpeciesProperty &x) {
  return ssprops[x];
}

Numeric &AtmPoint::operator[](AtmKey x) {
  switch (x) {
    case AtmKey::t:      return temperature;
    case AtmKey::p:      return pressure;
    case AtmKey::wind_u: return wind[0];
    case AtmKey::wind_v: return wind[1];
    case AtmKey::wind_w: return wind[2];
    case AtmKey::mag_u:  return mag[0];
    case AtmKey::mag_v:  return mag[1];
    case AtmKey::mag_w:  return mag[2];
  }
  std::unreachable();
}

std::pair<Numeric, Numeric> AtmPoint::levels(
    const QuantumIdentifier &band) const {
  return {operator[](band.lower()), operator[](band.upper())};
}

namespace Atm {
Point::Point(const IsoRatioOption isots_key) {
  switch (isots_key) {
    case IsoRatioOption::Builtin: {
      const SpeciesIsotopologueRatios x =
          Species::isotopologue_ratiosInitFromBuiltin();
      for (Index i = 0; i < x.maxsize; i++) {
        if (Species::Isotopologues[i].is_joker()) continue;
        if (Species::Isotopologues[i].is_predefined()) continue;
        isots[Species::Isotopologues[i]] = x.data[i];
      }
    } break;
    case IsoRatioOption::Hitran: {
      const SpeciesIsotopologueRatios x = Hitran::isotopologue_ratios();
      for (Index i = 0; i < x.maxsize; i++) {
        if (Species::Isotopologues[i].is_joker()) continue;
        if (Species::Isotopologues[i].is_predefined()) continue;
        isots[Species::Isotopologues[i]] = x.data[i];
      }
    } break;
    case IsoRatioOption::None:
    default:                   break;
  }
}

Field::Field(const IsoRatioOption isots_key) {
  switch (isots_key) {
    case IsoRatioOption::Builtin: {
      const SpeciesIsotopologueRatios x =
          Species::isotopologue_ratiosInitFromBuiltin();
      for (Index i = 0; i < x.maxsize; i++) {
        if (Species::Isotopologues[i].is_joker()) continue;
        if (Species::Isotopologues[i].is_predefined()) continue;
        isots[Species::Isotopologues[i]] = x.data[i];
      }
    } break;
    case IsoRatioOption::Hitran: {
      const SpeciesIsotopologueRatios x = Hitran::isotopologue_ratios();
      for (Index i = 0; i < x.maxsize; i++) {
        if (Species::Isotopologues[i].is_joker()) continue;
        if (Species::Isotopologues[i].is_predefined()) continue;
        isots[Species::Isotopologues[i]] = x.data[i];
      }
    } break;
    case IsoRatioOption::None: break;
  }
}

Field::Field() : Field(IsoRatioOption::Builtin) {}
Point::Point() : Point(IsoRatioOption::Builtin) {}

Field::Field(const Field &)                = default;
Field::Field(Field &&) noexcept            = default;
Field &Field::operator=(const Field &)     = default;
Field &Field::operator=(Field &&) noexcept = default;
Point::Point(Numeric p, Numeric t) : pressure(p), temperature(t) {}
Point::Point(const Point &)                = default;
Point::Point(Point &&) noexcept            = default;
Point &Point::operator=(const Point &)     = default;
Point &Point::operator=(Point &&) noexcept = default;
Data::Data()                               = default;
Data::Data(const Data &)                   = default;
Data::Data(Data &&) noexcept               = default;
Data &Data::operator=(const Data &)        = default;
Data &Data::operator=(Data &&) noexcept    = default;

bool Field::contains(const AtmKey &key) const { return other.contains(key); }

bool Field::contains(const SpeciesEnum &key) const {
  return specs.contains(key);
}

bool Field::contains(const SpeciesIsotope &key) const {
  return isots.contains(key);
}

bool Field::contains(const QuantumLevelIdentifier &key) const {
  return nlte.contains(key);
}

bool Field::contains(const ScatteringSpeciesProperty &key) const {
  return ssprops.contains(key);
}

bool Field::contains(const KeyVal &key) const {
  return std::visit([this](auto &x) { return contains(x); }, key);
}

bool Point::contains(const KeyVal &key) const {
  return std::visit([this](auto &x) { return has(x); }, key);
}

Numeric Point::mean_mass(SpeciesEnum s) const {
  Numeric ratio = 0.0;
  Numeric mass  = 0.0;
  for (auto &[isot, this_ratio] : isots) {
    if (isot.spec == s and not(isot.is_predefined() or isot.is_joker())) {
      ratio += this_ratio;
      mass  += this_ratio * isot.mass;
    }
  }

  ARTS_USER_ERROR_IF(ratio == 0,
                     "Cannot find a ratio for the mean mass of species \"{}\"",
                     toString<1>(s))

  return mass / ratio;
}

Numeric Point::mean_mass() const {
  Numeric vmr  = 0.0;
  Numeric mass = 0.0;
  for (auto &[spec, this_vmr] : specs) {
    vmr += this_vmr;
    if (this_vmr != 0.0) {
      mass += this_vmr * mean_mass(spec);
    }
  }

  ARTS_USER_ERROR_IF(vmr == 0,
                     "Cannot find a ratio for the mean mass of the atmosphere")

  return mass / vmr;
}

std::vector<KeyVal> Point::keys(bool keep_basic,
                                bool keep_specs,
                                bool keep_isots,
                                bool keep_nlte,
                                bool keep_ssprops) const {
  std::vector<KeyVal> out;
  out.reserve(size());

  if (keep_basic) {
    for (auto &a : enumtyps::AtmKeyTypes) out.emplace_back(a);
  }

  if (keep_specs) {
    for (auto &a : specs) out.emplace_back(a.first);
  }

  if (keep_nlte) {
    for (auto &a : nlte) out.emplace_back(a.first);
  }

  if (keep_ssprops) {
    for (auto &a : ssprops) out.emplace_back(a.first);
  }

  if (keep_isots) {
    for (auto &a : isots) out.emplace_back(a.first);
  }
  return out;
}

Index Point::nspec() const { return static_cast<Index>(specs.size()); }

Index Point::npart() const { return static_cast<Index>(ssprops.size()); }

Index Point::nisot() const { return static_cast<Index>(isots.size()); }

Index Point::nnlte() const { return static_cast<Index>(nlte.size()); }

Index Point::size() const {
  return nspec() + nnlte() + nother() + npart() + nisot();
}

Size Field::nspec() const { return specs.size(); }

Size Field::nisot() const { return isots.size(); }

Size Field::npart() const { return ssprops.size(); }

Size Field::nnlte() const { return nlte.size(); }

Size Field::nother() const { return other.size(); }

Size Field::size() const {
  return nspec() + nnlte() + nother() + npart() + nisot();
}

String Data::data_type() const {
  if (std::holds_alternative<SortedGriddedField3>(data))
    return "SortedGriddedField3";
  if (std::holds_alternative<Numeric>(data)) return "Numeric";
  if (std::holds_alternative<FunctionalData>(data)) return "FunctionalData";

  ARTS_USER_ERROR("Cannot understand data type; is this a new type")
}

namespace {
struct Limits {
  Numeric alt_low{std::numeric_limits<Numeric>::lowest()};
  Numeric alt_upp{std::numeric_limits<Numeric>::max()};
  Numeric lat_low{std::numeric_limits<Numeric>::lowest()};
  Numeric lat_upp{std::numeric_limits<Numeric>::max()};
  Numeric lon_low{std::numeric_limits<Numeric>::lowest()};
  Numeric lon_upp{std::numeric_limits<Numeric>::max()};
};

struct ComputeLimit {
  InterpolationExtrapolation type{InterpolationExtrapolation::Linear};
  Numeric alt, lat, lon;
};

Limits find_limits(const Numeric &) { return {}; }

Limits find_limits(const FunctionalData &) { return {}; }

Limits find_limits(const SortedGriddedField3 &gf3) {
  return {.alt_low = gf3.grid<0>().front(),
          .alt_upp = gf3.grid<0>().back(),
          .lat_low = gf3.grid<1>().front(),
          .lat_upp = gf3.grid<1>().back(),
          .lon_low = gf3.grid<2>().front(),
          .lon_upp = gf3.grid<2>().back()};
}

template <Index poly_alt, Index poly_lat, Index poly_lon>
struct interp_helper {
  using AltLag = lagrange_interp::lag_t<poly_alt>;
  using LatLag = lagrange_interp::lag_t<poly_lat>;
  using LonLag = lagrange_interp::lag_t<poly_lon, lagrange_interp::loncross>;

  Array<AltLag> lags_alt;
  Array<LatLag> lags_lat;
  Array<LonLag> lags_lon;

  constexpr interp_helper(const Vector &alt_grid,
                          const Vector &lat_grid,
                          const Vector &lon_grid,
                          const Vector &alt,
                          const Vector &lat,
                          const Vector &lon)
      : lags_alt(lagrange_interp::make_lags<poly_alt>(
            alt_grid, alt, -1, "Altitude")),
        lags_lat(lagrange_interp::make_lags<poly_lat>(
            lat_grid, lat, -1, "Latitude")),
        lags_lon(
            lagrange_interp::make_lags<poly_lon, lagrange_interp::loncross>(
                lon_grid, lon, -1, "Longitude")) {}

  Vector operator()(const Tensor3 &data) const {
    return flat_interp(data, lags_alt, lags_lat, lags_lon);
  }
};

template <Index poly_alt, Index poly_lat, Index poly_lon>
interp_helper<poly_alt, poly_lat, poly_lon> tvec_interpgrid(
    const Vector &alt_grid,
    const Vector &lat_grid,
    const Vector &lon_grid,
    const Vector &alt,
    const Vector &lat,
    const Vector &lon) {
  using Interpolater = interp_helper<poly_alt, poly_lat, poly_lon>;
  return Interpolater(alt_grid, lat_grid, lon_grid, alt, lat, lon);
}

template <Index poly_alt, Index poly_lat, Index poly_lon>
Array<std::array<std::pair<Index, Numeric>, 8>> tvec_interpgrid_weights(
    const Vector &alt_grid,
    const Vector &lat_grid,
    const Vector &lon_grid,
    const Vector &alt,
    const Vector &lat,
    const Vector &lon) {
  const auto interpolater = tvec_interpgrid<poly_alt, poly_lat, poly_lon>(
      alt_grid, lat_grid, lon_grid, alt, lat, lon);

  constexpr auto v0 = std::pair<Index, Numeric>{0, 0.};
  constexpr std::array v{v0, v0, v0, v0, v0, v0, v0, v0};
  Array<std::array<std::pair<Index, Numeric>, 8>> out(alt.size(), v);

  const Size n = lat_grid.size();
  const Size m = lon_grid.size();

  for (Size i = 0; i < alt.size(); i++) {
    const Size alt0 = interpolater.lags_alt[i].pos * m * n;
    const Size lat0 = interpolater.lags_lat[i].pos * m;
    const Size lon0 = interpolater.lags_lon[i].pos;

    Size j = 0;
    for (Size idx0 = 0; idx0 < 1 + poly_alt; idx0++) {
      for (Size idx1 = 0; idx1 < 1 + poly_lat; idx1++) {
        for (Size idx2 = 0; idx2 < 1 + poly_lon; idx2++) {
          out[i][j].first = alt0 + lat0 + lon0 + idx0 * n * m + idx1 * m + idx2;
          out[i][j].second = interpolater.lags_alt[i].lx[idx0] *
                             interpolater.lags_lat[i].lx[idx1] *
                             interpolater.lags_lon[i].lx[idx2];
          ++j;
        }
      }
    }
  }

  return out;
}

template <Index poly_alt, Index poly_lat, Index poly_lon>
Vector tvec_interp(const Tensor3 &v,
                   const Vector &alt_grid,
                   const Vector &lat_grid,
                   const Vector &lon_grid,
                   const Vector &alt,
                   const Vector &lat,
                   const Vector &lon) {
  const auto interpolater = tvec_interpgrid<poly_alt, poly_lat, poly_lon>(
      alt_grid, lat_grid, lon_grid, alt, lat, lon);
  return interpolater(v);
}

constexpr InterpolationExtrapolation combine(InterpolationExtrapolation a,
                                             InterpolationExtrapolation b) {
  using enum InterpolationExtrapolation;
  switch (a) {
    case None: return None;
    case Zero: {
      switch (b) {
        case None:    return None;
        case Zero:    return Zero;
        case Nearest: return Zero;
        case Linear:  return Zero;
      }
      std::unreachable();
    }
    case Nearest: {
      switch (b) {
        case None:    return None;
        case Zero:    return Zero;
        case Nearest: return Nearest;
        case Linear:  return Nearest;
      }
      std::unreachable();
    }
    case Linear: return b;
  }

  std::unreachable();
}

constexpr InterpolationExtrapolation combine(InterpolationExtrapolation a,
                                             InterpolationExtrapolation b,
                                             InterpolationExtrapolation c) {
  return combine(combine(a, b), c);
}

void select(InterpolationExtrapolation lowt,
            InterpolationExtrapolation uppt,
            Numeric lowv,
            Numeric uppv,
            Numeric v,
            Numeric &outv,
            InterpolationExtrapolation &outt) {
  if (v < lowv) {
    outt = lowt;
    if (outt == InterpolationExtrapolation::Nearest) v = lowv;
  } else if (uppv < v) {
    outt = uppt;
    if (outt == InterpolationExtrapolation::Nearest) v = uppv;
  }

  outv = v;
}

ComputeLimit find_limit(const Data &data,
                        const Limits &lim,
                        Numeric alt,
                        Numeric lat,
                        Numeric lon) {
  ComputeLimit out;
  InterpolationExtrapolation a{InterpolationExtrapolation::Linear},
      b{InterpolationExtrapolation::Linear},
      c{InterpolationExtrapolation::Linear};

  select(data.alt_low, data.alt_upp, lim.alt_low, lim.alt_upp, alt, out.alt, a);
  select(data.lat_low, data.lat_upp, lim.lat_low, lim.lat_upp, lat, out.lat, b);
  select(data.lon_low, data.lon_upp, lim.lon_low, lim.lon_upp, lon, out.lon, c);
  out.type = combine(a, b, c);

  return out;
}
}  // namespace

bool Point::is_lte() const noexcept { return nlte.empty(); }

namespace {
template <class Key, class T, class Hash, class KeyEqual, class Allocator>
std::vector<AtmKey> get_keys(
    const std::unordered_map<Key, T, Hash, KeyEqual, Allocator> &map) {
  std::vector<Key> out(map.size());
  std::transform(
      map.begin(), map.end(), out.begin(), [](auto &v) { return v.first; });
  return out;
}
}  // namespace

std::vector<KeyVal> Field::keys() const {
  std::vector<KeyVal> out;
  out.reserve(size());

  for (const auto &a : other) out.emplace_back(a.first);
  for (const auto &a : specs) out.emplace_back(a.first);
  for (const auto &a : isots) out.emplace_back(a.first);
  for (const auto &a : nlte) out.emplace_back(a.first);
  for (const auto &a : ssprops) out.emplace_back(a.first);

  return out;
}

void Point::check_and_fix() try {
  ARTS_USER_ERROR_IF(nonstd::isnan(pressure), "Pressure is NaN")
  ARTS_USER_ERROR_IF(nonstd::isnan(temperature), "Temperature is NaN")

  if (std::ranges::all_of(wind, [](auto v) { return nonstd::isnan(v); })) {
    wind = {0., 0., 0.};
  } else {
    ARTS_USER_ERROR_IF(
        std::ranges::any_of(wind, [](auto v) { return nonstd::isnan(v); }),
        "Cannot have partially missing wind field.  Consider setting the missing field to zero or add it completely.\n"
        "Wind field [wind_u, wind_v, wind_w] is: {:B,}",
        wind)
  }

  if (std::ranges::all_of(mag, [](auto v) { return nonstd::isnan(v); })) {
    mag = {0., 0., 0.};
  } else {
    ARTS_USER_ERROR_IF(
        std::ranges::any_of(mag, [](auto v) { return nonstd::isnan(v); }),
        "Cannot have partially missing magnetic field.  Consider setting the missing field to zero or add it completely.\n"
        "Magnetic field [mag_u, mag_v, mag_w] is: {:B,}",
        mag)
  }

  for (auto &spec : specs) {
    ARTS_USER_ERROR_IF(nonstd::isnan(spec.second) or spec.second < 0.0,
                       "VMR for \"{}\" is {}",
                       toString<1>(spec.first),
                       spec.second)
  }

  for (auto &isot : isots) {
    //! Cannot check isnan because it is a valid state for isotopologue ratios
    ARTS_USER_ERROR_IF(isot.second < 0.0,
                       "Isotopologue ratio for \"{}\" is {}",
                       isot.first.FullName(),
                       isot.second)
  }

  for (auto &nl : nlte) {
    ARTS_USER_ERROR_IF(nonstd::isnan(nl.second) or nl.second < 0.0,
                       "Non-LTE ratio for \"{}\" is {}",
                       nl.first,
                       nl.second)
  }

  for (auto &pp : ssprops) {
    ARTS_USER_ERROR_IF(nonstd::isnan(pp.second),
                       "Scattering Species Property value for \"",
                       pp.first,
                       pp.second)
  }
}
ARTS_METHOD_ERROR_CATCH

Data &Field::operator[](const AtmKey &key) { return other[key]; }

Data &Field::operator[](const SpeciesEnum &key) { return specs[key]; }

Data &Field::operator[](const SpeciesIsotope &key) { return isots[key]; }

Data &Field::operator[](const QuantumLevelIdentifier &key) { return nlte[key]; }

Data &Field::operator[](const ScatteringSpeciesProperty &key) {
  return ssprops[key];
}

Data &Field::operator[](const KeyVal &key) {
  return std::visit(
      [this](auto &key) -> Data & { return this->operator[](key); }, key);
}

const Data &Field::operator[](const AtmKey &key) const { return other.at(key); }

const Data &Field::operator[](const SpeciesEnum &key) const {
  return specs.at(key);
}

const Data &Field::operator[](const SpeciesIsotope &key) const {
  return isots.at(key);
}

const Data &Field::operator[](const QuantumLevelIdentifier &key) const {
  return nlte.at(key);
}

const Data &Field::operator[](const ScatteringSpeciesProperty &key) const {
  return ssprops.at(key);
}

const Data &Field::operator[](const KeyVal &key) const {
  return std::visit(
      [this](auto &key) -> const Data & { return this->operator[](key); }, key);
}

Numeric Point::operator[](const KeyVal &k) const {
  return std::visit([this](auto &key) { return this->operator[](key); }, k);
}

Numeric &Point::operator[](const KeyVal &k) {
  return std::visit(
      [this](auto &key) -> Numeric & {
        return const_cast<Point *>(this)->operator[](key);
      },
      k);
}

ConstVectorView Data::flat_view() const {
  return std::visit(
      [](auto &X) -> ConstVectorView {
        using T = std::remove_cvref_t<decltype(X)>;
        if constexpr (std::same_as<T, SortedGriddedField3>)
          return X.data.view_as(X.data.size());
        else if constexpr (std::same_as<T, Numeric>)
          return ConstVectorView{X};
        else if constexpr (std::same_as<T, FunctionalData>)
          return ternary::get_xc(X.f);
        std::unreachable();
      },
      data);
}

VectorView Data::flat_view() {
  return std::visit(
      [](auto &X) -> VectorView {
        using T = std::remove_cvref_t<decltype(X)>;
        if constexpr (std::same_as<T, SortedGriddedField3>)
          return X.data.view_as(X.data.size());
        else if constexpr (std::same_as<T, Numeric>)
          return VectorView{X};
        else if constexpr (std::same_as<T, FunctionalData>)
          return ternary::get_xm(X.f);
        std::unreachable();
      },
      data);
}

namespace {
std::vector<std::pair<Index, Numeric>> flat_weight_(
    const SortedGriddedField3 &data, Numeric alt, Numeric lat, Numeric lon) {
  return interp::flat_weight(data, alt, lat, lon);
}

std::vector<std::pair<Index, Numeric>> flat_weight_(Numeric,
                                                    Numeric,
                                                    Numeric,
                                                    Numeric) {
  constexpr std::pair<Index, Numeric> v1{0, 1.0};
  return {v1};
}

std::vector<std::pair<Index, Numeric>> flat_weight_(const FunctionalData &f,
                                                    Numeric alt,
                                                    Numeric lat,
                                                    Numeric lon) {
  return ternary::get_w(f.f, alt, lat, lon);
}
}  // namespace

//! Flat weights for the positions in an atmosphere
std::vector<std::pair<Index, Numeric>> Data::flat_weight(
    const Numeric alt, const Numeric lat, const Numeric lon) const {
  return std::visit([&](auto &v) { return flat_weight_(v, alt, lat, lon); },
                    data);
}

std::vector<std::pair<Index, Numeric>> Data::flat_weight(
    const Vector3 pos) const {
  return flat_weight(pos[0], pos[1], pos[2]);
}

bool Data::ok() const {
  if (std::holds_alternative<SortedGriddedField3>(data)) {
    auto &v = *std::get_if<SortedGriddedField3>(&data);
    return v.ok() and
           lagrange_interp::loncross::cycle(v.grid<2>().front()) ==
               v.grid<2>().front() and
           lagrange_interp::loncross::cycle(v.grid<2>().back()) ==
               v.grid<2>().back() and
           v.grid<1>().front() >= -90 and v.grid<1>().back() <= 90;
  }

  return true;
}

void Data::fix_cyclicity() {
  if (std::holds_alternative<SortedGriddedField3>(data)) {
    if (ok()) return;

    auto &v              = *std::get_if<SortedGriddedField3>(&data);
    const auto &lon_grid = v.grid<2>();

    std::vector<std::pair<Index, Numeric>> lon{lon_grid.size()};
    for (Size i = 0; i < lon.size(); ++i) {
      lon[i] = {i, lon_grid[i]};
      while (lon[i].second != lagrange_interp::loncross::cycle(lon[i].second)) {
        lon[i].second = lagrange_interp::loncross::cycle(lon[i].second);
      }
    }

    stdr::sort(lon, {}, &std::pair<Index, Numeric>::second);
    auto [end, _] = stdr::unique(lon, {}, &std::pair<Index, Numeric>::second);
    lon.erase(end, lon.end());

    Tensor3 data(v.shape()[0], v.shape()[1], lon.size());
    Vector lon_grid_new(lon.size());

    for (auto [i, l] : lon) {
      lon_grid_new[i] = l;
      data[joker, i]  = v.data[joker, joker, i];
    }

    v.data      = std::move(data);
    v.grid<2>() = std::move(lon_grid_new);
  }
}

Data::Data(Numeric x) : data(x) { adjust_interpolation_extrapolation(); }

Data::Data(SortedGriddedField3 x) : data(std::move(x)) {
  adjust_interpolation_extrapolation();
}

Data::Data(FunctionalData x) : data(std::move(x)) {
  adjust_interpolation_extrapolation();
}

Data &Data::operator=(Numeric x) {
  data = x;
  adjust_interpolation_extrapolation();
  return *this;
}

Data &Data::operator=(SortedGriddedField3 x) {
  data = std::move(x);
  adjust_interpolation_extrapolation();
  return *this;
}

Data &Data::operator=(FunctionalData x) {
  data = std::move(x);
  adjust_interpolation_extrapolation();
  return *this;
}
}  // namespace Atm

namespace {
template <Atm::KeyType T>
constexpr bool cmp(const AtmKeyVal &keyval, const T &key) {
  const auto *ptr = std::get_if<T>(&keyval);
  return ptr and *ptr == key;
}
}  // namespace

bool operator==(const AtmKeyVal &keyval, AtmKey key) {
  return cmp(keyval, key);
}

bool operator==(AtmKey key, const AtmKeyVal &keyval) {
  return cmp(keyval, key);
}

bool operator==(const AtmKeyVal &keyval, const SpeciesEnum &key) {
  return cmp(keyval, key);
}

bool operator==(const SpeciesEnum &key, const AtmKeyVal &keyval) {
  return cmp(keyval, key);
}

bool operator==(const AtmKeyVal &keyval, const QuantumLevelIdentifier &key) {
  return cmp(keyval, key);
}

bool operator==(const QuantumLevelIdentifier &key, const AtmKeyVal &keyval) {
  return cmp(keyval, key);
}

bool operator==(const ScatteringSpeciesProperty &key, const AtmKeyVal &keyval) {
  return cmp(keyval, key);
}

namespace Atm {
namespace interp {
namespace {
struct PositionalNumeric {
  const FieldData &data;
  const Numeric alt;
  const Numeric lat;
  const Numeric lon;

  operator Numeric() const {
    return std::visit([&](auto &d) { return get(d, alt, lat, lon); }, data);
  }
};

std::optional<Numeric> get_optional_limit(const Data &data,
                                          const Numeric alt,
                                          const Numeric lat,
                                          const Numeric lon) {
  const auto lim = Atm::find_limit(
      data,
      std::visit([](auto &d) { return Atm::find_limits(d); }, data.data),
      alt,
      lat,
      lon);

  ARTS_USER_ERROR_IF(
      lim.type == InterpolationExtrapolation::None,
      "Limit breached.  Position ({}, {}, {}) is out-of-bounds when no extrapolation is wanted",
      lim.alt,
      lim.lat,
      lim.lon)

  if (lim.type == InterpolationExtrapolation::Zero) return 0.0;

  if (lim.type == InterpolationExtrapolation::Nearest)
    return PositionalNumeric{
        .data = data.data, .alt = lim.alt, .lat = lim.lat, .lon = lim.lon};

  return std::nullopt;
}
}  // namespace
}  // namespace interp

Numeric Data::at(const Numeric alt,
                 const Numeric lat,
                 const Numeric lon) const {
  return interp::get_optional_limit(*this, alt, lat, lon)
      .value_or(interp::PositionalNumeric{data, alt, lat, lon});
}

Numeric Data::at(const Vector3 pos) const { return at(pos[0], pos[1], pos[2]); }

Point Field::at(const Numeric alt, const Numeric lat, const Numeric lon) const
    try {
  ARTS_USER_ERROR_IF(
      alt > top_of_atmosphere,
      "Cannot get values above the top of the atmosphere, which is at: {}"
      " m.\nYour max input altitude is: {} m.",
      top_of_atmosphere,
      alt)

  Point out;
  for (auto &&key : keys()) out[key] = operator[](key).at(alt, lat, lon);
  out.check_and_fix();
  return out;
}
ARTS_METHOD_ERROR_CATCH

Point Field::at(const Vector3 pos) const try {
  return at(pos[0], pos[1], pos[2]);
}
ARTS_METHOD_ERROR_CATCH
}  // namespace Atm

std::string std::formatter<AtmKeyVal>::to_string(const AtmKeyVal &v) const {
  std::string out;
  return std::visit(
      [fmt = tags.get_format_args()](const auto &val) {
        return std::vformat(fmt.c_str(), std::make_format_args(val));
      },
      v);
}

bool Atm::pressure_is_increasing(const std::span<const Point> &atm) {
  ARTS_USER_ERROR_IF(atm.size() < 2,
                     "Must have at least two atmospheric points")

  const auto begin = atm.begin();
  const auto end   = atm.end();

  const bool pressure_increasing = begin->pressure < next(begin)->pressure;

  if (pressure_increasing) {
    for (auto it = begin; it != end - 1; it = next(it)) {
      ARTS_USER_ERROR_IF(it->pressure >= next(it)->pressure,
                         "Pressure is not consistently increasing")
    }
  } else {
    for (auto it = begin; it != end - 1; it = next(it)) {
      ARTS_USER_ERROR_IF(not(it->pressure >= next(it)->pressure),
                         "Pressure is not consistently decreasing")
    }
  }

  return pressure_increasing;
}

AtmField Atm::atm_from_profile(const std::span<const Point> &atm,
                               const AscendingGrid &altitudes,
                               const InterpolationExtrapolation &alt_extrap,
                               const Numeric &top_of_atmosphere) {
  using std::distance;
  using std::next;

  const Size N = altitudes.size();

  ARTS_USER_ERROR_IF(N != atm.size(),
                     "Inconsistent altitude grid size: {} != {}",
                     N,
                     atm.size())

  SortedGriddedField3 gf3{};
  gf3.grid<0>()     = altitudes;
  gf3.grid<1>()     = Vector{0.0};
  gf3.grid<2>()     = Vector{0.0};
  gf3.gridname<0>() = "Altitude";
  gf3.gridname<1>() = "Latitude";
  gf3.gridname<2>() = "Longitude";
  gf3.data.resize(N, 1, 1);

  Atm::Data data;
  data.alt_low = alt_extrap;
  data.alt_upp = alt_extrap;
  data.lat_low = InterpolationExtrapolation::Nearest;
  data.lat_upp = InterpolationExtrapolation::Nearest;
  data.lon_low = InterpolationExtrapolation::Nearest;
  data.lon_upp = InterpolationExtrapolation::Nearest;

  AtmField out;

  if (pressure_is_increasing(atm)) {
    for (auto &key : atm.front().keys()) {
      for (Size i = 0; i < atm.size(); i++) {
        // Add from ending
        gf3.data[N - 1 - i, 0, 0] = atm[i][key];
      }

      gf3.data_name = std::format("{}", key);
      data.data     = gf3;
      out[key]      = data;
    }
  } else {
    for (auto &key : atm.front().keys()) {
      for (Size i = 0; i < atm.size(); i++) {
        // Add from beginning
        gf3.data[i, 0, 0] = atm[i][key];
      }

      gf3.data_name = std::format("{}", key);
      data.data     = gf3;
      out[key]      = data;
    }
  }

  if (std::isnan(top_of_atmosphere)) {
    out.top_of_atmosphere = altitudes.back();
  } else {
    out.top_of_atmosphere = top_of_atmosphere;
  }

  return out;
}

void Atm::extend_in_pressure(
    ArrayOfAtmPoint &atm,
    const Numeric &new_pressure,
    const InterpolationExtrapolation extrapolation_type,
    const bool logarithmic) {
  const bool pressure_increasing = pressure_is_increasing(atm);

  std::span<const AtmPoint> prof;
  std::array<Numeric, 2> bounds;
  ArrayOfAtmPoint::const_iterator pos;

  if (pressure_increasing) {
    auto up =
        std::ranges::lower_bound(atm, new_pressure, {}, &AtmPoint::pressure);
    auto lo = up - 1;

    if (lo <= atm.begin()) {
      lo = atm.begin();
      up = atm.begin() + 1;
    } else if (up >= atm.end()) {
      lo = atm.end() - 2;
      up = atm.end() - 1;
    }

    prof   = {lo, up + 1};
    bounds = {lo->pressure, up->pressure};

    if (new_pressure < bounds[0]) {
      pos = lo;
    } else if (new_pressure > bounds[1]) {
      pos = up + 1;
    } else {
      pos = up;
    }
  } else {
    auto lo =
        (std::ranges::lower_bound(
             atm | std::views::reverse, new_pressure, {}, &AtmPoint::pressure) +
         1)
            .base();
    auto up = lo - 1;

    if (lo >= atm.end()) {
      lo = atm.end() - 1;
      up = atm.end() - 2;
    } else if (up <= atm.begin()) {
      lo = atm.begin() + 1;
      up = atm.begin();
    }

    prof   = {up, lo + 1};
    bounds = {lo->pressure, up->pressure};

    if (new_pressure < bounds[0]) {
      pos = lo + 1;
    } else if (new_pressure > bounds[1]) {
      pos = up;
    } else {
      pos = lo;
    }
  }

  if (std::ranges::any_of(prof, Cmp::eq(new_pressure), &AtmPoint::pressure))
    return;

  const AscendingGrid altitudes =
      logarithmic ? AscendingGrid{std::log(bounds[0]), std::log(bounds[1])}
                  : AscendingGrid{bounds[0], bounds[1]};
  const Numeric palt = logarithmic ? std::log(new_pressure) : new_pressure;
  const Numeric top_of_atmosphere = std::max(palt, altitudes.back());

  //! Use a fake 1D atmosphere to extend all the data
  auto p = atm.insert(
      pos,
      atm_from_profile(prof, altitudes, extrapolation_type, top_of_atmosphere)
          .at(palt, 0, 0));
  p->pressure = new_pressure;
}

AtmField AtmField::gridded(const AscendingGrid &alt,
                           const AscendingGrid &lat,
                           const AscendingGrid &lon) const {
  AtmField out{IsoRatioOption::None};
  out.top_of_atmosphere = top_of_atmosphere;

  const Index nalt = alt.size();
  const Index nlat = lat.size();
  const Index nlon = lon.size();

  ARTS_USER_ERROR_IF(
      lagrange_interp::loncross::cycle(lon.front()) != lon.front() or
          lagrange_interp::loncross::cycle(lon.back()) != lon.back(),
      "The longitude grid is incorrect.  It needs to be [-180, 180), its: {:Bs,}",
      lon);

  ARTS_USER_ERROR_IF(
      lat.front() < -90 or lat.back() > 90,
      "The latitude grid is incorrect.  It needs to be [-90, 90], its: {:Bs,}",
      lat);

  ARTS_USER_ERROR_IF(nalt * nlat * nlon == 0, "Must have a grid")

  ARTS_USER_ERROR_IF(alt.back() > top_of_atmosphere or lat.front() < -90 or
                         lat.back() > 90 or lon.front() < -180 or
                         lon.back() > 180,
                     R"(Error gridding the atmospheric field.

Top of atmosphere: {}
Altitude grid:     {:B,}
Latitude grid:     {:B,}
Longitude grid:    {:B,}

Altitude grid must be within the top of the atmosphere.
Latitude must be within -90 to 90 degrees.
Longitude must be within -180 to 180 degrees.
)",
                     top_of_atmosphere,
                     alt,
                     lat,
                     lon)

  SortedGriddedField3 gf3{
      .data_name  = "",
      .data       = Tensor3(nalt, nlat, nlon),
      .grid_names = {"Altitude", "Latitude", "Longitude"},
      .grids      = {alt, lat, lon},
  };

  for (auto key : keys()) {
    const auto &atm_data = operator[](key);

    gf3.data_name = std::format("{}", key);

    std::string error;

#pragma omp parallel for collapse(3)
    for (Index i = 0; i < nalt; i++) {
      for (Index j = 0; j < nlat; j++) {
        for (Index k = 0; k < nlon; k++) {
          try {
            gf3[i, j, k] = atm_data.at(alt[i], lat[j], lon[k]);
          } catch (std::exception &e) {
#pragma omp critical
            if (error.empty()) error = e.what();
          }
        }
      }
    }

    ARTS_USER_ERROR_IF(not error.empty(), "Error for key {}:\n\n{}", key, error)

    auto &val = out[key];
    val.data  = gf3;
    val.adjust_interpolation_extrapolation();
  }

  return out;
}

Atm::Xarr::Xarr(const AtmField &atm, std::vector<Atm::Field::KeyVal> keys_)
    : toa(atm.top_of_atmosphere), keys(std::move(keys_)) {
  if (keys.empty()) keys = atm.keys();

  if (keys.empty()) {
    altitudes  = AscendingGrid{};
    latitudes  = AscendingGrid{};
    longitudes = AscendingGrid{};
    data       = Tensor4{};
    return;
  }

  {
    const auto &key      = keys.front();
    const auto &atm_data = atm[key].data;
    ARTS_USER_ERROR_IF(
        not std::holds_alternative<SortedGriddedField3>(atm_data),
        "Data for key {} is not a SortedGriddedField3",
        key)
    const auto &gf3 = std::get<SortedGriddedField3>(atm_data);

    altitudes  = gf3.grid<0>();
    latitudes  = gf3.grid<1>();
    longitudes = gf3.grid<2>();
    data.resize(
        keys.size(), altitudes.size(), latitudes.size(), longitudes.size());

    data[0] = gf3.data;
  }

  for (Size i = 1; i < keys.size(); i++) {
    const auto &key      = keys[i];
    const auto &atm_data = atm[key].data;
    ARTS_USER_ERROR_IF(
        not std::holds_alternative<SortedGriddedField3>(atm_data),
        "Data for key {} is not a SortedGriddedField3",
        key)
    const auto &gf3 = std::get<SortedGriddedField3>(atm_data);
    ARTS_USER_ERROR_IF(
        gf3.grid<0>() != altitudes.vec() or gf3.grid<1>() != latitudes.vec() or
            gf3.grid<2>() != longitudes.vec(),
        R"(Grids for key {0} is not the same as the first key (first key: {1})

Grids for {0}:
  Altitude:  {2:B,}
  Latitude:  {3:B,}
  Longitude: {4:B,}

Grids for {1}:
  Altitude:  {5:B,}
  Latitude:  {6:B,}
  Longitude: {7:B,}
)",
        key,
        keys.front(),
        altitudes,
        latitudes,
        longitudes,
        gf3.grid<0>(),
        gf3.grid<1>(),
        gf3.grid<2>())

    data[i] = gf3.data;
  }
}

std::string std::formatter<AtmPoint>::to_string(const AtmPoint &v) const {
  const std::string_view sep = tags.sep(true);

  std::string out = tags.vformat(R"("pressure": )"sv,
                                 v.pressure,
                                 sep,
                                 R"("temperature": )"sv,
                                 v.temperature,
                                 sep,
                                 R"("mag" :)"sv,
                                 v.mag,
                                 sep,
                                 R"("wind": )"sv,
                                 v.wind,
                                 sep);

  if (tags.short_str) {
    out += tags.vformat(R"("SpeciesEnum": )"sv,
                        v.specs.size(),
                        sep,
                        R"("SpeciesIsotope": )"sv,
                        v.isots.size(),
                        sep,
                        R"("QuantumLevelIdentifier": )"sv,
                        v.nlte.size(),
                        sep,
                        R"("ScatteringSpeciesProperty": )"sv,
                        v.ssprops.size());

  } else {
    out += tags.vformat(R"("SpeciesEnum": )"sv,
                        v.specs,
                        sep,
                        R"("SpeciesIsotope": )"sv,
                        v.isots,
                        sep,
                        R"("QuantumLevelIdentifier": )"sv,
                        v.nlte,
                        sep,
                        R"("ScatteringSpeciesProperty": )"sv,
                        v.ssprops);
  }

  return tags.bracket ? ("{" + out + "}") : out;
}

std::string std::formatter<AtmField>::to_string(const AtmField &v) const {
  std::string out;

  if (tags.short_str) {
    const std::string_view sep = tags.sep();

    out = tags.vformat(R"("top_of_atmosphere": )"sv,
                       v.top_of_atmosphere,
                       sep,
                       R"("Base": )"sv,
                       v.other.size(),
                       sep,
                       R"("SpeciesEnum": )"sv,
                       v.specs.size(),
                       sep,
                       R"("SpeciesIsotope": )"sv,
                       v.isots.size(),
                       sep,
                       R"("QuantumLevelIdentifier": )"sv,
                       v.nlte.size(),
                       sep,
                       R"("ScatteringSpeciesProperty": )"sv,
                       v.ssprops.size());
  } else {
    const std::string_view sep = tags.sep(true);

    out = tags.vformat(R"("top_of_atmosphere": )"sv,
                       v.top_of_atmosphere,
                       sep,
                       R"("Base": )"sv,
                       v.other,
                       sep,
                       R"("SpeciesEnum": )"sv,
                       v.specs,
                       sep,
                       R"("SpeciesIsotope": )"sv,
                       v.isots,
                       sep,
                       R"("QuantumLevelIdentifier": )"sv,
                       v.nlte,
                       sep,
                       R"("ScatteringSpeciesProperty": )"sv,
                       v.ssprops);
  }

  return tags.bracket ? ("{" + out + "}") : out;
}

static_assert(
    std::same_as<typename AtmField::KeyVal, AtmKeyVal>,
    "The order of arguments in the template of which Field inherits from is "
    "wrong.  KeyVal must be defined in the same way for this to work.");
