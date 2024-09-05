#include "atm.h"

#include <matpack.h>

#include <algorithm>
#include <iomanip>
#include <limits>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <type_traits>
#include <variant>
#include <vector>

#include "compare.h"
#include "configtypes.h"
#include "debug.h"
#include "enumsFieldComponent.h"
#include "hitran_species.h"
#include "interp.h"
#include "isotopologues.h"

AtmKey to_wind(const String& x) {
  switch(to<FieldComponent>(x)) {
    case FieldComponent::u:
      return AtmKey::wind_u;
    case FieldComponent::v:
      return AtmKey::wind_v;
    case FieldComponent::w:
      return AtmKey::wind_w;
  }
}

AtmKey to_mag(const String& x) {
  switch(to<FieldComponent>(x)) {
    case FieldComponent::u:
      return AtmKey::mag_u;
    case FieldComponent::v:
      return AtmKey::mag_v;
    case FieldComponent::w:
      return AtmKey::mag_w;
  }
}

Numeric AtmPoint::operator[](SpeciesEnum x) const try {
  return specs.at(x);
} catch (std::out_of_range &) {
  ARTS_USER_ERROR("Species VMR not found: \"", toString<1>(x), '"')
}

Numeric AtmPoint::operator[](const SpeciesIsotope &x) const try {
  return isots.at(x);
} catch (std::out_of_range &) {
  ARTS_USER_ERROR("Isotopologue ratio not found: \"", x, '"')
}

Numeric AtmPoint::operator[](const QuantumIdentifier &x) const try {
  return nlte.at(x);
} catch (std::out_of_range &) {
  ARTS_USER_ERROR("QuantumIdentifier not found: \"", x, '"')
}

Numeric AtmPoint::operator[](const ParticulatePropertyTag &x) const try {
  return partp.at(x);
} catch (std::out_of_range &) {
  ARTS_USER_ERROR("ParticulatePropertyTag not found: \"", x, '"')
}

Numeric AtmPoint::operator[](AtmKey x) const {
  switch (x) {
    case AtmKey::t:
      return temperature;
    case AtmKey::p:
      return pressure;
    case AtmKey::wind_u:
      return wind[0];
    case AtmKey::wind_v:
      return wind[1];
    case AtmKey::wind_w:
      return wind[2];
    case AtmKey::mag_u:
      return mag[0];
    case AtmKey::mag_v:
      return mag[1];
    case AtmKey::mag_w:
      return mag[2];
  }
  ARTS_USER_ERROR("Cannot reach")
}

Numeric &AtmPoint::operator[](SpeciesEnum x) { return specs[x]; }

Numeric &AtmPoint::operator[](const SpeciesIsotope &x) { return isots[x]; }

Numeric &AtmPoint::operator[](const QuantumIdentifier &x) { return nlte[x]; }

Numeric &AtmPoint::operator[](const ParticulatePropertyTag &x) {
  return partp[x];
}

Numeric &AtmPoint::operator[](AtmKey x) {
  switch (x) {
    case AtmKey::t:
      return temperature;
    case AtmKey::p:
      return pressure;
    case AtmKey::wind_u:
      return wind[0];
    case AtmKey::wind_v:
      return wind[1];
    case AtmKey::wind_w:
      return wind[2];
    case AtmKey::mag_u:
      return mag[0];
    case AtmKey::mag_v:
      return mag[1];
    case AtmKey::mag_w:
      return mag[2];
  }
  ARTS_USER_ERROR("Cannot reach")
}

std::pair<Numeric, Numeric> AtmPoint::levels(
    const QuantumIdentifier &band) const {
  return {operator[](band.LowerLevel()), operator[](band.UpperLevel())};
}

std::ostream &operator<<(std::ostream &os, const AtmKeyVal &key) {
  std::visit([&os](const auto &k) { os << k; }, key);
  return os;
}

template <>
const Atm::Data &
FieldMap::Map<Atm::Data,
              AtmKey,
              SpeciesEnum,
              SpeciesIsotope,
              QuantumIdentifier,
              ParticulatePropertyTag>::operator[](const KeyVal &k) const try {
  return std::visit(
      [this](auto &key) -> const Atm::Data & {
        return this->map<decltype(key)>().at(key);
      },
      k);
} catch (std::out_of_range &) {
  throw std::out_of_range(var_string("Key not found in map: \"", k, '\"'));
} catch (...) {
  throw;
}

template <>
Atm::Data &
FieldMap::Map<Atm::Data,
              AtmKey,
              SpeciesEnum,
              SpeciesIsotope,
              QuantumIdentifier,
              ParticulatePropertyTag>::operator[](const KeyVal &k) try {
  return std::visit(
      [this](auto &key) -> Atm::Data & {
        return const_cast<Map *>(this)->map<decltype(key)>()[key];
      },
      k);
}
ARTS_METHOD_ERROR_CATCH

template <>
bool FieldMap::Map<Atm::Data,
                   AtmKey,
                   SpeciesEnum,
                   SpeciesIsotope,
                   QuantumIdentifier,
                   ParticulatePropertyTag>::contains(const KeyVal &key) const {
  return std::visit(
      [this](auto &k) -> bool { return this->map<decltype(k)>().contains(k); },
      key);
}

namespace Atm {
Point::Point(const IsoRatioOption isots_key) {
  switch (isots_key) {
    case IsoRatioOption::Builtin: {
      const SpeciesIsotopologueRatios x =
          Species::isotopologue_ratiosInitFromBuiltin();
      for (Index i = 0; i < x.maxsize; i++) {
        if (Species::Isotopologues[i].joker()) continue;
        if (Species::is_predefined_model(Species::Isotopologues[i])) continue;
        isots[Species::Isotopologues[i]] = x.data[i];
      }
    } break;
    case IsoRatioOption::Hitran: {
      const SpeciesIsotopologueRatios x = Hitran::isotopologue_ratios();
      for (Index i = 0; i < x.maxsize; i++) {
        if (Species::Isotopologues[i].joker()) continue;
        if (Species::is_predefined_model(Species::Isotopologues[i])) continue;
        isots[Species::Isotopologues[i]] = x.data[i];
      }
    } break;
    case IsoRatioOption::None:
    default:
      break;
  }
}

Field::Field(const IsoRatioOption isots_key) {
  switch (isots_key) {
    case IsoRatioOption::Builtin: {
      const SpeciesIsotopologueRatios x =
          Species::isotopologue_ratiosInitFromBuiltin();
      for (Index i = 0; i < x.maxsize; i++) {
        if (Species::Isotopologues[i].joker()) continue;
        if (Species::is_predefined_model(Species::Isotopologues[i])) continue;
        isots()[Species::Isotopologues[i]] = x.data[i];
      }
    } break;
    case IsoRatioOption::Hitran: {
      const SpeciesIsotopologueRatios x = Hitran::isotopologue_ratios();
      for (Index i = 0; i < x.maxsize; i++) {
        if (Species::Isotopologues[i].joker()) continue;
        if (Species::is_predefined_model(Species::Isotopologues[i])) continue;
        isots()[Species::Isotopologues[i]] = x.data[i];
      }
    } break;
    case IsoRatioOption::None:
    default:
      break;
  }
}

const std::unordered_map<QuantumIdentifier, Data> &Field::nlte() const {
  return map<QuantumIdentifier>();
}

const std::unordered_map<SpeciesEnum, Data> &Field::specs() const {
  return map<SpeciesEnum>();
}

const std::unordered_map<SpeciesIsotope, Data> &Field::isots() const {
  return map<SpeciesIsotope>();
}

const std::unordered_map<AtmKey, Data> &Field::other() const {
  return map<AtmKey>();
}

const std::unordered_map<ParticulatePropertyTag, Data> &Field::partp() const {
  return map<ParticulatePropertyTag>();
}

std::unordered_map<QuantumIdentifier, Data> &Field::nlte() {
  return map<QuantumIdentifier>();
}

std::unordered_map<SpeciesEnum, Data> &Field::specs() {
  return map<SpeciesEnum>();
}

std::unordered_map<SpeciesIsotope, Data> &Field::isots() {
  return map<SpeciesIsotope>();
}

std::unordered_map<AtmKey, Data> &Field::other() { return map<AtmKey>(); }

std::unordered_map<ParticulatePropertyTag, Data> &Field::partp() {
  return map<ParticulatePropertyTag>();
}

std::ostream &operator<<(std::ostream &os, const Point &atm) {
  os << "Temperature: " << atm.temperature << " K,\n";
  os << "Pressure: " << atm.pressure << " Pa,\n";
  os << "Wind Field: [u: " << atm.wind[0] << ", v: " << atm.wind[1]
     << ", w: " << atm.wind[2] << "] m/s,\n";
  os << "Magnetic Field: [u: " << atm.mag[0] << ", v: " << atm.mag[1]
     << ", w: " << atm.mag[2] << "] T";

  for (auto &spec : atm.specs) {
    os << ",\n" << toString<1>(spec.first) << ": " << spec.second;
  }
  for (auto &spec : atm.isots) {
    os << ",\n" << spec.first << ": " << spec.second;
  }

  for (auto &vals : atm.nlte) {
    os << ",\n" << vals.first << ": " << vals.second;
  }

  return os;
}

std::ostream &operator<<(std::ostream &os, const Field &atm) {
  const auto printer = [&](auto &data) {
    using T = decltype(data);
    if constexpr (isFunctionalDataType<T>)
      os << " Functional";
    else if constexpr (isNumeric<T>)
      os << ' ' << data;
    else
      os << '\n' << data;
  };

  std::string_view space = "";
  for (auto &&key : atm.keys()) {
    os << std::exchange(space, ",\n") << key << ":";
    std::visit(printer, atm[key].data);
  }

  return os;
}

Numeric Point::mean_mass(SpeciesEnum s) const {
  Numeric ratio = 0.0;
  Numeric mass  = 0.0;
  for (auto &[isot, this_ratio] : isots) {
    if (isot.spec == s and not(is_predefined_model(isot) or isot.joker())) {
      ratio += this_ratio;
      mass  += this_ratio * isot.mass;
    }
  }

  ARTS_USER_ERROR_IF(ratio == 0,
                     "Cannot find a ratio for the mean mass of species \"",
                     toString<1>(s),
                     '"')

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

std::vector<KeyVal> Point::keys() const {
  std::vector<KeyVal> out;
  out.reserve(size());
  for (auto &a : enumtyps::AtmKeyTypes) out.emplace_back(a);
  for (auto &a : specs) out.emplace_back(a.first);
  for (auto &a : nlte) out.emplace_back(a.first);
  for (auto &a : partp) out.emplace_back(a.first);
  for (auto &a : isots) out.emplace_back(a.first);
  return out;
}

Index Point::nspec() const { return static_cast<Index>(specs.size()); }

Index Point::npart() const { return static_cast<Index>(partp.size()); }

Index Point::nisot() const { return static_cast<Index>(isots.size()); }

Index Point::nnlte() const { return static_cast<Index>(nlte.size()); }

Index Point::size() const {
  return nspec() + nnlte() + nother() + npart() + nisot();
}

Index Field::nspec() const { return static_cast<Index>(specs().size()); }

Index Field::nisot() const { return static_cast<Index>(isots().size()); }

Index Field::npart() const { return static_cast<Index>(partp().size()); }

Index Field::nnlte() const { return static_cast<Index>(nlte().size()); }

Index Field::nother() const { return static_cast<Index>(other().size()); }

String Data::data_type() const {
  if (std::holds_alternative<GriddedField3>(data)) return "GriddedField3";
  if (std::holds_alternative<Numeric>(data)) return "Numeric";
  if (std::holds_alternative<FunctionalData>(data)) return "FunctionalData";
  ARTS_ASSERT(
      false,
      "Cannot be reached, you have added a new type but not doen the plumbing...")
  ARTS_USER_ERROR("Cannot understand data type; is this a new type")
}

namespace detail {
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

Limits find_limits(const GriddedField3 &gf3) {
  return {gf3.grid<0>().front(),
          gf3.grid<0>().back(),
          gf3.grid<1>().front(),
          gf3.grid<1>().back(),
          gf3.grid<2>().front(),
          gf3.grid<2>().back()};
}

Vector vec_interp(const Numeric &v,
                  const Vector &alt,
                  const Vector &,
                  const Vector &) {
  Vector out(alt.size(), v);
  return out;
}

Vector vec_interp(const FunctionalData &v,
                  const Vector &alt,
                  const Vector &lat,
                  const Vector &lon) {
  const Index n = alt.size();
  Vector out(n);
  for (Index i = 0; i < n; i++) out[i] = v(alt[i], lat[i], lon[i]);
  return out;
}

template <Index poly_alt,
          Index poly_lat,
          Index poly_lon,
          bool precompute = false>
struct interp_helper {
  using AltLag = my_interp::Lagrange<poly_alt>;
  using LatLag = my_interp::Lagrange<poly_lat>;
  using LonLag =
      std::conditional_t<poly_lon == 0,
                         my_interp::Lagrange<0>,
                         my_interp::Lagrange<poly_lon,
                                             false,
                                             GridType::Cyclic,
                                             my_interp::cycle_m180_p180>>;

  Array<AltLag> lags_alt;
  Array<LatLag> lags_lat;
  Array<LonLag> lags_lon;
  [[no_unique_address]] std::conditional_t<
      precompute,
      decltype(my_interp::flat_interpweights(lags_alt, lags_lat, lags_lon)),
      my_interp::Empty> iws{};

  constexpr interp_helper(const Vector &alt_grid,
                          const Vector &lat_grid,
                          const Vector &lon_grid,
                          const Vector &alt,
                          const Vector &lat,
                          const Vector &lon)
    requires(not precompute)
      : lags_alt(
            my_interp::lagrange_interpolation_list<AltLag>(alt, alt_grid, -1)),
        lags_lat(
            my_interp::lagrange_interpolation_list<LatLag>(lat, lat_grid, -1)),
        lags_lon(my_interp::lagrange_interpolation_list<LonLag>(
            lon, lon_grid, -1)) {}

  constexpr interp_helper(const Vector &alt_grid,
                          const Vector &lat_grid,
                          const Vector &lon_grid,
                          const Vector &alt,
                          const Vector &lat,
                          const Vector &lon)
    requires(precompute)
      : lags_alt(
            my_interp::lagrange_interpolation_list<AltLag>(alt, alt_grid, -1)),
        lags_lat(
            my_interp::lagrange_interpolation_list<LatLag>(lat, lat_grid, -1)),
        lags_lon(
            my_interp::lagrange_interpolation_list<LonLag>(lon, lon_grid, -1)),
        iws(flat_interpweights(lags_alt, lags_lat, lags_lon)) {}

  Vector operator()(const Tensor3 &data) const {
    if constexpr (precompute)
      return flat_interp(data, iws, lags_alt, lags_lat, lags_lon);
    else
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

  const Index n = lat_grid.size();
  const Index m = lon_grid.size();

  for (Index i = 0; i < alt.size(); i++) {
    const Index alt0 = interpolater.lags_alt[i].pos * m * n;
    const Index lat0 = interpolater.lags_lat[i].pos * m;
    const Index lon0 = interpolater.lags_lon[i].pos;

    Index j = 0;
    for (Index idx0 = 0; idx0 < 1 + poly_alt; idx0++) {
      for (Index idx1 = 0; idx1 < 1 + poly_lat; idx1++) {
        for (Index idx2 = 0; idx2 < 1 + poly_lon; idx2++) {
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

Vector vec_interp(const GriddedField3 &v,
                  const Vector &alt,
                  const Vector &lat,
                  const Vector &lon) {
  ARTS_ASSERT(v.shape()[0] > 0)
  ARTS_ASSERT(v.shape()[1] > 0)
  ARTS_ASSERT(v.shape()[2] > 0)

  const bool d1 = v.shape()[0] == 1;
  const bool d2 = v.shape()[1] == 1;
  const bool d3 = v.shape()[2] == 1;

  const Index n = alt.size();

  if (d1 and d2 and d3) return Vector(n, v.data(0, 0, 0));
  if (d1 and d2)
    return tvec_interp<0, 0, 1>(
        v.data, v.grid<0>(), v.grid<1>(), v.grid<2>(), alt, lat, lon);
  if (d1 and d3)
    return tvec_interp<0, 1, 0>(
        v.data, v.grid<0>(), v.grid<1>(), v.grid<2>(), alt, lat, lon);
  if (d2 and d3)
    return tvec_interp<1, 0, 0>(
        v.data, v.grid<0>(), v.grid<1>(), v.grid<2>(), alt, lat, lon);
  if (d1)
    return tvec_interp<0, 1, 1>(
        v.data, v.grid<0>(), v.grid<1>(), v.grid<2>(), alt, lat, lon);
  if (d2)
    return tvec_interp<1, 0, 1>(
        v.data, v.grid<0>(), v.grid<1>(), v.grid<2>(), alt, lat, lon);
  if (d3)
    return tvec_interp<1, 1, 0>(
        v.data, v.grid<0>(), v.grid<1>(), v.grid<2>(), alt, lat, lon);
  return tvec_interp<1, 1, 1>(
      v.data, v.grid<0>(), v.grid<1>(), v.grid<2>(), alt, lat, lon);
}

Numeric limit(const Data &data, ComputeLimit lim, Numeric orig) {
  ARTS_USER_ERROR_IF(lim.type == InterpolationExtrapolation::None,
                     "Altitude limit breaced.  Position (",
                     lim.alt,
                     ", ",
                     lim.lat,
                     ", ",
                     lim.lon,
                     ") is out-of-bounds when no extrapolation is wanted")

  if (lim.type == InterpolationExtrapolation::Zero) return 0;

  if (lim.type == InterpolationExtrapolation::Nearest)
    return std::visit(
        [&](auto &d) {
          return vec_interp(d, {lim.alt}, {lim.lat}, {lim.lon})[0];
        },
        data.data);

  return orig;
}

constexpr InterpolationExtrapolation combine(InterpolationExtrapolation a,
                                             InterpolationExtrapolation b) {
  using enum InterpolationExtrapolation;
  switch (a) {
    case None:
      return None;
    case Zero: {
      switch (b) {
        case None:
          return None;
        case Zero:
          return Zero;
        case Nearest:
          return Zero;
        case Linear:
          return Zero;
      }
      std::unreachable();
    }
    case Nearest: {
      switch (b) {
        case None:
          return None;
        case Zero:
          return Zero;
        case Nearest:
          return Nearest;
        case Linear:
          return Nearest;
      }
      std::unreachable();
    }
    case Linear:
      return b;
  }

  return a;
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

Vector vec_interp(const Data &data,
                  const Vector &alt,
                  const Vector &lat,
                  const Vector &lon) {
  const auto compute = [&](auto &d) { return vec_interp(d, alt, lat, lon); };

  // Perform the interpolation
  Vector out = std::visit(compute, data.data);

  // Fix the extrapolations for ZERO and NONE and NEAREST
  const auto lim =
      std::visit([](auto &d) { return find_limits(d); }, data.data);
  const Index n = alt.size();
  for (Index i = 0; i < n; i++) {
    out[i] = limit(data, find_limit(data, lim, alt[i], lat[i], lon[i]), out[i]);
  }

  return out;
}
}  // namespace detail

bool Point::is_lte() const noexcept { return nlte.empty(); }

template <class Key, class T, class Hash, class KeyEqual, class Allocator>
std::vector<AtmKey> get_keys(
    const std::unordered_map<Key, T, Hash, KeyEqual, Allocator> &map) {
  std::vector<Key> out(map.size());
  std::transform(
      map.begin(), map.end(), out.begin(), [](auto &v) { return v.first; });
  return out;
}

ArrayOfQuantumIdentifier Field::nlte_keys() const {
  return keys<QuantumIdentifier>();
}

void Data::rescale(Numeric x) {
  std::visit(
      [x](auto &v) {
        using T = decltype(v);
        if constexpr (isFunctionalDataType<T>) {
          v = FunctionalData{[x, f = v](Numeric alt, Numeric lat, Numeric lon) -> Numeric {
            return x * f(alt, lat, lon);
          }};
        } else if constexpr (isGriddedField3<T>) {
          v.data *= x;
        } else {
          v *= x;
        }
      },
      data);
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
        "Wind field [wind_u wind_v wind_w] is: ",
        wind)
  }

  if (std::ranges::all_of(mag, [](auto v) { return nonstd::isnan(v); })) {
    mag = {0., 0., 0.};
  } else {
    ARTS_USER_ERROR_IF(
        std::ranges::any_of(mag, [](auto v) { return nonstd::isnan(v); }),
        "Cannot have partially missing magnetic field.  Consider setting the missing field to zero or add it completely.\n"
        "Magnetic field [mag_u mag_v mag_w] is: ",
        mag)
  }

  for (auto &spec : specs) {
    ARTS_USER_ERROR_IF(nonstd::isnan(spec.second) or spec.second < 0.0,
                       "VMR for \"",
                       toString<1>(spec.first),
                       "\" is ",
                       spec.second)
  }

  for (auto &isot : isots) {
    //! Cannot check isnan because it is a valid state for isotopologue ratios
    ARTS_USER_ERROR_IF(isot.second < 0.0,
                       "Isotopologue ratio for \"",
                       isot.first.FullName(),
                       "\" is ",
                       isot.second)
  }

  for (auto &nl : nlte) {
    ARTS_USER_ERROR_IF(nonstd::isnan(nl.second) or nl.second < 0.0,
                       "Non-LTE ratio for \"",
                       nl.first,
                       "\" is ",
                       nl.second)
  }

  for (auto &pp : partp) {
    ARTS_USER_ERROR_IF(nonstd::isnan(pp.second),
                       "Particulate Property Tag value for \"",
                       pp.first,
                       "\" is ",
                       pp.second)
  }
}
ARTS_METHOD_ERROR_CATCH

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

std::ostream &operator<<(std::ostream &os, const Array<Point> &a) {
  for (auto &x : a) os << x << '\n';
  return os;
}

ExhaustiveConstVectorView Data::flat_view() const {
  return std::visit(
      [](auto &X) -> ExhaustiveConstVectorView {
        using T = std::remove_cvref_t<decltype(X)>;
        if constexpr (std::same_as<T, GriddedField3>)
          return X.data.flat_view();
        else if constexpr (std::same_as<T, Numeric>)
          return ExhaustiveConstVectorView{X};
        else if constexpr (std::same_as<T, FunctionalData>)
          return ExhaustiveConstVectorView{};
        else
          static_assert(
              RawDataType<T>,
              "Cannot be reached, you have added a new type but not done the plumbing...");
      },
      data);
}

ExhaustiveVectorView Data::flat_view() {
  return std::visit(
      [](auto &X) -> ExhaustiveVectorView {
        using T = std::remove_cvref_t<decltype(X)>;
        if constexpr (std::same_as<T, GriddedField3>)
          return X.data.flat_view();
        else if constexpr (std::same_as<T, Numeric>)
          return ExhaustiveVectorView{X};
        else if constexpr (std::same_as<T, FunctionalData>)
          return ExhaustiveVectorView{};
        else
          static_assert(
              RawDataType<T>,
              "Cannot be reached, you have added a new type but not done the plumbing...");
      },
      data);
}

namespace interp {
using altlag1 = my_interp::Lagrange<1>;
using altlag0 = my_interp::Lagrange<0>;

using latlag1 = my_interp::Lagrange<1>;
using latlag0 = my_interp::Lagrange<0>;

using lonlag1 =
    my_interp::Lagrange<1, false, GridType::Cyclic, my_interp::cycle_m180_p180>;
using lonlag0 =
    my_interp::Lagrange<0, false, GridType::Cyclic, my_interp::cycle_m180_p180>;

using altlags = std::variant<altlag0, altlag1>;
using latlags = std::variant<latlag0, latlag1>;
using lonlags = std::variant<lonlag0, lonlag1>;

std::array<std::pair<Index, Numeric>, 8> flat_weight_(const GriddedField3 &gf3,
                                                      const Numeric alt,
                                                      const Numeric lat,
                                                      const Numeric lon) {
  if (not gf3.ok()) throw std::runtime_error("bad field");

  const Index nalt = gf3.grid<0>().size();
  const Index nlat = gf3.grid<1>().size();
  const Index nlon = gf3.grid<2>().size();

  return std::visit(
      [NN = nlat * nlon, N = nlon]<typename ALT, typename LAT, typename LON>(
          const ALT &al, const LAT &la, const LON &lo) {
        const auto x = interpweights(al, la, lo);
        constexpr std::pair<Index, Numeric> v0{0, 0.0};
        std::array<std::pair<Index, Numeric>, 8> out{
            v0, v0, v0, v0, v0, v0, v0, v0};

        Index m = 0;
        for (Index i = 0; i < al.size(); i++) {
          for (Index j = 0; j < la.size(); j++) {
            for (Index k = 0; k < lo.size(); k++, ++m) {
              out[m] = {{(al.pos + i) * NN + (la.pos + j) * N + lo.pos + k},
                        x(i, j, k)};
            }
          }
        }
        return out;
      },
      nalt == 1 ? altlags{gf3.lag<0, altlag0>(alt)}
                : altlags{gf3.lag<0, altlag1>(alt)},
      nlat == 1 ? latlags{gf3.lag<1, latlag0>(lat)}
                : latlags{gf3.lag<1, latlag1>(lat)},
      nlon == 1 ? lonlags{gf3.lag<2, lonlag0>(lon)}
                : lonlags{gf3.lag<2, lonlag1>(lon)});
}
}  // namespace interp

std::array<std::pair<Index, Numeric>, 8> flat_weight_(const GriddedField3 &data,
                                                      Numeric alt,
                                                      Numeric lat,
                                                      Numeric lon) {
  return interp::flat_weight_(data, alt, lat, lon);
}

std::array<std::pair<Index, Numeric>, 8> flat_weight_(Numeric,
                                                      Numeric,
                                                      Numeric,
                                                      Numeric) {
  constexpr std::pair<Index, Numeric> v0{0, 0.0};
  constexpr std::pair<Index, Numeric> v1{0, 1.0};
  return {v1, v0, v0, v0, v0, v0, v0, v0};
}

std::array<std::pair<Index, Numeric>, 8> flat_weight_(const FunctionalData &,
                                                      Numeric,
                                                      Numeric,
                                                      Numeric) {
  constexpr std::pair<Index, Numeric> v0{0, 0.0};
  return {v0, v0, v0, v0, v0, v0, v0, v0};
}

//! Flat weights for the positions in an atmosphere
std::array<std::pair<Index, Numeric>, 8> Data::flat_weight(
    const Numeric alt, const Numeric lat, const Numeric lon) const {
  return std::visit([&](auto &v) { return flat_weight_(v, alt, lat, lon); },
                    data);
}

std::array<std::pair<Index, Numeric>, 8> Data::flat_weight(
    const Vector3 pos) const {
  return flat_weight(pos[0], pos[1], pos[2]);
}
}  // namespace Atm

template <Atm::KeyType T>
constexpr bool cmp(const AtmKeyVal &keyval, const T &key) {
  const auto *ptr = std::get_if<T>(&keyval);
  return ptr and *ptr == key;
}

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

bool operator==(const AtmKeyVal &keyval, const SpeciesIsotope &key) {
  return cmp(keyval, key);
}

bool operator==(const SpeciesIsotope &key, const AtmKeyVal &keyval) {
  return cmp(keyval, key);
}

bool operator==(const AtmKeyVal &keyval, const QuantumIdentifier &key) {
  return cmp(keyval, key);
}

bool operator==(const QuantumIdentifier &key, const AtmKeyVal &keyval) {
  return cmp(keyval, key);
}

bool operator==(const AtmKeyVal &keyval, const ParticulatePropertyTag &key) {
  return cmp(keyval, key);
}

bool operator==(const ParticulatePropertyTag &key, const AtmKeyVal &keyval) {
  return cmp(keyval, key);
}

std::ostream &operator<<(std::ostream &os, const ParticulatePropertyTag &ppt) {
  return os << ppt.name;
}

namespace Atm {
namespace interp {
Numeric get(const GriddedField3 &gf3,
            const Numeric alt,
            const Numeric lat,
            const Numeric lon) {
  if (not gf3.ok()) throw std::runtime_error("bad field");

  return std::visit(
      [&data = gf3.data](auto &&al, auto &&la, auto &&lo) {
        return my_interp::interp(data, al, la, lo);
      },
      gf3.grid<0>().size() == 1 ? altlags{gf3.lag<0, altlag0>(alt)}
                                : altlags{gf3.lag<0, altlag1>(alt)},
      gf3.grid<1>().size() == 1 ? latlags{gf3.lag<1, latlag0>(lat)}
                                : latlags{gf3.lag<1, latlag1>(lat)},
      gf3.grid<2>().size() == 1 ? lonlags{gf3.lag<2, lonlag0>(lon)}
                                : lonlags{gf3.lag<2, lonlag1>(lon)});
}

constexpr Numeric get(const Numeric num,
                      const Numeric,
                      const Numeric,
                      const Numeric) {
  return num;
}

Numeric get(const FunctionalData &fd,
            const Numeric alt,
            const Numeric lat,
            const Numeric lon) {
  return fd(alt, lat, lon);
}

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
  const auto lim = detail::find_limit(
      data,
      std::visit([](auto &d) { return detail::find_limits(d); }, data.data),
      alt,
      lat,
      lon);

  ARTS_USER_ERROR_IF(lim.type == InterpolationExtrapolation::None,
                     "Altitude limit breaced.  Position (",
                     lim.alt,
                     ", ",
                     lim.lat,
                     ", ",
                     lim.lon,
                     ") is out-of-bounds when no extrapolation is wanted")

  if (lim.type == InterpolationExtrapolation::Zero) return 0.0;

  if (lim.type == InterpolationExtrapolation::Nearest)
    return PositionalNumeric{data.data, lim.alt, lim.lat, lim.lon};

  return std::nullopt;
}
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
      "Cannot get values above the top of the atmosphere, which is at: ",
      top_of_atmosphere,
      " m.\nYour max input altitude is: ",
      alt,
      " m.")

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
