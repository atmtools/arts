#include <enumsHydrostaticPressureOption.h>
#include <enumsIsoRatioOption.h>
#include <enumsMissingFieldComponentError.h>
#include <workspace.h>
#include <zconf.h>

#include <algorithm>
#include <iomanip>
#include <iterator>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <variant>

#include "atm.h"
#include "compare.h"
#include "configtypes.h"
#include "debug.h"
#include "igrf13.h"
#include "interp.h"
#include "interpolation.h"
#include "isotopologues.h"
#include "lbl_data.h"
#include "matpack_data.h"
#include "mc_interp.h"
#include "operators.h"
#include "predef_data.h"
#include "predefined_absorption_models.h"
#include "quantum_numbers.h"
#include "species.h"
#include "species_tags.h"
#include "xml_io.h"

void atmospheric_fieldInit(AtmField &atmospheric_field,
                           const Numeric &top_of_atmosphere,
                           const String &default_isotopologue) {
  atmospheric_field = AtmField{to<IsoRatioOption>(default_isotopologue)};
  atmospheric_field.top_of_atmosphere = top_of_atmosphere;
}

void atmospheric_pointInit(AtmPoint &atmospheric_point,
                           const String &default_isotopologue) {
  atmospheric_point = AtmPoint{to<IsoRatioOption>(default_isotopologue)};
}

template <Size I = 0, typename... T, Size N = sizeof...(T)>
std::variant<T...> read_variant(const std::variant<T...> &_,
                                const String &filename)
  requires(N != 0)
{
  try {
    typename std::tuple_element<I, std::tuple<T...>>::type x{};
    xml_read_from_file(filename, x);
    return x;
  } catch (...) {
  }

  if constexpr (I + 1 < N) {
    return read_variant<I + 1, T...>(_, filename);
  }

  ARTS_USER_ERROR("Could not read valid data from file: {}", filename)
}

template <typename T>
void append_data(
    AtmField &atmospheric_field,
    const String &basename,
    const String &extrapolation,
    const Index &missing_is_zero,
    const Index &replace_existing,
    const Index &ignore_missing,
    const std::unordered_map<T, Index> &keys,
    const auto &to_string = [](const auto &x) { return var_string(x); })
  requires(
      std::same_as<typename std::decay_t<decltype(to_string(T{}))>, String>)
{
  const String my_base = complete_basename(basename);

  Atm::Data x;
  x.alt_low = x.alt_upp = x.lat_low = x.lat_upp = x.lon_low = x.lon_upp =
      to<InterpolationExtrapolation>(extrapolation);

  for (const auto &[key, _] : keys) {
    String filename = var_string(my_base, to_string(key), ".xml");

    if (static_cast<bool>(replace_existing) or
        not atmospheric_field.contains(key)) {
      if (find_xml_file_existence(filename)) {
        x.data                 = read_variant(Atm::FieldData{1.0}, filename);
        atmospheric_field[key] = x;
      } else if (static_cast<bool>(missing_is_zero)) {
        x.data                 = 0.0;
        atmospheric_field[key] = x;
      } else if (static_cast<bool>(ignore_missing)) {
      } else {
        ARTS_USER_ERROR_IF(
            not atmospheric_field.contains(key),
            "Filename: \"{}\" does not exist"
            " and no options for workarounds are given.  Cannot populate atmospheric_field with key: {}",
            filename,
            to_string(key))
      }
    }
  }
}

void atmospheric_fieldAppendBaseData(AtmField &atmospheric_field,
                                     const String &basename,
                                     const String &extrapolation,
                                     const String &deal_with_field_component,
                                     const Index &replace_existing,
                                     const Index &allow_missing_pressure,
                                     const Index &allow_missing_temperature) {
  std::unordered_map<AtmKey, Index> keys;

  for (auto &key : enumtyps::AtmKeyTypes) {
    ++keys[key];
  }

  append_data(atmospheric_field,
              basename,
              extrapolation,
              0,
              replace_existing,
              1,
              keys,
              [](const AtmKey &x) { return var_string(x); });

  using enum AtmKey;

  ARTS_USER_ERROR_IF(
      not atmospheric_field.has(p) and
          not static_cast<bool>(allow_missing_pressure),
      "Pressure is missing from the read atmospheric field at \"{}\"",
      basename)

  ARTS_USER_ERROR_IF(
      not atmospheric_field.has(t) and
          not static_cast<bool>(allow_missing_temperature),
      "Temperature is missing from the read atmospheric field at \"{}\"",
      basename)

  switch (to<MissingFieldComponentError>(deal_with_field_component)) {
    case MissingFieldComponentError::Throw:
      if (atmospheric_field.has(wind_u) or atmospheric_field.has(wind_v) or
          atmospheric_field.has(wind_w)) {
        ARTS_USER_ERROR_IF(
            not atmospheric_field.has(wind_u, wind_v, wind_w),
            "Need all wind components, has [u: {}, v: {}, w: {}]",
            atmospheric_field.has(wind_u),
            atmospheric_field.has(wind_v),
            atmospheric_field.has(wind_w));
      }

      if (atmospheric_field.has(mag_u) or atmospheric_field.has(mag_v) or
          atmospheric_field.has(mag_w)) {
        ARTS_USER_ERROR_IF(not atmospheric_field.has(mag_u, mag_v, mag_w),
                           "Need all mag components, has [u: {}, v: {}, w: {}]",
                           atmospheric_field.has(mag_u),
                           atmospheric_field.has(mag_v),
                           atmospheric_field.has(mag_w));
      }
      break;
    case MissingFieldComponentError::Zero:
      if (atmospheric_field.has(wind_u) or atmospheric_field.has(wind_v) or
          atmospheric_field.has(wind_w)) {
        if (not atmospheric_field.has(wind_u)) atmospheric_field[wind_u] = 0.0;
        if (not atmospheric_field.has(wind_v)) atmospheric_field[wind_v] = 0.0;
        if (not atmospheric_field.has(wind_w)) atmospheric_field[wind_w] = 0.0;
      }

      if (atmospheric_field.has(mag_u) or atmospheric_field.has(mag_v) or
          atmospheric_field.has(mag_w)) {
        if (not atmospheric_field.has(mag_u)) atmospheric_field[mag_u] = 0.0;
        if (not atmospheric_field.has(mag_v)) atmospheric_field[mag_v] = 0.0;
        if (not atmospheric_field.has(mag_w)) atmospheric_field[mag_w] = 0.0;
      }
      break;

    case MissingFieldComponentError::Ignore:
      break;
  }
}

void keysSpecies(std::unordered_map<SpeciesEnum, Index> keys,
                 const ArrayOfAbsorptionBand &absorption_bands) {
  if (absorption_bands.empty()) return;

  for (auto &[key, value] : absorption_bands) {
    ++keys[key.Species()];

    for (auto &line : value.lines) {
      for (auto &ls : line.ls.single_models) {
        keys[ls.species];
      }
    }
  }
}

void atmospheric_fieldAppendLineSpeciesData(
    AtmField &atmospheric_field,
    const ArrayOfAbsorptionBand &absorption_bands,
    const String &basename,
    const String &extrapolation,
    const Index &missing_is_zero,
    const Index &replace_existing) {
  std::unordered_map<SpeciesEnum, Index> keys;
  keysSpecies(keys, absorption_bands);

  append_data(atmospheric_field,
              basename,
              extrapolation,
              missing_is_zero,
              replace_existing,
              0,
              keys,
              [](const SpeciesEnum &x) { return String{toString<1>(x)}; });
}

void keysIsotopologue(std::unordered_map<SpeciesIsotope, Index> keys,
                      const ArrayOfAbsorptionBand &absorption_bands) {
  if (absorption_bands.empty()) return;

  for (auto &[key, value] : absorption_bands) {
    ++keys[key.Isotopologue()];
  }
}

void atmospheric_fieldAppendLineIsotopologueData(
    AtmField &atmospheric_field,
    const ArrayOfAbsorptionBand &absorption_bands,
    const String &basename,
    const String &extrapolation,
    const Index &missing_is_zero,
    const Index &replace_existing) {
  std::unordered_map<SpeciesIsotope, Index> keys;
  keysIsotopologue(keys, absorption_bands);

  for (auto &[key, value] : absorption_bands) {
    ++keys[key.Isotopologue()];
  }

  append_data(atmospheric_field,
              basename,
              extrapolation,
              missing_is_zero,
              replace_existing,
              0,
              keys,
              [](const SpeciesIsotope &x) { return x.FullName(); });
}

void keysNLTE(std::unordered_map<QuantumIdentifier, Index> keys,
              const ArrayOfAbsorptionBand &absorption_bands) {
  if (absorption_bands.empty()) return;

  for (auto &[key, value] : absorption_bands) {
    ++keys[key.UpperLevel()];
    ++keys[key.LowerLevel()];
  }
}

void atmospheric_fieldAppendLineLevelData(
    AtmField &atmospheric_field,
    const ArrayOfAbsorptionBand &absorption_bands,
    const String &basename,
    const String &extrapolation,
    const Index &missing_is_zero,
    const Index &replace_existing) {
  std::unordered_map<QuantumIdentifier, Index> keys;
  keysNLTE(keys, absorption_bands);

  for (auto &[key, value] : absorption_bands) {
    ++keys[key.UpperLevel()];
    ++keys[key.LowerLevel()];
  }

  append_data(atmospheric_field,
              basename,
              extrapolation,
              missing_is_zero,
              replace_existing,
              0,
              keys,
              [](const QuantumIdentifier &x) { return var_string(x); });
}

void keysSpecies(std::unordered_map<SpeciesEnum, Index> &keys,
                 const ArrayOfArrayOfSpeciesTag &absorption_species) {
  if (absorption_species.empty()) return;

  for (auto &species_tags : absorption_species) {
    ++keys[species_tags.Species()];
  }
}

void atmospheric_fieldAppendTagsSpeciesData(
    AtmField &atmospheric_field,
    const ArrayOfArrayOfSpeciesTag &absorption_species,
    const String &basename,
    const String &extrapolation,
    const Index &missing_is_zero,
    const Index &replace_existing) {
  std::unordered_map<SpeciesEnum, Index> keys;
  keysSpecies(keys, absorption_species);

  append_data(atmospheric_field,
              basename,
              extrapolation,
              missing_is_zero,
              replace_existing,
              0,
              keys,
              [](const SpeciesEnum &x) { return String{toString<1>(x)}; });
}

void atmospheric_fieldAppendLookupTableSpeciesData(
    AtmField &atmospheric_field,
    const GasAbsLookup &absorption_lookup_table_data,
    const String &basename,
    const String &extrapolation,
    const Index &missing_is_zero,
    const Index &replace_existing) {
  atmospheric_fieldAppendTagsSpeciesData(atmospheric_field,
                                         absorption_lookup_table_data.Species(),
                                         basename,
                                         extrapolation,
                                         missing_is_zero,
                                         replace_existing);
}

void keysSpecies(std::unordered_map<SpeciesEnum, Index> &keys,
                 const ArrayOfCIARecord &absorption_cia_data) {
  if (absorption_cia_data.empty()) return;

  for (auto &cia_record : absorption_cia_data) {
    const auto [spec1, spec2] = cia_record.TwoSpecies();
    ++keys[spec1];
    ++keys[spec2];
  }
}

void atmospheric_fieldAppendCIASpeciesData(
    AtmField &atmospheric_field,
    const ArrayOfCIARecord &absorption_cia_data,
    const String &basename,
    const String &extrapolation,
    const Index &missing_is_zero,
    const Index &replace_existing) {
  std::unordered_map<SpeciesEnum, Index> keys;
  keysSpecies(keys, absorption_cia_data);

  append_data(atmospheric_field,
              basename,
              extrapolation,
              missing_is_zero,
              replace_existing,
              0,
              keys,
              [](const SpeciesEnum &x) { return String{toString<1>(x)}; });
}

void keysSpecies(std::unordered_map<SpeciesEnum, Index> &keys,
                 const ArrayOfXsecRecord &absorption_xsec_fit_data) {
  if (absorption_xsec_fit_data.empty()) return;

  for (auto &xsec_record : absorption_xsec_fit_data) {
    ++keys[xsec_record.Species()];
  }
}

void atmospheric_fieldAppendXsecSpeciesData(
    AtmField &atmospheric_field,
    const ArrayOfXsecRecord &absorption_xsec_fit_data,
    const String &basename,
    const String &extrapolation,
    const Index &missing_is_zero,
    const Index &replace_existing) {
  std::unordered_map<SpeciesEnum, Index> keys;
  keysSpecies(keys, absorption_xsec_fit_data);

  append_data(atmospheric_field,
              basename,
              extrapolation,
              missing_is_zero,
              replace_existing,
              0,
              keys,
              [](const SpeciesEnum &x) { return String{toString<1>(x)}; });
}

void keysSpecies(std::unordered_map<SpeciesEnum, Index> &keys,
                 const PredefinedModelData &absorption_predefined_model_data) {
  if (absorption_predefined_model_data.data.empty()) return;

  for (auto &predef_record : absorption_predefined_model_data.data) {
    ++keys[predef_record.first.spec];
  }

  for (auto &spec : Absorption::PredefinedModel::VMRS::species) {
    ++keys[spec];
  }
}

void atmospheric_fieldAppendPredefSpeciesData(
    AtmField &atmospheric_field,
    const PredefinedModelData &absorption_predefined_model_data,
    const String &basename,
    const String &extrapolation,
    const Index &missing_is_zero,
    const Index &replace_existing) {
  std::unordered_map<SpeciesEnum, Index> keys;
  keysSpecies(keys, absorption_predefined_model_data);

  append_data(atmospheric_field,
              basename,
              extrapolation,
              missing_is_zero,
              replace_existing,
              0,
              keys,
              [](const SpeciesEnum &x) { return String{toString<1>(x)}; });
}

void atmospheric_fieldAppendAbsorptionData(const Workspace &ws,
                                           AtmField &atmospheric_field,
                                           const String &basename,
                                           const String &extrapolation,
                                           const Index &missing_is_zero,
                                           const Index &replace_existing,
                                           const Index &load_isot,
                                           const Index &load_nlte) {
  std::unordered_map<SpeciesEnum, Index> keys;

  if (const String lines_str = "absorption_bands";
      ws.wsv_and_contains(lines_str)) {
    using lines_t    = ArrayOfAbsorptionBand;
    const auto &data = ws.get<lines_t>(lines_str);

    keysSpecies(keys, data);
    if (static_cast<bool>(load_isot)) {
      atmospheric_fieldAppendLineIsotopologueData(atmospheric_field,
                                                  data,
                                                  basename,
                                                  extrapolation,
                                                  missing_is_zero,
                                                  replace_existing);
    }

    if (static_cast<bool>(load_nlte)) {
      atmospheric_fieldAppendLineLevelData(atmospheric_field,
                                           data,
                                           basename,
                                           extrapolation,
                                           missing_is_zero,
                                           replace_existing);
    }
  }

  if (const String cia_str = "absorption_cia_data";
      ws.wsv_and_contains(cia_str)) {
    using cia_t = ArrayOfCIARecord;
    keysSpecies(keys, ws.get<cia_t>(cia_str));
  }

  if (const String xsec_str = "absorption_xsec_fit_data";
      ws.wsv_and_contains(xsec_str)) {
    using xsec_t = ArrayOfXsecRecord;
    keysSpecies(keys, ws.get<xsec_t>(xsec_str));
  }

  if (const String predef_str = "absorption_predefined_model_data";
      ws.wsv_and_contains(predef_str)) {
    using predef_t = PredefinedModelData;
    keysSpecies(keys, ws.get<predef_t>(predef_str));
  }

  if (const String species_str = "absorption_species";
      ws.wsv_and_contains(species_str)) {
    using aospec_t = ArrayOfArrayOfSpeciesTag;
    keysSpecies(keys, ws.get<aospec_t>(species_str));
  }

  append_data(atmospheric_field,
              basename,
              extrapolation,
              missing_is_zero,
              replace_existing,
              0,
              keys,
              [](const SpeciesEnum &x) { return String{toString<1>(x)}; });
}

void atmospheric_fieldIGRF(AtmField &atmospheric_field, const Time &time) {
  using IGRF::igrf;

  //! We need explicit planet-size as IGRF requires the radius
  //! This is the WGS84 version of that, with radius of equator and pole
  static constexpr Vector2 ell{6378137., 6356752.314245};

  //! This struct deals with the computations, it's internally saving a mutable
  //! state
  struct res {
    Time t;

    [[nodiscard]] Vector3 comp(Numeric al, Numeric la, Numeric lo) const {
      return igrf({al, la, lo}, ell, t);
    }

    [[nodiscard]] Numeric get_u(Numeric z, Numeric la, Numeric lo) const {
      return igrf({z, la, lo}, ell, t)[0];
    }

    [[nodiscard]] Numeric get_v(Numeric z, Numeric la, Numeric lo) const {
      return igrf({z, la, lo}, ell, t)[1];
    }

    [[nodiscard]] Numeric get_w(Numeric z, Numeric la, Numeric lo) const {
      return igrf({z, la, lo}, ell, t)[2];
    }
  };

  atmospheric_field[AtmKey::mag_u] = Atm::FunctionalData{
      [cpy = res(time)](Numeric h, Numeric lat, Numeric lon) {
        return cpy.get_u(h, lat, lon);
      }};
  atmospheric_field[AtmKey::mag_v] = Atm::FunctionalData{
      [cpy = res(time)](Numeric h, Numeric lat, Numeric lon) {
        return cpy.get_v(h, lat, lon);
      }};
  atmospheric_field[AtmKey::mag_w] = Atm::FunctionalData{
      [cpy = res(time)](Numeric h, Numeric lat, Numeric lon) {
        return cpy.get_w(h, lat, lon);
      }};
}

enum class atmospheric_fieldHydrostaticPressureDataOptions : char {
  Lat,
  Lon,
  Hypsometric,
  Hydrostatic,
};

template <atmospheric_fieldHydrostaticPressureDataOptions... input_opts>
struct atmospheric_fieldHydrostaticPressureData {
  using enum atmospheric_fieldHydrostaticPressureDataOptions;
  static constexpr std::array<atmospheric_fieldHydrostaticPressureDataOptions,
                              sizeof...(input_opts)>
      opts{input_opts...};
  static constexpr bool do_lat = std::ranges::any_of(opts, Cmp::eq(Lat));
  static constexpr bool do_lon = std::ranges::any_of(opts, Cmp::eq(Lon));

  using LatLag = my_interp::Lagrange<Index{do_lat}>;
  using LonLag = my_interp::Lagrange<Index{do_lon},
                                     false,
                                     GridType::Cyclic,
                                     my_interp::cycle_m180_p180>;

  template <atmospheric_fieldHydrostaticPressureDataOptions X>
  static constexpr bool is = std::ranges::any_of(opts, Cmp::eq(X));

  Tensor3 grad_p;
  Tensor3 pre;
  Vector alt;
  Vector lat;
  Vector lon;

  static Numeric step(Numeric p, Numeric h, Numeric d) {
    if constexpr (is<Hypsometric>) {
      return p * std::exp(-h * d);
    } else if constexpr (is<Hydrostatic>) {
      return std::max<Numeric>(0, p * (1.0 - h * d));
    }
  }

  atmospheric_fieldHydrostaticPressureData(Tensor3 in_grad_p,
                                           const GriddedField2 &pre0,
                                           Vector in_alt)
      : grad_p(std::move(in_grad_p)),
        pre(grad_p),  // Init sizes
        alt(std::move(in_alt)),
        lat(pre0.grid<0>()),
        lon(pre0.grid<1>()) {
    pre[0] = pre0.data;
    for (Index i = 1; i < alt.nelem(); i++) {
      for (Index j = 0; j < lat.nelem(); j++) {
        for (Index k = 0; k < lon.nelem(); k++) {
          const Numeric h  = alt[i] - alt[i - 1];
          const Numeric p0 = pre(i - 1, j, k);
          const Numeric d0 = grad_p(i - 1, j, k);

          pre(i, j, k) = step(p0, h, d0);
        }
      }
    }
  }

  [[nodiscard]] std::pair<Index, Numeric> find_alt(Numeric al) const {
    auto i  = std::distance(alt.begin(), std::ranges::upper_bound(alt, al));
    i      -= (i == alt.size());
    while (i > 0 and alt[i] > al) {
      i--;
    }

    return {i, al - alt[i]};
  }

  [[nodiscard]] std::pair<Numeric, Numeric> level(Index alt_ind,
                                                  Numeric la,
                                                  Numeric lo) const {
    const LatLag latlag(0, la, lat);
    const LonLag lonlag(0, lo, lon);
    const auto iw   = interpweights(latlag, lonlag);
    const Numeric p = interp(pre[alt_ind], iw, latlag, lonlag);
    const Numeric d = interp(grad_p[alt_ind], iw, latlag, lonlag);
    return {p, d};
  }

  Numeric operator()(Numeric al, Numeric la, Numeric lo) const {
    const auto [i, h] = find_alt(al);
    const auto [p, d] = level(i, la, lo);
    return step(p, h, d);
  }
};

void atmospheric_fieldHydrostaticPressure(
    AtmField &atmospheric_field,
    const NumericTernaryOperator &gravity_operator,
    const GriddedField2 &p0,
    const Vector &alts,
    const Numeric &fixed_specific_gas_constant,
    const Numeric &fixed_atm_temperature,
    const String &hydrostatic_option) {
  using enum atmospheric_fieldHydrostaticPressureDataOptions;
  using enum HydrostaticPressureOption;

  const Vector &lats = p0.grid<0>();
  const Vector &lons = p0.grid<1>();

  const Index nalt = alts.size();
  const Index nlat = lats.size();
  const Index nlon = lons.size();

  ARTS_USER_ERROR_IF(nalt < 1 or nlat < 1 or nlon < 1,
                     "Must have at least 1-sized alt, lat, and lon grids")

  ARTS_USER_ERROR_IF(
      p0.grid_names[0] not_eq "Latitude" or
          p0.grid_names[1] not_eq "Longitude" or not p0.ok(),
      "Bad gridded field, must have right size.\n"
      "Must also have \"Latitude\" as first grid and \"Longitude\" as second grid.\n"
      "Field:{}\n",
      p0)

  const bool has_def_t = fixed_atm_temperature > 0;
  const bool has_def_r = fixed_specific_gas_constant > 0;

  ARTS_USER_ERROR_IF(
      not has_def_t and not atmospheric_field.contains(AtmKey::t),
      "atmospheric_field lacks temperature and no default temperature given")

  ARTS_USER_ERROR_IF(
      not has_def_r and atmospheric_field.nspec() == 0,
      "atmospheric_field lacks species and no default specific gas constant given")

  const Tensor3 scale_factor = [&]() {
    Tensor3 scl(nalt, nlat, nlon);
    for (Index i = 0; i < nalt; i++) {
      for (Index j = 0; j < nlat; j++) {
        for (Index k = 0; k < nlon; k++) {
          const Numeric al = alts[i];
          const Numeric la = lats[j];
          const Numeric lo = lons[k];

          const Numeric g = gravity_operator(al, la, lo);
          const AtmPoint atmospheric_point{atmospheric_field.at(al, la, lo)};

          const Numeric inv_specific_gas_constant =
              has_def_r ? 1.0 / fixed_specific_gas_constant
                        : (1e-3 * atmospheric_point.mean_mass() / Constant::R);
          const Numeric inv_temp = has_def_t
                                       ? 1.0 / fixed_atm_temperature
                                       : 1.0 / atmospheric_point.temperature;

          // Partial rho, no pressure
          scl(i, j, k) = g * inv_specific_gas_constant * inv_temp;
        }
      }
    }

    return scl;
  }();

  switch (to<HydrostaticPressureOption>(hydrostatic_option)) {
    case HypsometricEquation:
      if (nlon > 1 and nlat > 1) {
        atmospheric_field[AtmKey::p] = Atm::FunctionalData{
            atmospheric_fieldHydrostaticPressureData<Hypsometric, Lat, Lon>(
                scale_factor, p0, alts)};
      } else if (nlat > 1) {
        atmospheric_field[AtmKey::p] = Atm::FunctionalData{
            atmospheric_fieldHydrostaticPressureData<Hypsometric, Lat>(
                scale_factor, p0, alts)};
      } else if (nlon > 1) {
        atmospheric_field[AtmKey::p] = Atm::FunctionalData{
            atmospheric_fieldHydrostaticPressureData<Hypsometric, Lon>(
                scale_factor, p0, alts)};
      } else {
        atmospheric_field[AtmKey::p] = Atm::FunctionalData{
            atmospheric_fieldHydrostaticPressureData<Hypsometric>(
                scale_factor, p0, alts)};
      }
      break;

    case HydrostaticEquation:
      if (nlon > 1 and nlat > 1) {
        atmospheric_field[AtmKey::p] = Atm::FunctionalData{
            atmospheric_fieldHydrostaticPressureData<Hydrostatic, Lat, Lon>(
                scale_factor, p0, alts)};
      } else if (nlat > 1) {
        atmospheric_field[AtmKey::p] = Atm::FunctionalData{
            atmospheric_fieldHydrostaticPressureData<Hydrostatic, Lat>(
                scale_factor, p0, alts)};
      } else if (nlon > 1) {
        atmospheric_field[AtmKey::p] = Atm::FunctionalData{
            atmospheric_fieldHydrostaticPressureData<Hydrostatic, Lon>(
                scale_factor, p0, alts)};
      } else {
        atmospheric_field[AtmKey::p] = Atm::FunctionalData{
            atmospheric_fieldHydrostaticPressureData<Hydrostatic>(
                scale_factor, p0, alts)};
      }
      break;
  }
}
