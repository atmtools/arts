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
#include "atm_mag_field.h"
#include "compare.h"
#include "configtypes.h"
#include "debug.h"
#include "enumsFieldComponent.h"
#include "igrf13.h"
#include "interp.h"
#include "interpolation.h"
#include "isotopologues.h"
#include "lbl_data.h"
#include "operator_io.h"
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
  ARTS_TIME_REPORT

  atmospheric_field = AtmField{to<IsoRatioOption>(default_isotopologue)};
  atmospheric_field.top_of_atmosphere = top_of_atmosphere;
}

void atmospheric_pointInit(AtmPoint &atmospheric_point,
                           const String &default_isotopologue) {
  ARTS_TIME_REPORT

  atmospheric_point = AtmPoint{to<IsoRatioOption>(default_isotopologue)};
}

template <Size I = 0, typename... T, Size N = sizeof...(T)>
std::variant<T...> xml_read_from_file_variant(const std::variant<T...> &_,
                                              const String &filename)
  requires(N != 0)
{
  using myT = typename std::tuple_element<I, std::tuple<T...>>::type;

  std::string error{};
  try {
    try {
      myT x{};
      xml_read_from_file(filename, x);
      return x;
    } catch (std::exception &e) {
      error = std::format(
          "Group {} reports:\n{}", xml_io_stream<myT>::type_name, e.what());
    }

    if constexpr (I + 1 < N) {
      return xml_read_from_file_variant<I + 1, T...>(_, filename);
    }

    throw std::runtime_error(std::format(
        "Could not read valid data from file: {}\nTried all groups in {}",
        filename,
        xml_io_stream<std::variant<T...>>::type_name));
  } catch (std::exception &e) {
    ARTS_USER_ERROR("{}\n\n{}", e.what(), error);
  }
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
    const auto &to_string = [](const auto &x) { return std::format("{}", x); })
  requires(
      std::same_as<typename std::decay_t<decltype(to_string(T{}))>, String>)
{
  const String my_base = complete_basename(basename);

  Atm::Data x;
  x.alt_low = x.alt_upp = x.lat_low = x.lat_upp = x.lon_low = x.lon_upp =
      to<InterpolationExtrapolation>(extrapolation);

  for (const auto &[key, _] : keys) {
    String filename = std::format("{}{}.xml", my_base, to_string(key));

    if (static_cast<bool>(replace_existing) or
        not atmospheric_field.contains(key)) {
      if (find_xml_file_existence(filename)) {
        x.data = xml_read_from_file_variant(Atm::FieldData{1.0}, filename);
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
  ARTS_TIME_REPORT

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
              [](const AtmKey &x) { return std::format("{}", x); });

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

    case MissingFieldComponentError::Ignore: break;
  }
}

void keysSpecies(std::unordered_map<SpeciesEnum, Index> &keys,
                 const AbsorptionBands &absorption_bands) {
  if (absorption_bands.empty()) return;

  for (auto &[key, value] : absorption_bands) {
    ++keys[key.Species()];

    for (auto &line : value.lines) {
      for (auto &ls : line.ls.single_models) {
        ++keys[ls.species];
      }
    }
  }
}

void atmospheric_fieldAppendLineSpeciesData(
    AtmField &atmospheric_field,
    const AbsorptionBands &absorption_bands,
    const String &basename,
    const String &extrapolation,
    const Index &missing_is_zero,
    const Index &replace_existing) {
  ARTS_TIME_REPORT

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
                      const AbsorptionBands &absorption_bands) {
  if (absorption_bands.empty()) return;

  for (auto &[key, value] : absorption_bands) {
    ++keys[key.Isotopologue()];
  }
}

void atmospheric_fieldAppendLineIsotopologueData(
    AtmField &atmospheric_field,
    const AbsorptionBands &absorption_bands,
    const String &basename,
    const String &extrapolation,
    const Index &missing_is_zero,
    const Index &replace_existing) {
  ARTS_TIME_REPORT

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
              const AbsorptionBands &absorption_bands) {
  if (absorption_bands.empty()) return;

  for (auto &[key, value] : absorption_bands) {
    ++keys[key.UpperLevel()];
    ++keys[key.LowerLevel()];
  }
}

void atmospheric_fieldAppendLineLevelData(
    AtmField &atmospheric_field,
    const AbsorptionBands &absorption_bands,
    const String &basename,
    const String &extrapolation,
    const Index &missing_is_zero,
    const Index &replace_existing) {
  ARTS_TIME_REPORT

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
              [](const QuantumIdentifier &x) { return std::format("{}", x); });
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
  ARTS_TIME_REPORT

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
  ARTS_TIME_REPORT

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
                 const AbsorptionLookupTables &absorption_lookup_table) {
  if (absorption_lookup_table.empty()) return;

  for (auto &&spec : absorption_lookup_table | stdv::keys) {
    ++keys[spec];
  }
}

void atmospheric_fieldAppendLookupTableSpeciesData(
    AtmField &atmospheric_field,
    const AbsorptionLookupTables &absorption_lookup_table,
    const String &basename,
    const String &extrapolation,
    const Index &missing_is_zero,
    const Index &replace_existing) {
  ARTS_TIME_REPORT

  std::unordered_map<SpeciesEnum, Index> keys;
  keysSpecies(keys, absorption_lookup_table);

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
  ARTS_TIME_REPORT

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
  if (absorption_predefined_model_data.empty()) return;

  for (auto &predef_record : absorption_predefined_model_data) {
    ++keys[predef_record.first.spec];
  }
}

void atmospheric_fieldAppendPredefSpeciesData(
    AtmField &atmospheric_field,
    const PredefinedModelData &absorption_predefined_model_data,
    const String &basename,
    const String &extrapolation,
    const Index &missing_is_zero,
    const Index &replace_existing) {
  ARTS_TIME_REPORT

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
  ARTS_TIME_REPORT

  if (const String lines_str = "absorption_bands";
      ws.wsv_and_contains(lines_str)) {
    using lines_t                = AbsorptionBands;
    const auto &absorption_bands = ws.get<lines_t>(lines_str);

    atmospheric_fieldAppendLineSpeciesData(atmospheric_field,
                                           absorption_bands,
                                           basename,
                                           extrapolation,
                                           missing_is_zero,
                                           replace_existing);

    if (static_cast<bool>(load_isot)) {
      atmospheric_fieldAppendLineIsotopologueData(atmospheric_field,
                                                  absorption_bands,
                                                  basename,
                                                  extrapolation,
                                                  missing_is_zero,
                                                  replace_existing);
    }

    if (static_cast<bool>(load_nlte)) {
      atmospheric_fieldAppendLineLevelData(atmospheric_field,
                                           absorption_bands,
                                           basename,
                                           extrapolation,
                                           missing_is_zero,
                                           replace_existing);
    }
  }

  if (const String cia_str = "absorption_cia_data";
      ws.wsv_and_contains(cia_str)) {
    using cia_t                     = ArrayOfCIARecord;
    const auto &absorption_cia_data = ws.get<cia_t>(cia_str);
    atmospheric_fieldAppendCIASpeciesData(atmospheric_field,
                                          absorption_cia_data,
                                          basename,
                                          extrapolation,
                                          missing_is_zero,
                                          replace_existing);
  }

  if (const String lookup_str = "absorption_lookup_table";
      ws.wsv_and_contains(lookup_str)) {
    using lookup_t                      = AbsorptionLookupTables;
    const auto &absorption_lookup_table = ws.get<lookup_t>(lookup_str);
    atmospheric_fieldAppendLookupTableSpeciesData(atmospheric_field,
                                                  absorption_lookup_table,
                                                  basename,
                                                  extrapolation,
                                                  missing_is_zero,
                                                  replace_existing);
  }

  if (const String xsec_str = "absorption_xsec_fit_data";
      ws.wsv_and_contains(xsec_str)) {
    using xsec_t                         = ArrayOfXsecRecord;
    const auto &absorption_xsec_fit_data = ws.get<xsec_t>(xsec_str);
    atmospheric_fieldAppendXsecSpeciesData(atmospheric_field,
                                           absorption_xsec_fit_data,
                                           basename,
                                           extrapolation,
                                           missing_is_zero,
                                           replace_existing);
  }

  if (const String predef_str = "absorption_predefined_model_data";
      ws.wsv_and_contains(predef_str)) {
    using predef_t                               = PredefinedModelData;
    const auto &absorption_predefined_model_data = ws.get<predef_t>(predef_str);
    atmospheric_fieldAppendPredefSpeciesData(atmospheric_field,
                                             absorption_predefined_model_data,
                                             basename,
                                             extrapolation,
                                             missing_is_zero,
                                             replace_existing);
  }

  if (const String species_str = "absorption_species";
      ws.wsv_and_contains(species_str)) {
    using aospec_t                 = ArrayOfArrayOfSpeciesTag;
    const auto &absorption_species = ws.get<aospec_t>(species_str);
    atmospheric_fieldAppendTagsSpeciesData(atmospheric_field,
                                           absorption_species,
                                           basename,
                                           extrapolation,
                                           missing_is_zero,
                                           replace_existing);
  }
}

void atmospheric_fieldIGRF(AtmField &atmospheric_field, const Time &time) {
  ARTS_TIME_REPORT

  atmospheric_field[AtmKey::mag_u] =
      NumericTernaryOperator{Atm::IGRF13(time, FieldComponent::u)};
  atmospheric_field[AtmKey::mag_v] =
      NumericTernaryOperator{Atm::IGRF13(time, FieldComponent::v)};
  atmospheric_field[AtmKey::mag_w] =
      NumericTernaryOperator{Atm::IGRF13(time, FieldComponent::w)};
}

void atmospheric_fieldHydrostaticPressure(
    AtmField &atmospheric_field,
    const NumericTernaryOperator &gravity_operator,
    const AscendingGrid &alts,
    const SortedGriddedField2 &p0,
    const Numeric &fixed_specific_gas_constant,
    const Numeric &fixed_atm_temperature,
    const String &hydrostatic_option) {
  ARTS_TIME_REPORT

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
          scl[i, j, k] = g * inv_specific_gas_constant * inv_temp;
        }
      }
    }

    return scl;
  }();

  atmospheric_field[AtmKey::p] = Atm::FunctionalData{Atm::HydrostaticPressure(
      scale_factor,
      p0,
      alts,
      to<HydrostaticPressureOption>(hydrostatic_option))};
}

void atmospheric_fieldHydrostaticPressure(
    AtmField &atmospheric_field,
    const NumericTernaryOperator &gravity_operator,
    const AscendingGrid &alts,
    const Numeric &p0,
    const Numeric &fixed_specific_gas_constant,
    const Numeric &fixed_atm_temperature,
    const String &hydrostatic_option) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      not atmospheric_field.has(AtmKey::t),
      "Must have a temperature field to call this workspace method with a single Numeric reference pressure, so that latitude and longitude grids can be extracted")

  const auto &t = atmospheric_field[AtmKey::t];

  ARTS_USER_ERROR_IF(
      not std::holds_alternative<SortedGriddedField3>(t.data),
      "Temperature field must be a SortedGriddedField3 to call this workspace method with a single Numeric reference pressure, so that latitude and longitude grids can be extracted")

  const auto &t0 = std::get<SortedGriddedField3>(t.data);
  const SortedGriddedField2 p0_field{
      .data_name  = "Pressure",
      .data       = Matrix(t0.grid<1>().size(), t0.grid<2>().size(), p0),
      .grid_names = {t0.gridname<1>(), t0.gridname<2>()},
      .grids      = {t0.grid<1>(), t0.grid<2>()}};

  atmospheric_fieldHydrostaticPressure(atmospheric_field,
                                       gravity_operator,
                                       alts,
                                       p0_field,
                                       fixed_specific_gas_constant,
                                       fixed_atm_temperature,
                                       hydrostatic_option);
}

template <Atm::KeyType T>
void atmospheric_fieldRegridTemplate(AtmData &data,
                                     const T &key,
                                     const AscendingGrid &alt,
                                     const AscendingGrid &lat,
                                     const AscendingGrid &lon,
                                     const String &extrapolation) {
  ARTS_USER_ERROR_IF(alt.size() * lat.size() * lon.size() == 0,
                     R"(Cannot regrid to empty grid
alt: {:Bs,} [{} elements]
lat: {:Bs,} [{} elements]
lon: {:Bs,} [{} elements]
)",
                     alt,
                     alt.size(),
                     lat,
                     lat.size(),
                     lon,
                     lon.size())

  const InterpolationExtrapolation extrap =
      to<InterpolationExtrapolation>(extrapolation);

  SortedGriddedField3 new_field{
      .data_name  = String{std::format("{}", key)},
      .data       = Tensor3(alt.size(), lat.size(), lon.size()),
      .grid_names = {String{"Altitude"},
                     String{"Latitude"},
                     String{"Longitude"}},
      .grids      = {alt, lat, lon},
  };

  for (Size i = 0; i < alt.size(); i++) {
    for (Size j = 0; j < lat.size(); j++) {
      for (Size k = 0; k < lon.size(); k++) {
        new_field[i, j, k] = data.at(alt[i], lat[j], lon[k]);
      }
    }
  }

  data.data    = std::move(new_field);
  data.alt_upp = extrap;
  data.alt_low = extrap;
  data.lat_upp = extrap;
  data.lat_low = extrap;
  data.lon_upp = extrap;
  data.lon_low = extrap;
}

void atmospheric_fieldRegrid(AtmField &atmospheric_field,
                             const AscendingGrid &alt,
                             const AscendingGrid &lat,
                             const AscendingGrid &lon,
                             const ScatteringSpeciesProperty &key,
                             const String &extrapolation) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(not atmospheric_field.contains(key),
                     R"(No scattering species key "{}" in atmospheric_field)",
                     key)

  atmospheric_fieldRegridTemplate(
      atmospheric_field[key], key, alt, lat, lon, extrapolation);
}

void atmospheric_fieldRegrid(AtmField &atmospheric_field,
                             const AscendingGrid &alt,
                             const AscendingGrid &lat,
                             const AscendingGrid &lon,
                             const SpeciesEnum &key,
                             const String &extrapolation) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(not atmospheric_field.contains(key),
                     R"(No VMR key "{}" in atmospheric_field)",
                     key)

  atmospheric_fieldRegridTemplate(
      atmospheric_field[key], key, alt, lat, lon, extrapolation);
}

void atmospheric_fieldRegrid(AtmField &atmospheric_field,
                             const AscendingGrid &alt,
                             const AscendingGrid &lat,
                             const AscendingGrid &lon,
                             const SpeciesIsotope &key,
                             const String &extrapolation) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(not atmospheric_field.contains(key),
                     R"(No isotopologue ratio key "{}" in atmospheric_field)",
                     key)

  atmospheric_fieldRegridTemplate(
      atmospheric_field[key], key, alt, lat, lon, extrapolation);
}

void atmospheric_fieldRegrid(AtmField &atmospheric_field,
                             const AscendingGrid &alt,
                             const AscendingGrid &lat,
                             const AscendingGrid &lon,
                             const QuantumIdentifier &key,
                             const String &extrapolation) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(not atmospheric_field.contains(key),
                     R"(No NLTE key "{}" in atmospheric_field)",
                     key)

  atmospheric_fieldRegridTemplate(
      atmospheric_field[key], key, alt, lat, lon, extrapolation);
}

void atmospheric_fieldRegrid(AtmField &atmospheric_field,
                             const AscendingGrid &alt,
                             const AscendingGrid &lat,
                             const AscendingGrid &lon,
                             const AtmKey &key,
                             const String &extrapolation) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(not atmospheric_field.contains(key),
                     R"(No atmospheric key "{}" in atmospheric_field)",
                     key)

  atmospheric_fieldRegridTemplate(
      atmospheric_field[key], key, alt, lat, lon, extrapolation);
}

void atmospheric_fieldRegridAll(AtmField &atmospheric_field,
                                const AscendingGrid &alt,
                                const AscendingGrid &lat,
                                const AscendingGrid &lon,
                                const String &extrapolation) {
  ARTS_TIME_REPORT

  for (auto &[key, data] : atmospheric_field.map<ScatteringSpeciesProperty>()) {
    atmospheric_fieldRegridTemplate(data, key, alt, lat, lon, extrapolation);
  }

  for (auto &[key, data] : atmospheric_field.map<SpeciesEnum>()) {
    atmospheric_fieldRegridTemplate(data, key, alt, lat, lon, extrapolation);
  }

  for (auto &[key, data] : atmospheric_field.map<SpeciesIsotope>()) {
    atmospheric_fieldRegridTemplate(data, key, alt, lat, lon, extrapolation);
  }

  for (auto &[key, data] : atmospheric_field.map<QuantumIdentifier>()) {
    atmospheric_fieldRegridTemplate(data, key, alt, lat, lon, extrapolation);
  }

  for (auto &[key, data] : atmospheric_field.map<AtmKey>()) {
    atmospheric_fieldRegridTemplate(data, key, alt, lat, lon, extrapolation);
  }
}

void atmospheric_fieldFromProfile(
    AtmField &atmospheric_field,
    const ArrayOfAtmPoint &atmospheric_profile,
    const AscendingGrid &altitude_grid,
    const InterpolationExtrapolation &altitude_extrapolation) {
  ARTS_USER_ERROR_IF(
      not arr::same_size(atmospheric_profile, altitude_grid),
      "Mismatch in size between atmospheric_profile and altitude_grid");
  ARTS_USER_ERROR_IF(altitude_grid.empty(), "Empty altitude grid")

  atmospheric_field = Atm::atm_from_profile(atmospheric_profile,
                                            altitude_grid,
                                            altitude_extrapolation,
                                            altitude_grid.back());
}
