#include <array_algo.h>
#include <debug.h>
#include <enumsHydrostaticPressureOption.h>
#include <enumsIsoRatioOption.h>
#include <enumsMissingFieldComponentError.h>
#include <geodetic.h>
#include <workspace.h>
#include <zconf.h>

#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <variant>

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

namespace {
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

void keysSpecies(std::unordered_map<SpeciesEnum, Index> &keys,
                 const AbsorptionBands &absorption_bands) {
  if (absorption_bands.empty()) return;

  for (auto &[key, value] : absorption_bands) {
    ++keys[key.isot.spec];

    for (auto &line : value.lines) {
      for (auto &ls : line.ls.single_models) {
        if (ls.first != "AIR"_spec) ++keys[ls.first];
      }
    }
  }
}

void keysIsotopologue(std::unordered_map<SpeciesIsotope, Index> &keys,
                      const AbsorptionBands &absorption_bands) {
  if (absorption_bands.empty()) return;

  for (auto &[key, value] : absorption_bands) {
    ++keys[key.isot];
  }
}

void keysNLTE(std::unordered_map<QuantumLevelIdentifier, Index> &keys,
              const AbsorptionBands &absorption_bands) {
  if (absorption_bands.empty()) return;

  for (auto &[key, value] : absorption_bands) {
    ++keys[key.upper()];
    ++keys[key.lower()];
  }
}

void keysSpecies(std::unordered_map<SpeciesEnum, Index> &keys,
                 const ArrayOfArrayOfSpeciesTag &absorption_species) {
  if (absorption_species.empty()) return;

  for (auto &species_tags : absorption_species) {
    if (species_tags.Species() != "AIR"_spec) ++keys[species_tags.Species()];
  }
}

void keysSpecies(std::unordered_map<SpeciesEnum, Index> &keys,
                 const ArrayOfCIARecord &absorption_cia_data) {
  if (absorption_cia_data.empty()) return;

  for (auto &cia_record : absorption_cia_data) {
    const auto [spec1, spec2] = cia_record.TwoSpecies();
    if (spec1 != "AIR"_spec) ++keys[spec1];
    if (spec2 != "AIR"_spec) ++keys[spec2];
  }
}

void keysSpecies(std::unordered_map<SpeciesEnum, Index> &keys,
                 const AbsorptionLookupTables &absorption_lookup_table) {
  if (absorption_lookup_table.empty()) return;

  for (auto &&spec : absorption_lookup_table | stdv::keys) {
    if (spec != "AIR"_spec) ++keys[spec];
  }
}

void keysSpecies(std::unordered_map<SpeciesEnum, Index> &keys,
                 const ArrayOfXsecRecord &absorption_xsec_fit_data) {
  if (absorption_xsec_fit_data.empty()) return;

  for (auto &xsec_record : absorption_xsec_fit_data) {
    if (xsec_record.Species() != "AIR"_spec) ++keys[xsec_record.Species()];
  }
}

void keysSpecies(std::unordered_map<SpeciesEnum, Index> &keys,
                 const PredefinedModelData &absorption_predefined_model_data) {
  if (absorption_predefined_model_data.empty()) return;

  for (auto &predef_record : absorption_predefined_model_data) {
    if (predef_record.first.spec != "AIR"_spec)
      ++keys[predef_record.first.spec];
  }
}
}  // namespace

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
      not atmospheric_field.contains(p) and
          not static_cast<bool>(allow_missing_pressure),
      "Pressure is missing from the read atmospheric field at \"{}\"",
      basename)

  ARTS_USER_ERROR_IF(
      not atmospheric_field.contains(t) and
          not static_cast<bool>(allow_missing_temperature),
      "Temperature is missing from the read atmospheric field at \"{}\"",
      basename)

  switch (to<MissingFieldComponentError>(deal_with_field_component)) {
    case MissingFieldComponentError::Throw:
      if (atmospheric_field.contains(wind_u) or
          atmospheric_field.contains(wind_v) or
          atmospheric_field.contains(wind_w)) {
        ARTS_USER_ERROR_IF(
            not atmospheric_field.contains(wind_u) or
                not atmospheric_field.contains(wind_v) or
                not atmospheric_field.contains(wind_w),
            "Need all wind components, has [u: {}, v: {}, w: {}]",
            atmospheric_field.contains(wind_u),
            atmospheric_field.contains(wind_v),
            atmospheric_field.contains(wind_w));
      }

      if (atmospheric_field.contains(mag_u) or
          atmospheric_field.contains(mag_v) or
          atmospheric_field.contains(mag_w)) {
        ARTS_USER_ERROR_IF(not atmospheric_field.contains(mag_u) or
                               not atmospheric_field.contains(mag_v) or
                               not atmospheric_field.contains(mag_w),
                           "Need all mag components, has [u: {}, v: {}, w: {}]",
                           atmospheric_field.contains(mag_u),
                           atmospheric_field.contains(mag_v),
                           atmospheric_field.contains(mag_w));
      }
      break;
    case MissingFieldComponentError::Zero:
      if (atmospheric_field.contains(wind_u) or
          atmospheric_field.contains(wind_v) or
          atmospheric_field.contains(wind_w)) {
        if (not atmospheric_field.contains(wind_u))
          atmospheric_field[wind_u] = 0.0;
        if (not atmospheric_field.contains(wind_v))
          atmospheric_field[wind_v] = 0.0;
        if (not atmospheric_field.contains(wind_w))
          atmospheric_field[wind_w] = 0.0;
      }

      if (atmospheric_field.contains(mag_u) or
          atmospheric_field.contains(mag_v) or
          atmospheric_field.contains(mag_w)) {
        if (not atmospheric_field.contains(mag_u))
          atmospheric_field[mag_u] = 0.0;
        if (not atmospheric_field.contains(mag_v))
          atmospheric_field[mag_v] = 0.0;
        if (not atmospheric_field.contains(mag_w))
          atmospheric_field[mag_w] = 0.0;
      }
      break;

    case MissingFieldComponentError::Ignore: break;
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
    ++keys[key.isot];
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

void atmospheric_fieldAppendLineLevelData(
    AtmField &atmospheric_field,
    const AbsorptionBands &absorption_bands,
    const String &basename,
    const String &extrapolation,
    const Index &missing_is_zero,
    const Index &replace_existing) {
  ARTS_TIME_REPORT

  std::unordered_map<QuantumLevelIdentifier, Index> keys;
  keysNLTE(keys, absorption_bands);

  for (auto &[key, value] : absorption_bands) {
    ++keys[key.upper()];
    ++keys[key.lower()];
  }

  append_data(
      atmospheric_field,
      basename,
      extrapolation,
      missing_is_zero,
      replace_existing,
      0,
      keys,
      [](const QuantumLevelIdentifier &x) { return std::format("{}", x); });
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

void atmospheric_fieldAppendPredefSpeciesData(
    AtmField &atmospheric_field,
    const PredefinedModelData &absorption_predefined_model_data,
    const String &basename,
    const String &extrapolation,
    const Index &missing_is_zero,
    const Index &replace_existing) {
  ARTS_TIME_REPORT

  if (absorption_predefined_model_data.empty()) return;

  const auto to_string = [](const SpeciesEnum &x) {
    return String{toString<1>(x)};
  };

  std::unordered_map<SpeciesEnum, Index> keys;
  keysSpecies(keys, absorption_predefined_model_data);

  append_data(atmospheric_field,
              basename,
              extrapolation,
              missing_is_zero,
              replace_existing,
              0,
              keys,
              to_string);

  // H2O might be used, so we will maybe read it
  constexpr auto water = "H2O"_spec;
  if (not keys.contains(water)) {
    keys = {{water, 1}};

    append_data(atmospheric_field,
                basename,
                extrapolation,
                missing_is_zero,
                replace_existing,
                1,
                keys,
                to_string);
  }
}

void atmospheric_fieldAppendAuto(const Workspace &ws,
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

  const NumericTernaryOperator magu{
      .f = NumericTernary{Atm::IGRF13(time, FieldComponent::u)}};
  const NumericTernaryOperator magv{
      .f = NumericTernary{Atm::IGRF13(time, FieldComponent::v)}};
  const NumericTernaryOperator magw{
      .f = NumericTernary{Atm::IGRF13(time, FieldComponent::w)}};

  atmospheric_field[AtmKey::mag_u] = magu;
  atmospheric_field[AtmKey::mag_v] = magv;
  atmospheric_field[AtmKey::mag_w] = magw;
}

void atmospheric_fieldSchmidthFieldFromIGRF(AtmField &atmospheric_field,
                                            const Time &time) {
  ARTS_TIME_REPORT

  const auto magu{Atm::IGRF13(time, FieldComponent::u)};
  const auto magv{Atm::IGRF13(time, FieldComponent::v)};
  const auto magw{Atm::IGRF13(time, FieldComponent::w)};

  atmospheric_field[AtmKey::mag_u] = NumericTernaryOperator{.f = from(magu)};
  atmospheric_field[AtmKey::mag_v] = NumericTernaryOperator{.f = from(magv)};
  atmospheric_field[AtmKey::mag_w] = NumericTernaryOperator{.f = from(magw)};
}

namespace {
Atm::MagnitudeField to_magnitude_field(const GeodeticField3 &u,
                                       const GeodeticField3 &v,
                                       const GeodeticField3 &w) {
  ARTS_USER_ERROR_IF(u.shape() != v.shape() or u.shape() != w.shape(),
                     R"(The field components must have the same shape

u.shape(): {:B,}
v.shape(): {:B,}
w.shape(): {:B,}                  )",
                     u.shape(),
                     v.shape(),
                     w.shape())

  Tensor3 Mag(u.data);
  Tensor3 Theta(u.data);
  Tensor3 Phi(u.data);
  for (Index i = 0; i < Mag.npages(); i++) {
    for (Index j = 0; j < Mag.nrows(); j++) {
      for (Index k = 0; k < Mag.ncols(); k++) {
        const Vector3 x = ecef2geocentric({u[i, j, k], v[i, j, k], w[i, j, k]});
        Mag[i, j, k]    = x[0];
        Theta[i, j, k]  = x[1];
        Phi[i, j, k]    = x[2];
      }
    }
  }

  Atm::MagnitudeField field;

  field.magnitude.data       = std::move(Mag);
  field.theta.data           = std::move(Theta);
  field.phi.data             = std::move(Phi);
  field.magnitude.grids      = u.grids;
  field.theta.grids          = u.grids;
  field.phi.grids            = u.grids;
  field.magnitude.grid_names = u.grid_names;
  field.theta.grid_names     = u.grid_names;
  field.phi.grid_names       = u.grid_names;
  field.magnitude.data_name  = "magnitude";
  field.theta.data_name      = "theta";
  field.phi.data_name        = "phi";

  return field;
}
}  // namespace

void atmospheric_fieldAbsoluteMagneticField(AtmField &atmospheric_field) try {
  ARTS_TIME_REPORT

  constexpr AtmKey u = AtmKey::mag_u;
  constexpr AtmKey v = AtmKey::mag_v;
  constexpr AtmKey w = AtmKey::mag_w;

  Atm::MagnitudeField field =
      to_magnitude_field(atmospheric_field[u].get<GeodeticField3>(),
                         atmospheric_field[v].get<GeodeticField3>(),
                         atmospheric_field[w].get<GeodeticField3>());

  field.component      = FieldComponent::u;
  atmospheric_field[u] = AtmFunctionalData{.f = field};

  field.component      = FieldComponent::v;
  atmospheric_field[v] = AtmFunctionalData{.f = field};

  field.component      = FieldComponent::w;
  atmospheric_field[w] = AtmFunctionalData{.f = std::move(field)};
}
ARTS_METHOD_ERROR_CATCH

void atmospheric_fieldAbsoluteWindField(AtmField &atmospheric_field) try {
  ARTS_TIME_REPORT

  constexpr AtmKey u = AtmKey::wind_u;
  constexpr AtmKey v = AtmKey::wind_v;
  constexpr AtmKey w = AtmKey::wind_w;

  Atm::MagnitudeField field =
      to_magnitude_field(atmospheric_field[u].get<GeodeticField3>(),
                         atmospheric_field[v].get<GeodeticField3>(),
                         atmospheric_field[w].get<GeodeticField3>());

  field.component      = FieldComponent::u;
  atmospheric_field[u] = AtmFunctionalData{.f = field};

  field.component      = FieldComponent::v;
  atmospheric_field[v] = AtmFunctionalData{.f = field};

  field.component      = FieldComponent::w;
  atmospheric_field[w] = AtmFunctionalData{.f = std::move(field)};
}
ARTS_METHOD_ERROR_CATCH

void atmospheric_fieldHydrostaticPressure(
    AtmField &atmospheric_field,
    const NumericTernaryOperator &gravity_operator,
    const AscendingGrid &alts,
    const GeodeticField2 &p0,
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
      not atmospheric_field.contains(AtmKey::t),
      "Must have a temperature field to call this workspace method with a single Numeric reference pressure, so that latitude and longitude grids can be extracted")

  const auto &t = atmospheric_field[AtmKey::t];

  ARTS_USER_ERROR_IF(
      not std::holds_alternative<GeodeticField3>(t.data),
      "Temperature field must be a GeodeticField3 to call this workspace method with a single Numeric reference pressure, so that latitude and longitude grids can be extracted")

  const auto &t0 = std::get<GeodeticField3>(t.data);
  const GeodeticField2 p0_field{
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
