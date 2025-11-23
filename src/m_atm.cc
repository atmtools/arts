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

void atm_fieldInit(AtmField &atm_field,
                   const Numeric &top_of_atmosphere,
                   const String &default_isotopologue) {
  ARTS_TIME_REPORT

  atm_field = AtmField{to<IsoRatioOption>(default_isotopologue)};
  atm_field.top_of_atmosphere = top_of_atmosphere;
}

void atm_pointInit(AtmPoint &atm_point, const String &default_isotopologue) {
  ARTS_TIME_REPORT

  atm_point = AtmPoint{to<IsoRatioOption>(default_isotopologue)};
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
    AtmField &atm_field,
    const String &basename,
    const String &extrapolation,
    const Index &missing_is_zero,
    const Index &replace_existing_,
    const Index &ignore_missing_,
    const std::unordered_map<T, Index> &keys,
    const auto &to_string = [](const auto &x) { return std::format("{}", x); })
  requires(
      std::same_as<typename std::decay_t<decltype(to_string(T{}))>, String>)
{
  const String my_base      = complete_basename(basename);
  const bool ignore_missing = static_cast<bool>(ignore_missing_);
  const bool replace        = static_cast<bool>(replace_existing_);

  Atm::Data x;
  x.alt_low = x.alt_upp = x.lat_low = x.lat_upp = x.lon_low = x.lon_upp =
      to<InterpolationExtrapolation>(extrapolation);

  std::vector<std::string> error{};

  for (const auto &[key, _] : keys) {
    String filename = std::format("{}{}.xml", my_base, to_string(key));

    if (replace or not atm_field.contains(key)) {
      if (find_xml_file_existence(filename)) {
        try {
          xml_extend_from_file(filename, atm_field);
          const bool has = atm_field.contains(key);
          if (not has and not ignore_missing) goto error;
        } catch (...) {
          x.data         = xml_read_from_file_variant(x.data, filename);
          atm_field[key] = x;
        }
      } else if (static_cast<bool>(missing_is_zero)) {
        x.data         = 0.0;
        atm_field[key] = x;
      } else if (ignore_missing) {
      } else {
      error:
        if (error.empty()) error.emplace_back("Missing species / files:");
        error.emplace_back(std::format("  {}", filename));
      }
    }
  }

  ARTS_USER_ERROR_IF(not error.empty(), "{:n}", error)
}

void keysSpecies(std::unordered_map<SpeciesEnum, Index> &keys,
                 const AbsorptionBands &abs_bands) {
  if (abs_bands.empty()) return;

  for (auto &[key, value] : abs_bands) {
    ++keys[key.isot.spec];

    for (auto &line : value.lines) {
      for (auto &ls : line.ls.single_models) {
        if (ls.first != "AIR"_spec) ++keys[ls.first];
      }
    }
  }
}

void keysIsotopologue(std::unordered_map<SpeciesIsotope, Index> &keys,
                      const AbsorptionBands &abs_bands) {
  if (abs_bands.empty()) return;

  for (auto &[key, value] : abs_bands) {
    ++keys[key.isot];
  }
}

void keysNLTE(std::unordered_map<QuantumLevelIdentifier, Index> &keys,
              const AbsorptionBands &abs_bands) {
  if (abs_bands.empty()) return;

  for (auto &[key, value] : abs_bands) {
    ++keys[key.upper()];
    ++keys[key.lower()];
  }
}

void keysSpecies(std::unordered_map<SpeciesEnum, Index> &keys,
                 const ArrayOfSpeciesTag &abs_species) {
  if (abs_species.empty()) return;

  for (auto &species_tags : abs_species) {
    if (species_tags.Spec() != "AIR"_spec) ++keys[species_tags.Spec()];
  }
}

void keysSpecies(std::unordered_map<SpeciesEnum, Index> &keys,
                 const CIARecords &abs_cia_data) {
  if (abs_cia_data.empty()) return;

  for (const auto &[species_pair, cia_record] : abs_cia_data) {
    const auto spec1 = species_pair.spec1;
    const auto spec2 = species_pair.spec2;
    if (spec1 != "AIR"_spec) ++keys[spec1];
    if (spec2 != "AIR"_spec) ++keys[spec2];
  }
}

void keysSpecies(std::unordered_map<SpeciesEnum, Index> &keys,
                 const AbsorptionLookupTables &abs_lookup_data) {
  if (abs_lookup_data.empty()) return;

  for (auto &&spec : abs_lookup_data | stdv::keys) {
    if (spec != "AIR"_spec) ++keys[spec];
  }
}

void keysSpecies(std::unordered_map<SpeciesEnum, Index> &keys,
                 const XsecRecords &abs_xfit_data) {
  if (abs_xfit_data.empty()) return;

  for (auto &[species, xsec_record] : abs_xfit_data) {
    if (species != "AIR"_spec) ++keys[species];
  }
}

void keysSpecies(std::unordered_map<SpeciesEnum, Index> &keys,
                 const PredefinedModelData &abs_predef_data) {
  if (abs_predef_data.empty()) return;

  for (auto &predef_record : abs_predef_data) {
    if (predef_record.first.spec != "AIR"_spec)
      ++keys[predef_record.first.spec];
  }
}
}  // namespace

void atm_fieldAppendBaseData(AtmField &atm_field,
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

  append_data(atm_field,
              basename,
              extrapolation,
              0,
              replace_existing,
              1,
              keys,
              [](const AtmKey &x) { return std::format("{}", x); });

  using enum AtmKey;

  ARTS_USER_ERROR_IF(
      not atm_field.contains(p) and
          not static_cast<bool>(allow_missing_pressure),
      "Pressure is missing from the read atmospheric field at \"{}\"",
      basename)

  ARTS_USER_ERROR_IF(
      not atm_field.contains(t) and
          not static_cast<bool>(allow_missing_temperature),
      "Temperature is missing from the read atmospheric field at \"{}\"",
      basename)

  switch (to<MissingFieldComponentError>(deal_with_field_component)) {
    case MissingFieldComponentError::Throw:
      if (atm_field.contains(wind_u) or atm_field.contains(wind_v) or
          atm_field.contains(wind_w)) {
        ARTS_USER_ERROR_IF(
            not atm_field.contains(wind_u) or not atm_field.contains(wind_v) or
                not atm_field.contains(wind_w),
            "Need all wind components, has [u: {}, v: {}, w: {}]",
            atm_field.contains(wind_u),
            atm_field.contains(wind_v),
            atm_field.contains(wind_w));
      }

      if (atm_field.contains(mag_u) or atm_field.contains(mag_v) or
          atm_field.contains(mag_w)) {
        ARTS_USER_ERROR_IF(not atm_field.contains(mag_u) or
                               not atm_field.contains(mag_v) or
                               not atm_field.contains(mag_w),
                           "Need all mag components, has [u: {}, v: {}, w: {}]",
                           atm_field.contains(mag_u),
                           atm_field.contains(mag_v),
                           atm_field.contains(mag_w));
      }
      break;
    case MissingFieldComponentError::Zero:
      if (atm_field.contains(wind_u) or atm_field.contains(wind_v) or
          atm_field.contains(wind_w)) {
        if (not atm_field.contains(wind_u)) atm_field[wind_u] = 0.0;
        if (not atm_field.contains(wind_v)) atm_field[wind_v] = 0.0;
        if (not atm_field.contains(wind_w)) atm_field[wind_w] = 0.0;
      }

      if (atm_field.contains(mag_u) or atm_field.contains(mag_v) or
          atm_field.contains(mag_w)) {
        if (not atm_field.contains(mag_u)) atm_field[mag_u] = 0.0;
        if (not atm_field.contains(mag_v)) atm_field[mag_v] = 0.0;
        if (not atm_field.contains(mag_w)) atm_field[mag_w] = 0.0;
      }
      break;

    case MissingFieldComponentError::Ignore: break;
  }
}

void atm_fieldAppendLineSpeciesData(AtmField &atm_field,
                                    const AbsorptionBands &abs_bands,
                                    const String &basename,
                                    const String &extrapolation,
                                    const Index &missing_is_zero,
                                    const Index &replace_existing) {
  ARTS_TIME_REPORT

  std::unordered_map<SpeciesEnum, Index> keys;
  keysSpecies(keys, abs_bands);

  append_data(atm_field,
              basename,
              extrapolation,
              missing_is_zero,
              replace_existing,
              0,
              keys,
              [](const SpeciesEnum &x) { return String{toString<1>(x)}; });
}

void atm_fieldAppendLineIsotopologueData(AtmField &atm_field,
                                         const AbsorptionBands &abs_bands,
                                         const String &basename,
                                         const String &extrapolation,
                                         const Index &missing_is_zero,
                                         const Index &replace_existing) {
  ARTS_TIME_REPORT

  std::unordered_map<SpeciesIsotope, Index> keys;
  keysIsotopologue(keys, abs_bands);

  for (auto &[key, value] : abs_bands) {
    ++keys[key.isot];
  }

  append_data(atm_field,
              basename,
              extrapolation,
              missing_is_zero,
              replace_existing,
              0,
              keys,
              [](const SpeciesIsotope &x) { return x.FullName(); });
}

void atm_fieldAppendLineLevelData(AtmField &atm_field,
                                  const AbsorptionBands &abs_bands,
                                  const String &basename,
                                  const String &extrapolation,
                                  const Index &missing_is_zero,
                                  const Index &replace_existing) {
  ARTS_TIME_REPORT

  std::unordered_map<QuantumLevelIdentifier, Index> keys;
  keysNLTE(keys, abs_bands);

  for (auto &[key, value] : abs_bands) {
    ++keys[key.upper()];
    ++keys[key.lower()];
  }

  append_data(
      atm_field,
      basename,
      extrapolation,
      missing_is_zero,
      replace_existing,
      0,
      keys,
      [](const QuantumLevelIdentifier &x) { return std::format("{}", x); });
}

void atm_fieldAppendTagsSpeciesData(AtmField &atm_field,
                                    const ArrayOfSpeciesTag &abs_species,
                                    const String &basename,
                                    const String &extrapolation,
                                    const Index &missing_is_zero,
                                    const Index &replace_existing) {
  ARTS_TIME_REPORT

  std::unordered_map<SpeciesEnum, Index> keys;
  keysSpecies(keys, abs_species);

  append_data(atm_field,
              basename,
              extrapolation,
              missing_is_zero,
              replace_existing,
              0,
              keys,
              [](const SpeciesEnum &x) { return String{toString<1>(x)}; });
}

void atm_fieldAppendCIASpeciesData(AtmField &atm_field,
                                   const CIARecords &abs_cia_data,
                                   const String &basename,
                                   const String &extrapolation,
                                   const Index &missing_is_zero,
                                   const Index &replace_existing) {
  ARTS_TIME_REPORT

  std::unordered_map<SpeciesEnum, Index> keys;
  keysSpecies(keys, abs_cia_data);

  append_data(atm_field,
              basename,
              extrapolation,
              missing_is_zero,
              replace_existing,
              0,
              keys,
              [](const SpeciesEnum &x) { return String{toString<1>(x)}; });
}

void atm_fieldAppendLookupTableSpeciesData(
    AtmField &atm_field,
    const AbsorptionLookupTables &abs_lookup_data,
    const String &basename,
    const String &extrapolation,
    const Index &missing_is_zero,
    const Index &replace_existing) {
  ARTS_TIME_REPORT

  std::unordered_map<SpeciesEnum, Index> keys;
  keysSpecies(keys, abs_lookup_data);

  append_data(atm_field,
              basename,
              extrapolation,
              missing_is_zero,
              replace_existing,
              0,
              keys,
              [](const SpeciesEnum &x) { return String{toString<1>(x)}; });
}

void atm_fieldAppendXsecSpeciesData(AtmField &atm_field,
                                    const XsecRecords &abs_xfit_data,
                                    const String &basename,
                                    const String &extrapolation,
                                    const Index &missing_is_zero,
                                    const Index &replace_existing) {
  ARTS_TIME_REPORT

  std::unordered_map<SpeciesEnum, Index> keys;
  keysSpecies(keys, abs_xfit_data);

  append_data(atm_field,
              basename,
              extrapolation,
              missing_is_zero,
              replace_existing,
              0,
              keys,
              [](const SpeciesEnum &x) { return String{toString<1>(x)}; });
}

void atm_fieldAppendPredefSpeciesData(
    AtmField &atm_field,
    const PredefinedModelData &abs_predef_data,
    const String &basename,
    const String &extrapolation,
    const Index &missing_is_zero,
    const Index &replace_existing) {
  ARTS_TIME_REPORT

  if (abs_predef_data.empty()) return;

  const auto to_string = [](const SpeciesEnum &x) {
    return String{toString<1>(x)};
  };

  std::unordered_map<SpeciesEnum, Index> keys;
  keysSpecies(keys, abs_predef_data);

  append_data(atm_field,
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

    append_data(atm_field,
                basename,
                extrapolation,
                missing_is_zero,
                replace_existing,
                1,
                keys,
                to_string);
  }
}

void atm_fieldAppendAuto(const Workspace &ws,
                         AtmField &atm_field,
                         const String &basename,
                         const String &extrapolation,
                         const Index &missing_is_zero,
                         const Index &replace_existing,
                         const Index &load_isot,
                         const Index &load_nlte) {
  ARTS_TIME_REPORT

  if (const String lines_str = "abs_bands"; ws.wsv_and_contains(lines_str)) {
    using lines_t         = AbsorptionBands;
    const auto &abs_bands = ws.get<lines_t>(lines_str);

    atm_fieldAppendLineSpeciesData(atm_field,
                                   abs_bands,
                                   basename,
                                   extrapolation,
                                   missing_is_zero,
                                   replace_existing);

    if (static_cast<bool>(load_isot)) {
      atm_fieldAppendLineIsotopologueData(atm_field,
                                          abs_bands,
                                          basename,
                                          extrapolation,
                                          missing_is_zero,
                                          replace_existing);
    }

    if (static_cast<bool>(load_nlte)) {
      atm_fieldAppendLineLevelData(atm_field,
                                   abs_bands,
                                   basename,
                                   extrapolation,
                                   missing_is_zero,
                                   replace_existing);
    }
  }

  if (const String cia_str = "abs_cia_data"; ws.wsv_and_contains(cia_str)) {
    using cia_t              = CIARecords;
    const auto &abs_cia_data = ws.get<cia_t>(cia_str);
    atm_fieldAppendCIASpeciesData(atm_field,
                                  abs_cia_data,
                                  basename,
                                  extrapolation,
                                  missing_is_zero,
                                  replace_existing);
  }

  if (const String lookup_str = "abs_lookup_data";
      ws.wsv_and_contains(lookup_str)) {
    using lookup_t              = AbsorptionLookupTables;
    const auto &abs_lookup_data = ws.get<lookup_t>(lookup_str);
    atm_fieldAppendLookupTableSpeciesData(atm_field,
                                          abs_lookup_data,
                                          basename,
                                          extrapolation,
                                          missing_is_zero,
                                          replace_existing);
  }

  if (const String xsec_str = "abs_xfit_data"; ws.wsv_and_contains(xsec_str)) {
    using xsec_t              = XsecRecords;
    const auto &abs_xfit_data = ws.get<xsec_t>(xsec_str);
    atm_fieldAppendXsecSpeciesData(atm_field,
                                   abs_xfit_data,
                                   basename,
                                   extrapolation,
                                   missing_is_zero,
                                   replace_existing);
  }

  if (const String predef_str = "abs_predef_data";
      ws.wsv_and_contains(predef_str)) {
    using predef_t              = PredefinedModelData;
    const auto &abs_predef_data = ws.get<predef_t>(predef_str);
    atm_fieldAppendPredefSpeciesData(atm_field,
                                     abs_predef_data,
                                     basename,
                                     extrapolation,
                                     missing_is_zero,
                                     replace_existing);
  }

  if (const String species_str = "abs_species";
      ws.wsv_and_contains(species_str)) {
    using aospec_t          = ArrayOfSpeciesTag;
    const auto &abs_species = ws.get<aospec_t>(species_str);
    atm_fieldAppendTagsSpeciesData(atm_field,
                                   abs_species,
                                   basename,
                                   extrapolation,
                                   missing_is_zero,
                                   replace_existing);
  }
}

void atm_fieldIGRF(AtmField &atm_field, const Time &time) {
  ARTS_TIME_REPORT

  const NumericTernaryOperator magu{
      .f = NumericTernary{Atm::IGRF13(time, FieldComponent::u)}};
  const NumericTernaryOperator magv{
      .f = NumericTernary{Atm::IGRF13(time, FieldComponent::v)}};
  const NumericTernaryOperator magw{
      .f = NumericTernary{Atm::IGRF13(time, FieldComponent::w)}};

  atm_field[AtmKey::mag_u] = magu;
  atm_field[AtmKey::mag_v] = magv;
  atm_field[AtmKey::mag_w] = magw;
}

void atm_fieldSchmidthFieldFromIGRF(AtmField &atm_field, const Time &time) {
  ARTS_TIME_REPORT

  const auto magu{Atm::IGRF13(time, FieldComponent::u)};
  const auto magv{Atm::IGRF13(time, FieldComponent::v)};
  const auto magw{Atm::IGRF13(time, FieldComponent::w)};

  atm_field[AtmKey::mag_u] = NumericTernaryOperator{.f = from(magu)};
  atm_field[AtmKey::mag_v] = NumericTernaryOperator{.f = from(magv)};
  atm_field[AtmKey::mag_w] = NumericTernaryOperator{.f = from(magw)};
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

void atm_fieldAbsoluteMagneticField(AtmField &atm_field) try {
  ARTS_TIME_REPORT

  constexpr AtmKey u = AtmKey::mag_u;
  constexpr AtmKey v = AtmKey::mag_v;
  constexpr AtmKey w = AtmKey::mag_w;

  Atm::MagnitudeField field =
      to_magnitude_field(atm_field[u].get<GeodeticField3>(),
                         atm_field[v].get<GeodeticField3>(),
                         atm_field[w].get<GeodeticField3>());

  field.component = FieldComponent::u;
  atm_field[u]    = AtmFunctionalData{.f = field};

  field.component = FieldComponent::v;
  atm_field[v]    = AtmFunctionalData{.f = field};

  field.component = FieldComponent::w;
  atm_field[w]    = AtmFunctionalData{.f = std::move(field)};
}
ARTS_METHOD_ERROR_CATCH

void atm_fieldAbsoluteWindField(AtmField &atm_field) try {
  ARTS_TIME_REPORT

  constexpr AtmKey u = AtmKey::wind_u;
  constexpr AtmKey v = AtmKey::wind_v;
  constexpr AtmKey w = AtmKey::wind_w;

  Atm::MagnitudeField field =
      to_magnitude_field(atm_field[u].get<GeodeticField3>(),
                         atm_field[v].get<GeodeticField3>(),
                         atm_field[w].get<GeodeticField3>());

  field.component = FieldComponent::u;
  atm_field[u]    = AtmFunctionalData{.f = field};

  field.component = FieldComponent::v;
  atm_field[v]    = AtmFunctionalData{.f = field};

  field.component = FieldComponent::w;
  atm_field[w]    = AtmFunctionalData{.f = std::move(field)};
}
ARTS_METHOD_ERROR_CATCH

void atm_fieldHydrostaticPressure(
    AtmField &atm_field,
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
      not has_def_t and not atm_field.contains(AtmKey::t),
      "atm_field lacks temperature and no default temperature given")

  ARTS_USER_ERROR_IF(
      not has_def_r and atm_field.nspec() == 0,
      "atm_field lacks species and no default specific gas constant given")

  const Tensor3 scale_factor = [&]() {
    Tensor3 scl(nalt, nlat, nlon);
    for (Index i = 0; i < nalt; i++) {
      for (Index j = 0; j < nlat; j++) {
        for (Index k = 0; k < nlon; k++) {
          const Numeric al = alts[i];
          const Numeric la = lats[j];
          const Numeric lo = lons[k];

          const Numeric g = gravity_operator(al, la, lo);
          const AtmPoint atm_point{atm_field.at(al, la, lo)};

          const Numeric inv_specific_gas_constant =
              has_def_r ? 1.0 / fixed_specific_gas_constant
                        : (1e-3 * atm_point.mean_mass() / Constant::R);
          const Numeric inv_temp = has_def_t ? 1.0 / fixed_atm_temperature
                                             : 1.0 / atm_point.temperature;

          // Partial rho, no pressure
          scl[i, j, k] = g * inv_specific_gas_constant * inv_temp;
        }
      }
    }

    return scl;
  }();

  atm_field[AtmKey::p] = Atm::FunctionalData{Atm::HydrostaticPressure(
      scale_factor,
      p0,
      alts,
      to<HydrostaticPressureOption>(hydrostatic_option))};
}

void atm_fieldHydrostaticPressure(
    AtmField &atm_field,
    const NumericTernaryOperator &gravity_operator,
    const AscendingGrid &alts,
    const Numeric &p0,
    const Numeric &fixed_specific_gas_constant,
    const Numeric &fixed_atm_temperature,
    const String &hydrostatic_option) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      not atm_field.contains(AtmKey::t),
      "Must have a temperature field to call this workspace method with a single Numeric reference pressure, so that latitude and longitude grids can be extracted")

  const auto &t = atm_field[AtmKey::t];

  ARTS_USER_ERROR_IF(
      not std::holds_alternative<GeodeticField3>(t.data),
      "Temperature field must be a GeodeticField3 to call this workspace method with a single Numeric reference pressure, so that latitude and longitude grids can be extracted")

  const auto &t0 = std::get<GeodeticField3>(t.data);
  const GeodeticField2 p0_field{
      .data_name  = "Pressure",
      .data       = Matrix(t0.grid<1>().size(), t0.grid<2>().size(), p0),
      .grid_names = {t0.gridname<1>(), t0.gridname<2>()},
      .grids      = {t0.grid<1>(), t0.grid<2>()}};

  atm_fieldHydrostaticPressure(atm_field,
                               gravity_operator,
                               alts,
                               p0_field,
                               fixed_specific_gas_constant,
                               fixed_atm_temperature,
                               hydrostatic_option);
}

void atm_fieldFromProfile(
    AtmField &atm_field,
    const ArrayOfAtmPoint &atm_profile,
    const AscendingGrid &alt_grid,
    const InterpolationExtrapolation &altitude_extrapolation) {
  ARTS_USER_ERROR_IF(not arr::same_size(atm_profile, alt_grid),
                     "Mismatch in size between atm_profile and alt_grid");
  ARTS_USER_ERROR_IF(alt_grid.empty(), "Empty altitude grid")

  atm_field = Atm::atm_from_profile(
      atm_profile, alt_grid, altitude_extrapolation, alt_grid.back());
}
