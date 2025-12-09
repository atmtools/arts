#include <arts_omp.h>
#include <configtypes.h>
#include <debug.h>
#include <enumsAbsorptionBandSortingOption.h>
#include <enumsLineByLineCutoffType.h>
#include <isotopologues.h>
#include <jacobian.h>
#include <lbl.h>
#include <lbl_data.h>
#include <lbl_lineshape_linemixing.h>
#include <minimize.h>
#include <partfun.h>
#include <path_point.h>
#include <quantum.h>
#include <rtepack.h>
#include <sorting.h>
#include <species_tags.h>
#include <time_report.h>
#include <workspace.h>
#include <xml_io.h>
#include <xml_io_old.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iterator>
#include <ranges>
#include <unordered_map>

void abs_bandsSelectFrequencyByLine(AbsorptionBands& abs_bands,
                                    const Numeric& fmin,
                                    const Numeric& fmax) try {
  ARTS_TIME_REPORT

  std::vector<QuantumIdentifier> to_remove;

  for (auto& [key, band] : abs_bands) {
    band.sort();

    auto& lines = band.lines;

    lines.erase(std::remove_if(lines.begin(),
                               lines.end(),
                               [fmin, fmax](const lbl::line& l) {
                                 return l.f0 < fmin or l.f0 > fmax;
                               }),
                lines.end());

    if (lines.empty()) to_remove.push_back(key);
  }

  for (const auto& key : to_remove) abs_bands.erase(key);
}
ARTS_METHOD_ERROR_CATCH

void abs_bandsSelectFrequencyByBand(AbsorptionBands& abs_bands,
                                    const Numeric& fmin,
                                    const Numeric& fmax) try {
  ARTS_TIME_REPORT

  std::vector<QuantumIdentifier> to_remove;

  for (auto& [key, band] : abs_bands) {
    band.sort();

    const bool criteria = band.lines.empty() or band.lines.front().f0 > fmax or
                          band.lines.back().f0 < fmin;

    if (criteria) to_remove.push_back(key);
  }

  for (const auto& key : to_remove) abs_bands.erase(key);
}
ARTS_METHOD_ERROR_CATCH

void abs_bandsKeepID(AbsorptionBands& abs_bands,
                     const QuantumIdentifier& id,
                     const Index& line) try {
  ARTS_TIME_REPORT

  if (auto ptr = abs_bands.find(id); ptr != abs_bands.end()) {
    const auto& [key, band] = *ptr;

    QuantumIdentifier newk = key;
    AbsorptionBand newb    = band;
    abs_bands              = {};
    AbsorptionBand& data   = abs_bands[newk];
    data                   = std::move(newb);

    if (line >= 0) {
      ARTS_USER_ERROR_IF(static_cast<Size>(line) >= band.lines.size(),
                         "Line index out of range: {}",
                         line)
      data.lines = {data.lines[line]};
    }
  } else {
    abs_bands = {};
  }
}
ARTS_METHOD_ERROR_CATCH

void abs_bandsReadSpeciesSplitCatalog(
    AbsorptionBands& abs_bands,
    const ArrayOfSpeciesTag& absorbtion_species,
    const String& basename,
    const Index& ignore_missing_) try {
  ARTS_TIME_REPORT

  abs_bands.clear();

  const bool ignore_missing = static_cast<bool>(ignore_missing_);

  const String my_base = complete_basename(basename);

  std::set<SpeciesIsotope> isotopologues;
  for (auto& spec : absorbtion_species) {
    if (spec.type == SpeciesTagType::Plain) {
      if (spec.is_joker()) {
        for (auto&& isot : Species::isotopologues(spec.Spec())) {
          if (isot.is_predefined()) continue;
          if (isot.is_joker()) continue;
          isotopologues.insert(isot);
        }
      } else {
        isotopologues.insert(spec.Isotopologue());
      }
    }
  }

  std::vector<std::string> read_errors;
  std::vector<std::string> file_errors;
  std::vector<SpeciesIsotope> visot(isotopologues.begin(), isotopologues.end());
#pragma omp parallel for schedule(dynamic) if (!arts_omp_in_parallel())
  for (std::size_t iisot = 0; iisot < isotopologues.size(); iisot++) {
    try {
      const auto& isot{visot[iisot]};
      String filename{my_base + isot.FullName() + ".xml"};
      if (find_xml_file_existence(filename)) {
        AbsorptionBands other;
        xml_read_from_file(filename, other);
#pragma omp critical(abs_bandsReadSpeciesSplitCatalogInsert)
        abs_bands.insert(std::make_move_iterator(other.begin()),
                         std::make_move_iterator(other.end()));
      } else if (not ignore_missing) {
#pragma omp critical(abs_bandsReadSpeciesSplitCatalogFileError)
        file_errors.push_back(filename);
      }
    } catch (const std::exception& e) {
#pragma omp critical(abs_bandsReadSpeciesSplitCatalogReadError)
      read_errors.emplace_back(e.what());
    }
  }
  if (file_errors.size() and read_errors.size())
    ARTS_USER_ERROR(
        "Files not found:\n{}\n\nReading errors:\n{}", file_errors, read_errors)

  if (file_errors.size()) ARTS_USER_ERROR("Files not found:\n{}", file_errors)
  if (read_errors.size()) ARTS_USER_ERROR("{}", read_errors);
}
ARTS_METHOD_ERROR_CATCH

void abs_bandsReadSplit(AbsorptionBands& abs_bands, const String& dir) try {
  ARTS_TIME_REPORT

  abs_bands = {};

  std::vector<std::filesystem::path> paths;
  stdr::copy_if(std::filesystem::directory_iterator(std::filesystem::path(dir)),
                std::back_inserter(paths),
                [](auto& entry) {
                  return entry.is_regular_file() and
                         entry.path().extension() == ".xml";
                });
  stdr::sort(paths);

  std::string error{};

#pragma omp parallel for schedule(dynamic) if (!arts_omp_in_parallel())
  for (Size i = 0; i < paths.size(); i++) {
    try {
      AbsorptionBands bands;
      xml_read_from_file(paths[i].string(), bands);
#pragma omp critical
      abs_bands.insert(std::make_move_iterator(bands.begin()),
                       std::make_move_iterator(bands.end()));
    } catch (std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  if (not error.empty()) throw std::runtime_error(error);
}
ARTS_METHOD_ERROR_CATCH

void abs_bandsSaveSplit(const AbsorptionBands& abs_bands,
                        const String& dir) try {
  ARTS_TIME_REPORT

  auto create_if_not = [](const std::filesystem::path& path) {
    if (not std::filesystem::exists(path)) {
      std::filesystem::create_directories(path);
    }
    return path;
  };

  const auto p = create_if_not(dir);

  std::unordered_map<SpeciesIsotope, AbsorptionBands> isotopologues_data;
  for (auto& [key, band] : abs_bands) {
    isotopologues_data[key.isot][key] = band;
  }

  for (const auto& [isot, bands] : isotopologues_data) {
    xml_write_to_file(
        (p / std::format("{}.xml", isot)).string(), bands, FileType::ascii, 0);
  }
}
ARTS_METHOD_ERROR_CATCH

void abs_bandsSetZeeman(AbsorptionBands& abs_bands,
                        const SpeciesIsotope& species,
                        const Numeric& fmin,
                        const Numeric& fmax,
                        const Index& _on) try {
  ARTS_TIME_REPORT

  const bool on = static_cast<bool>(_on);

  for (auto& [key, band] : abs_bands) {
    if (key.isot != species) continue;

    for (auto& line : band.lines) {
      if (line.f0 >= fmin and line.f0 <= fmax) {
        line.z.on = on;
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_propmatAddLines(PropmatVector& pm,
                              StokvecVector& sv,
                              PropmatMatrix& dpm,
                              StokvecMatrix& dsv,
                              const AscendingGrid& f_grid,
                              const JacobianTargets& jac_targets,
                              const SpeciesEnum& species,
                              const AbsorptionBands& abs_bands,
                              const LinemixingEcsData& abs_ecs_data,
                              const AtmPoint& atm_point,
                              const PropagationPathPoint& path_point,
                              const Index& no_negative_absorption) try {
  ARTS_TIME_REPORT

  const Size n = arts_omp_get_max_threads();
  if (n == 1 or static_cast<Size>(arts_omp_in_parallel()) or
      n > f_grid.size()) {
    lbl::calculate(pm,
                   sv,
                   dpm,
                   dsv,
                   f_grid,
                   Range(0, f_grid.size()),
                   jac_targets,
                   species,
                   abs_bands,
                   abs_ecs_data,
                   atm_point,
                   path_point.los,
                   no_negative_absorption);
  } else {
    const auto f_ranges = matpack::omp_offset_count(f_grid.size(), n);
    std::string error;
#pragma omp parallel for
    for (Size i = 0; i < n; i++) {
      try {
        lbl::calculate(pm,
                       sv,
                       dpm,
                       dsv,
                       f_grid,
                       f_ranges[i],
                       jac_targets,
                       species,
                       abs_bands,
                       abs_ecs_data,
                       atm_point,
                       path_point.los,
                       no_negative_absorption);
      } catch (std::exception& e) {
#pragma omp critical
        if (error.empty()) error = e.what();
      }
    }

    if (not error.empty()) throw std::runtime_error(error);
  }
}
ARTS_METHOD_ERROR_CATCH

void abs_bandsReadHITRAN(AbsorptionBands& abs_bands,
                         const String& filename,
                         const Vector2& frequency_range,
                         const String& line_strength_option,
                         const Index& compute_zeeman_parameters) try {
  ARTS_TIME_REPORT

  using namespace Quantum;

  const AbsorptionBand default_band{.lines        = {},
                                    .lineshape    = LineByLineLineshape::VP_LTE,
                                    .cutoff       = LineByLineCutoffType::None,
                                    .cutoff_value = NAN};

  const auto selection = to<HitranLineStrengthOption>(line_strength_option);

  const bool do_zeeman = static_cast<bool>(compute_zeeman_parameters);

  const auto data =
      lbl::read_hitran_par(open_input_file(filename), frequency_range);

  abs_bands = {};
  for (auto& line : data) {
    auto [mapped_band, _] = abs_bands.try_emplace(
        global_state(global_types, line.qid), default_band);

    mapped_band->second.lines.emplace_back(
        line.from(selection, local_state(local_types, line.qid), do_zeeman));
  }
}
ARTS_METHOD_ERROR_CATCH

namespace {
template <class Key, class T>
const T& get_value(const Key& k, const std::unordered_map<Key, T>& m) {
  auto it = m.find(k);
  ARTS_USER_ERROR_IF(it == m.end(), "Key not found: {}", k)
  return it->second;
}

template <class Key, class T>
T& get_value(const Key& k, std::unordered_map<Key, T>& m) {
  auto it = m.find(k);
  ARTS_USER_ERROR_IF(it == m.end(), "Key not found: {}", k)
  return it->second;
}
}  // namespace

void abs_bandsLineMixingAdaptation(AbsorptionBands& abs_bands,
                                   const LinemixingEcsData& abs_ecs_data,
                                   const AtmPoint& atm_point,
                                   const AscendingGrid& temperatures,
                                   const QuantumIdentifier& band_key,
                                   const Index& rosenkranz_fit_order,
                                   const Index& polynomial_fit_degree) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(temperatures.empty() or temperatures.front() <= 0.0,
                     "Need a positive temperature grid")

  ARTS_USER_ERROR_IF(rosenkranz_fit_order != 1 and rosenkranz_fit_order != 2,
                     "Only 1 or 2 is supported for the ordered fit")
  ARTS_USER_ERROR_IF(
      polynomial_fit_degree < 1 or
          polynomial_fit_degree > static_cast<Index>(temperatures.size()) - 1,
      "Polynomial degree must be between 1 and the number of temperatures - 1")

  auto& band = get_value(band_key, abs_bands);

  if (band.lines.empty()) return;

  for (auto& line : band.lines | stdv::drop(1)) {
    ARTS_USER_ERROR_IF(
        stdr::any_of(line.ls.single_models | stdv::keys,
                     [&other = band.lines.front().ls.single_models](
                         SpeciesEnum x) { return not other.contains(x); }),
        "Inconsistent line shape models, all lines must have the same broadening species")
  }

  const Size K = band.lines.front().ls.single_models.size();
  const Size M = temperatures.size();
  const Size N = band.lines.size();

  lbl::voigt::ecs::ComputeData com_data({}, atm_point);
  ComplexTensor3 eqv_str(M, K, N), eqv_val(M, K, N);

  lbl::voigt::ecs::equivalent_values(eqv_str,
                                     eqv_val,
                                     com_data,
                                     band_key,
                                     band,
                                     abs_ecs_data.at(band_key.isot),
                                     atm_point,
                                     temperatures);
  using enum LineShapeModelVariable;
  using enum LineShapeModelType;
  using Constant::c, Constant::pi;

  band.sort(LineByLineVariable::f0);

  Matrix lbl_str(M, N);
  for (Size i = 0; i < M; i++) {
    const Numeric Q = PartitionFunctions::Q(temperatures[i], band_key.isot);
    for (Size k = 0; k < N; k++) {
      auto& line    = band.lines[k];
      lbl_str[i, k] = line.s(temperatures[i], Q) * Math::pow2(c) / (8 * pi);
    }
  }

  for (auto& line : band.lines) {
    for (auto& lsm : line.ls.single_models | stdv::values) {
      lsm.remove_variables<Y, G, DV>();
    }
  }

  const ArrayOfSpeciesEnum specs{
      std::from_range, band.lines.front().ls.single_models | stdv::keys};

  ComplexTensor3 lbl_val(M, K, N);
  for (Size i = 0; i < M; i++) {
    for (Size j = 0; j < K; j++) {
      for (Size k = 0; k < N; k++) {
        auto& line       = band.lines[k];
        lbl_val[i, j, k] = Complex{
            line.f0 + line.ls.single_models.at(specs[j]).D0(
                          line.ls.T0, temperatures[i], atm_point.pressure),
            line.ls.single_models.at(specs[j]).G0(
                line.ls.T0, temperatures[i], atm_point.pressure)};
      }
    }
  }

  for (Size i = 0; i < M; i++) {
    for (Size j = 0; j < K; j++) {
      auto s = eqv_str[i, j, joker];
      auto v = eqv_val[i, j, joker];
      bubble_sort_by(
          [&](auto I1, auto I2) { return v[I1].real() > v[I2].real(); }, s, v);
    }
  }

  eqv_val -= lbl_val;

  for (Size j = 0; j < K; j++) eqv_str[joker, j, joker] /= lbl_str;

  eqv_val.real() /= Math::pow2(atm_point.pressure);
  eqv_val.imag() /= Math::pow3(atm_point.pressure);

  eqv_str.real() -= 1.0;
  eqv_str.real() /= Math::pow2(atm_point.pressure);
  eqv_str.imag() /= atm_point.pressure;

  using namespace Minimize;
  for (Size j = 0; j < K; j++) {
    for (Size k = 0; k < N; k++) {
      auto& spec = band.lines[k].ls.single_models.at(specs[j]);

      if (rosenkranz_fit_order >= 1) {
        auto yfit = polyfit(
            temperatures, eqv_str[joker, j, k].imag(), polynomial_fit_degree);
        ARTS_USER_ERROR_IF(
            not yfit, "Cannot fit y for line {} of band {}", k, band_key)
        spec.data[Y] = {POLY, *yfit};
      }

      if (rosenkranz_fit_order >= 2) {
        auto gfit = polyfit(
            temperatures, eqv_str[joker, j, k].real(), polynomial_fit_degree);
        ARTS_USER_ERROR_IF(
            not gfit, "Cannot fit g for line {} of band {}", k, band_key)
        spec.data[G] = {POLY, *gfit};

        auto dfit = polyfit(
            temperatures, eqv_val[joker, j, k].real(), polynomial_fit_degree);
        ARTS_USER_ERROR_IF(
            not dfit, "Cannot fit dv for line {} of band {}", k, band_key)
        spec.data[DV] = {POLY, *dfit};
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void abs_bandsReadSpeciesSplitARTSCAT(
    AbsorptionBands& abs_bands,
    const ArrayOfSpeciesTag& absorbtion_species,
    const String& basename,
    const Index& ignore_missing_,
    const Index& pure_species_) try {
  ARTS_TIME_REPORT

  Array<ArrayOfArtscatMeta> meta;

  const bool ignore_missing = static_cast<bool>(ignore_missing_);
  const bool pure_species   = static_cast<bool>(pure_species_);

  abs_bands.clear();

  const String my_base = complete_basename(basename);

  std::vector<std::string> file_errors;

  if (pure_species) {
    std::set<SpeciesEnum> specieses;
    for (auto& spec : absorbtion_species) {
      if (spec.type == SpeciesTagType::Plain) {
        specieses.insert(spec.Spec());
      }
    }

    std::vector<SpeciesEnum> vspecieses(specieses.begin(), specieses.end());
    for (std::size_t ispec = 0; ispec < specieses.size(); ispec++) {
      const auto& spec{vspecieses[ispec]};
      String filename{my_base + String{toString(spec)} + ".xml"};
      if (find_xml_file_existence(filename)) {
        xml_read_from_file(filename, meta.emplace_back());
      } else if (not ignore_missing) {
        file_errors.push_back(filename);
      }
    }
  } else {
    std::set<SpeciesIsotope> isotopologues;
    for (auto& spec : absorbtion_species) {
      if (spec.type == SpeciesTagType::Plain) {
        if (spec.is_joker()) {
          for (auto&& isot : Species::isotopologues(spec.Spec())) {
            if (isot.is_predefined()) continue;
            if (isot.is_joker()) continue;
            isotopologues.insert(isot);
          }
        } else {
          isotopologues.insert(spec.Isotopologue());
        }
      }
    }

    std::vector<SpeciesIsotope> visot(isotopologues.begin(),
                                      isotopologues.end());
    for (std::size_t iisot = 0; iisot < isotopologues.size(); iisot++) {
      const auto& isot{visot[iisot]};
      String filename{my_base + isot.FullName() + ".xml"};
      if (find_xml_file_existence(filename)) {
        xml_read_from_file(filename, meta.emplace_back());
      } else if (not ignore_missing) {
        file_errors.push_back(filename);
      }
    }
  }

  ARTS_USER_ERROR_IF(
      not file_errors.empty(), "Files not found:\n{:B,}", file_errors)

  for (auto& am : meta) {
    for (auto& m : am) {
      abs_bands[std::move(m.quantumidentity)].emplace_back(std::move(m.data));
    }
  }
}
ARTS_METHOD_ERROR_CATCH

namespace {
void sumup_zeeman(PropmatVectorView spectral_propmat,
                  PropmatMatrixView spectral_propmat_jac,
                  VectorView dispersion,
                  MatrixView dispersion_jac,
                  ComplexVectorView pm,
                  ComplexMatrixView dpm,
                  const AbsorptionBands& abs_bands,
                  const ConstVectorView& freq_grid,
                  const JacobianTargets& jac_targets,
                  const AtmPoint& atm_point,
                  const SpeciesEnum& species,
                  const Index& no_negative_absorption,
                  const Vector3& mag,
                  const Vector2& los,
                  const ZeemanPolarization pol) {
  using namespace lbl::zeeman;

  pm  = 0;
  dpm = 0;
  lbl::voigt::lte::calculate(pm,
                             dpm,
                             abs_bands,
                             freq_grid,
                             jac_targets,
                             atm_point,
                             pol,
                             species,
                             no_negative_absorption);
  const Propmat npm     = norm_view(pol, mag, los);
  const Propmat dnpm_du = dnorm_view_du(pol, mag, los);
  const Propmat dnpm_dv = dnorm_view_dv(pol, mag, los);
  const Propmat dnpm_dw = dnorm_view_dw(pol, mag, los);

  for (Size i = 0; i < pm.size(); i++) {
    spectral_propmat[i] += scale(npm, pm[i]);
    dispersion[i]       -= npm.A() * pm[i].imag();

    for (auto& atm_target : jac_targets.atm) {
      const Size j = atm_target.target_pos;

      std::visit(
          [&,
           &x = spectral_propmat_jac[j, i],
           &y = dispersion_jac[j, i]]<typename T>(const T& type) {
            if constexpr (std::same_as<T, AtmKey>) {
              if (type == AtmKey::mag_u) {
                x += scale(npm, dnpm_du, pm[i], dpm[j, i]);
                return;
              }
              if (type == AtmKey::mag_v) {
                x += scale(npm, dnpm_dv, pm[i], dpm[j, i]);
                return;
              }
              if (type == AtmKey::mag_w) {
                x += scale(npm, dnpm_dw, pm[i], dpm[j, i]);
                return;
              }
            }
            x += scale(npm, dpm[j, i]);
            y -= npm.A() * dpm[j, i].imag();
          },
          atm_target.type);
    }
  }
}
}  // namespace

void spectral_propmatAddVoigtLTE(PropmatVector& spectral_propmat,
                                 PropmatMatrix& spectral_propmat_jac,
                                 Vector& dispersion,
                                 Matrix& dispersion_jac,
                                 const AscendingGrid& freq_grid,
                                 const JacobianTargets& jac_targets,
                                 const SpeciesEnum& species,
                                 const AbsorptionBands& abs_bands,
                                 const AtmPoint& atm_point,
                                 const PropagationPathPoint& path_point,
                                 const Index& no_negative_absorption) try {
  ARTS_TIME_REPORT

  using enum ZeemanPolarization;

  //! FIXME: these should be part of workspace once things work?
  dispersion.resize(freq_grid.size());
  dispersion_jac.resize(jac_targets.target_count(), freq_grid.size());
  dispersion     = 0;
  dispersion_jac = 0;

  ARTS_USER_ERROR_IF(spectral_propmat.shape() != dispersion.shape() or
                         spectral_propmat.shape() != freq_grid.shape() or
                         spectral_propmat_jac.shape() != dispersion_jac.shape(),
                     R"(Inconsistent shapes:

spectral_propmat:     {:B,}
dispersion:           {:B,}
freq_grid:            {:B,}
spectral_propmat_jac: {:B,}
dispersion_jac:       {:B,}
)",
                     spectral_propmat.shape(),
                     dispersion.shape(),
                     freq_grid.shape(),
                     spectral_propmat_jac.shape(),
                     dispersion_jac.shape());

  const Size nf = freq_grid.size();
  const Size nt = spectral_propmat_jac.nrows();
  if (nf == 0) return;

  ComplexVector pm(nf, 0.0);
  ComplexMatrix dpm(nt, nf, 0.0);

  bool has_zeeman = lbl::voigt::lte::calculate(pm,
                                               dpm,
                                               abs_bands,
                                               freq_grid,
                                               jac_targets,
                                               atm_point,
                                               no,
                                               species,
                                               no_negative_absorption);

  spectral_propmat     += pm.real();
  spectral_propmat_jac += dpm.real();
  dispersion           -= pm.imag();
  dispersion_jac       -= dpm.imag();

  if (has_zeeman) {
    const auto mag = atm_point.mag;
    const auto los = path_point.los;
    for (auto pol : {sm, sp, pi}) {
      sumup_zeeman(spectral_propmat,
                   spectral_propmat_jac,
                   dispersion,
                   dispersion_jac,
                   pm,
                   dpm,
                   abs_bands,
                   freq_grid,
                   jac_targets,
                   atm_point,
                   species,
                   no_negative_absorption,
                   mag,
                   los,
                   pol);
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void single_propmatInit(Propmat& single_propmat,
                        PropmatVector& single_propmat_jac,
                        Stokvec& single_nlte_srcvec,
                        StokvecVector& single_nlte_srcvec_jac,
                        Numeric& single_dispersion,
                        Vector& single_dispersion_jac,
                        const JacobianTargets& jac_targets) try {
  ARTS_TIME_REPORT

  const Size nt = jac_targets.target_count();

  single_propmat = Propmat{};
  single_propmat_jac.resize(nt);
  single_propmat_jac = Propmat{};

  single_nlte_srcvec = Stokvec{};
  single_nlte_srcvec_jac.resize(nt);
  single_nlte_srcvec_jac = Stokvec{};

  single_dispersion = 0;
  single_dispersion_jac.resize(nt);
  single_dispersion_jac = 0;
}
ARTS_METHOD_ERROR_CATCH

void single_propmatAddVoigtLTE(Propmat& single_propmat,
                               PropmatVector& single_propmat_jac,
                               Numeric& single_dispersion,
                               Vector& single_dispersion_jac,
                               const Numeric& frequency,
                               const JacobianTargets& jac_targets,
                               const SpeciesEnum& species,
                               const AbsorptionBands& abs_bands,
                               const AtmPoint& atm_point,
                               const PropagationPathPoint& path_point,
                               const Index& no_negative_absorption) try {
  ARTS_TIME_REPORT

  using enum ZeemanPolarization;

  const Size nt = jac_targets.target_count();

  //! FIXME: these should be part of workspace once things work
  single_dispersion = 0;
  single_dispersion_jac.resize(nt);
  single_dispersion_jac = 0;

  ARTS_USER_ERROR_IF(
      single_propmat_jac.shape() != single_dispersion_jac.shape(),
      R"(Inconsistent shapes:

single_propmat_jac:    {:B,}
single_dispersion_jac: {:B,}
)",
      single_propmat_jac.shape(),
      single_dispersion_jac.shape());

  Complex pm_(0.0);
  ComplexVector dpm_(nt, 0.0);

  ComplexVectorView pm(pm_);
  ComplexMatrixView dpm = dpm_.view_as(nt, 1);
  const ConstVectorView freq_grid(frequency);
  PropmatVectorView spectral_propmat{single_propmat};
  PropmatMatrixView spectral_propmat_jac{single_propmat_jac.view_as(nt, 1)};
  VectorView dispersion{single_dispersion};
  MatrixView dispersion_jac{single_dispersion_jac.view_as(nt, 1)};

  bool has_zeeman = lbl::voigt::lte::calculate(pm,
                                               dpm,
                                               abs_bands,
                                               freq_grid,
                                               jac_targets,
                                               atm_point,
                                               no,
                                               species,
                                               no_negative_absorption);

  spectral_propmat     += pm.real();
  spectral_propmat_jac += dpm.real();
  dispersion           -= pm.imag();
  dispersion_jac       -= dpm.imag();

  if (has_zeeman) {
    const auto mag = atm_point.mag;
    const auto los = path_point.los;
    for (auto pol : {sm, sp, pi}) {
      sumup_zeeman(spectral_propmat,
                   spectral_propmat_jac,
                   dispersion,
                   dispersion_jac,
                   pm,
                   dpm,
                   abs_bands,
                   freq_grid,
                   jac_targets,
                   atm_point,
                   species,
                   no_negative_absorption,
                   mag,
                   los,
                   pol);
    }
  }
}
ARTS_METHOD_ERROR_CATCH
