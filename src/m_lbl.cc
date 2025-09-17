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
#include <functional>
#include <iterator>
#include <ranges>
#include <span>
#include <unordered_map>

namespace {
std::vector<std::pair<Index, Index>> omp_offset_count(const Index N,
                                                      const Index n) {
  std::vector<std::pair<Index, Index>> result(n, {0, 0});
  const Index dn        = N / n;
  result.front().second = dn;

  for (Index i = 1; i < n - 1; i++) {
    result[i].first  = result[i - 1].first + dn;
    result[i].second = dn;
  }

  result.back().first  = result[n - 2].first + dn;
  result.back().second = N - result.back().first;

  return result;
}
}  // namespace

void absorption_bandsSelectFrequencyByLine(AbsorptionBands& absorption_bands,
                                           const Numeric& fmin,
                                           const Numeric& fmax) try {
  ARTS_TIME_REPORT

  std::vector<QuantumIdentifier> to_remove;

  for (auto& [key, band] : absorption_bands) {
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

  for (const auto& key : to_remove) absorption_bands.erase(key);
}
ARTS_METHOD_ERROR_CATCH

void absorption_bandsSelectFrequencyByBand(AbsorptionBands& absorption_bands,
                                           const Numeric& fmin,
                                           const Numeric& fmax) try {
  ARTS_TIME_REPORT

  std::vector<QuantumIdentifier> to_remove;

  for (auto& [key, band] : absorption_bands) {
    band.sort();

    const bool criteria = band.lines.empty() or band.lines.front().f0 > fmax or
                          band.lines.back().f0 < fmin;

    if (criteria) to_remove.push_back(key);
  }

  for (const auto& key : to_remove) absorption_bands.erase(key);
}
ARTS_METHOD_ERROR_CATCH

void absorption_bandsRemoveID(AbsorptionBands& absorption_bands,
                              const QuantumIdentifier& id) try {
  ARTS_TIME_REPORT

  absorption_bands.erase(id);
}
ARTS_METHOD_ERROR_CATCH

void sortedIndexOfBands(ArrayOfIndex& sorted_idxs,
                        const AbsorptionBands& absorption_bands,
                        const String& criteria,
                        const Index& reverse,
                        const Numeric& temperature) try {
  ARTS_TIME_REPORT

  struct order {
    QuantumIdentifier qid;
    Numeric value;
    Index idx;
  };

  std::vector<order> qid_sorter;
  qid_sorter.reserve(absorption_bands.size());

  const auto sort_opt = to<AbsorptionBandSortingOption>(criteria);

  for (auto& [key, band] : absorption_bands) {
    auto& v = qid_sorter.emplace_back(key, 0.0, qid_sorter.size()).value;

    switch (sort_opt) {
      case AbsorptionBandSortingOption::IntegratedIntensity:
        v = std::transform_reduce(
            band.lines.begin(),
            band.lines.end(),
            Numeric{0},
            std::plus<>{},
            [T = temperature](const lbl::line& l) {
              return -l.f0 *
                     std::expm1(-Constant::h * l.f0 / (Constant::k * T)) *
                     l.s(T, 1);
            });
        break;
      case AbsorptionBandSortingOption::FrontFrequency:
        if (band.size()) v = band.lines.front().f0;
        break;
      case AbsorptionBandSortingOption::None: break;
    }
  }

  std::ranges::sort(qid_sorter, {}, &order::qid);

  Size i = 0;
  while (i < qid_sorter.size()) {
    auto [first, last] =
        std::ranges::equal_range(qid_sorter | std::views::drop(i),
                                 qid_sorter[i],
                                 [](const auto& p1, const auto& p2) {
                                   return p1.qid.isot < p2.qid.isot;
                                 });

    auto span  = std::span{first, last};
    i         += span.size();

    if (reverse) {
      std::ranges::sort(span | std::views::reverse, {}, &order::value);
    } else {
      std::ranges::sort(span, {}, &order::value);
    }
  }

  sorted_idxs.resize(0);
  sorted_idxs.reserve(qid_sorter.size());
  std::ranges::move(qid_sorter | std::views::transform(&order::idx),
                    std::back_inserter(sorted_idxs));
}
ARTS_METHOD_ERROR_CATCH

void absorption_bandsKeepID(AbsorptionBands& absorption_bands,
                            const QuantumIdentifier& id,
                            const Index& line) try {
  ARTS_TIME_REPORT

  if (auto ptr = absorption_bands.find(id); ptr != absorption_bands.end()) {
    const auto& [key, band] = *ptr;

    QuantumIdentifier newk = key;
    AbsorptionBand newb    = band;
    absorption_bands       = {};
    AbsorptionBand& data   = absorption_bands[newk];
    data                   = std::move(newb);

    if (line >= 0) {
      ARTS_USER_ERROR_IF(static_cast<Size>(line) >= band.lines.size(),
                         "Line index out of range: {}",
                         line)
      data.lines = {data.lines[line]};
    }
  } else {
    absorption_bands = {};
  }
}
ARTS_METHOD_ERROR_CATCH

void absorption_bandsReadSpeciesSplitCatalog(
    AbsorptionBands& absorption_bands,
    const ArrayOfArrayOfSpeciesTag& absorbtion_species,
    const String& basename,
    const Index& ignore_missing_) try {
  ARTS_TIME_REPORT

  absorption_bands.clear();

  const bool ignore_missing = static_cast<bool>(ignore_missing_);

  const String my_base = complete_basename(basename);

  std::set<SpeciesIsotope> isotopologues;
  for (auto& specs : absorbtion_species) {
    for (auto& spec : specs) {
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
#pragma omp critical(absorption_bandsReadSpeciesSplitCatalogInsert)
        absorption_bands.insert(std::make_move_iterator(other.begin()),
                                std::make_move_iterator(other.end()));
      } else if (not ignore_missing) {
#pragma omp critical(absorption_bandsReadSpeciesSplitCatalogFileError)
        file_errors.push_back(filename);
      }
    } catch (const std::exception& e) {
#pragma omp critical(absorption_bandsReadSpeciesSplitCatalogReadError)
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

void absorption_bandsReadSplit(AbsorptionBands& absorption_bands,
                               const String& dir) try {
  ARTS_TIME_REPORT

  absorption_bands = {};

  std::vector<std::filesystem::path> paths;
  std::ranges::copy_if(
      std::filesystem::directory_iterator(std::filesystem::path(dir)),
      std::back_inserter(paths),
      [](auto& entry) {
        return entry.is_regular_file() and entry.path().extension() == ".xml";
      });
  std::ranges::sort(paths);

  std::vector<AbsorptionBands> splitbands(paths.size());
  std::string error{};

#pragma omp parallel for schedule(dynamic)
  for (Size i = 0; i < paths.size(); i++) {
    try {
      xml_read_from_file(paths[i].string(), splitbands[i]);
    } catch (std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  if (not error.empty()) throw std::runtime_error(error);

  absorption_bands.reserve(std::transform_reduce(
      splitbands.begin(),
      splitbands.end(),
      Size{0},
      std::plus<>{},
      [](const AbsorptionBands& bands) { return bands.size(); }));
  for (auto& bands : splitbands) {
    for (auto& [key, data] : bands) {
      ARTS_USER_ERROR_IF(absorption_bands.find(key) != absorption_bands.end(),
                         "Read multiple bands of ID: {}",
                         key)
      absorption_bands[key] = std::move(data);
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void absorption_bandsSaveSplit(const AbsorptionBands& absorption_bands,
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
  for (auto& [key, band] : absorption_bands) {
    isotopologues_data[key.isot][key] = band;
  }

  for (const auto& [isot, bands] : isotopologues_data) {
    xml_write_to_file(
        (p / std::format("{}.xml", isot)).string(), bands, FileType::ascii, 0);
  }
}
ARTS_METHOD_ERROR_CATCH

void absorption_bandsSetZeeman(AbsorptionBands& absorption_bands,
                               const SpeciesIsotope& species,
                               const Numeric& fmin,
                               const Numeric& fmax,
                               const Index& _on) try {
  ARTS_TIME_REPORT

  const bool on = static_cast<bool>(_on);

  for (auto& [key, band] : absorption_bands) {
    if (key.isot != species) continue;

    for (auto& line : band.lines) {
      if (line.f0 >= fmin and line.f0 <= fmax) {
        line.z.on = on;
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void propagation_matrixAddLines(PropmatVector& pm,
                                StokvecVector& sv,
                                PropmatMatrix& dpm,
                                StokvecMatrix& dsv,
                                const AscendingGrid& f_grid,
                                const JacobianTargets& jacobian_targets,
                                const SpeciesEnum& species,
                                const AbsorptionBands& absorption_bands,
                                const LinemixingEcsData& ecs_data,
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
                   jacobian_targets,
                   species,
                   absorption_bands,
                   ecs_data,
                   atm_point,
                   path_point.los,
                   no_negative_absorption);
  } else {
    const auto ompv = omp_offset_count(f_grid.size(), n);
    std::string error;
#pragma omp parallel for
    for (Size i = 0; i < n; i++) {
      try {
        const Index offset = ompv[i].first;
        const Index nelem  = ompv[i].second;
        const Range f_range(offset, nelem);
        lbl::calculate(pm,
                       sv,
                       dpm,
                       dsv,
                       f_grid,
                       f_range,
                       jacobian_targets,
                       species,
                       absorption_bands,
                       ecs_data,
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

void absorption_bandsReadHITRAN(AbsorptionBands& absorption_bands,
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

  absorption_bands = {};
  for (auto& line : data) {
    auto [mapped_band, _] = absorption_bands.try_emplace(
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

void absorption_bandsLineMixingAdaptation(
    AbsorptionBands& absorption_bands,
    const LinemixingEcsData& ecs_data,
    const AtmPoint& atmospheric_point,
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

  auto& band = get_value(band_key, absorption_bands);

  if (band.lines.empty()) return;

  const auto& orig_front_ls  = band.lines.front().ls;
  const bool orig_one_by_one = orig_front_ls.one_by_one;

  for (auto& line : band.lines | std::ranges::views::drop(1)) {
    ARTS_USER_ERROR_IF(
        orig_one_by_one != line.ls.one_by_one,
        "Inconsistent line shape models - all lines must be consistently set to use one-by-one or not")

    ARTS_USER_ERROR_IF(
        not std::ranges::equal(line.ls.single_models,
                               orig_front_ls.single_models,
                               [](const auto& lsl, const auto& lsr) {
                                 return lsl.species == lsr.species;
                               }),
        "Inconsistent line shape models, all lines must have the same broadening species")

    line.ls.one_by_one = true;
  }
  band.lines.front().ls.one_by_one = true;

  const Size K = band.lines.front().ls.single_models.size();
  const Size M = temperatures.size();
  const Size N = band.lines.size();

  lbl::voigt::ecs::ComputeData com_data({}, atmospheric_point);
  ComplexTensor3 eqv_str(M, K, N), eqv_val(M, K, N);

  lbl::voigt::ecs::equivalent_values(eqv_str,
                                     eqv_val,
                                     com_data,
                                     band_key,
                                     band,
                                     ecs_data.at(band_key.isot),
                                     atmospheric_point,
                                     temperatures);

  band.sort(LineByLineVariable::f0);

  Matrix lbl_str(M, N);
  for (Size i = 0; i < M; i++) {
    const Numeric Q = PartitionFunctions::Q(temperatures[i], band_key.isot);
    for (Size k = 0; k < N; k++) {
      auto& line    = band.lines[k];
      lbl_str[i, k] = line.s(temperatures[i], Q) * Math::pow2(Constant::c) /
                      (8 * Constant::pi);
    }
  }

  for (auto& line : band.lines) {
    for (auto& lsm : line.ls.single_models) {
      lsm.remove_variables<LineShapeModelVariable::Y,
                           LineShapeModelVariable::G,
                           LineShapeModelVariable::DV>();
    }
  }

  ComplexTensor3 lbl_val(M, K, N);
  for (Size i = 0; i < M; i++) {
    for (Size j = 0; j < K; j++) {
      for (Size k = 0; k < N; k++) {
        auto& line       = band.lines[k];
        lbl_val[i, j, k] = Complex{
            line.f0 + line.ls.single_models[j].D0(line.ls.T0,
                                                  temperatures[i],
                                                  atmospheric_point.pressure),
            line.ls.single_models[j].G0(
                line.ls.T0, temperatures[i], atmospheric_point.pressure)};
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

  eqv_val.real() /= Math::pow2(atmospheric_point.pressure);
  eqv_val.imag() /= Math::pow3(atmospheric_point.pressure);

  eqv_str.real() -= 1.0;
  eqv_str.real() /= Math::pow2(atmospheric_point.pressure);
  eqv_str.imag() /= atmospheric_point.pressure;

  using namespace Minimize;
  for (Size j = 0; j < K; j++) {
    for (Size k = 0; k < N; k++) {
      if (rosenkranz_fit_order >= 1) {
        auto yfit = polyfit(
            temperatures, eqv_str[joker, j, k].imag(), polynomial_fit_degree);
        ARTS_USER_ERROR_IF(
            not yfit, "Cannot fit y for line {} of band {}", k, band_key)
        band.lines[k].ls.single_models[j].data[LineShapeModelVariable::Y] =
            lbl::temperature::data{LineShapeModelType::POLY, *yfit};
      }

      if (rosenkranz_fit_order >= 2) {
        auto gfit = polyfit(
            temperatures, eqv_str[joker, j, k].real(), polynomial_fit_degree);
        ARTS_USER_ERROR_IF(
            not gfit, "Cannot fit g for line {} of band {}", k, band_key)
        band.lines[k].ls.single_models[j].data[LineShapeModelVariable::G] =
            lbl::temperature::data{LineShapeModelType::POLY, *gfit};

        auto dfit = polyfit(
            temperatures, eqv_val[joker, j, k].real(), polynomial_fit_degree);
        ARTS_USER_ERROR_IF(
            not dfit, "Cannot fit dv for line {} of band {}", k, band_key)
        band.lines[k].ls.single_models[j].data[LineShapeModelVariable::DV] =
            lbl::temperature::data{LineShapeModelType::POLY, *dfit};
      }
    }
  }

  for (auto& line : band.lines) line.ls.one_by_one = orig_one_by_one;
}
ARTS_METHOD_ERROR_CATCH

void absorption_bandsReadSpeciesSplitARTSCAT(
    AbsorptionBands& absorption_bands,
    const ArrayOfArrayOfSpeciesTag& absorbtion_species,
    const String& basename,
    const Index& ignore_missing_,
    const Index& pure_species_) try {
  ARTS_TIME_REPORT

  Array<ArrayOfArtscatMeta> meta;

  const bool ignore_missing = static_cast<bool>(ignore_missing_);
  const bool pure_species   = static_cast<bool>(pure_species_);

  absorption_bands.clear();

  const String my_base = complete_basename(basename);

  std::vector<std::string> file_errors;

  if (pure_species) {
    std::set<SpeciesEnum> specieses;
    for (auto& specs : absorbtion_species) {
      for (auto& spec : specs) {
        if (spec.type == SpeciesTagType::Plain) {
          specieses.insert(spec.Spec());
        }
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
    for (auto& specs : absorbtion_species) {
      for (auto& spec : specs) {
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
      absorption_bands[std::move(m.quantumidentity)].emplace_back(
          std::move(m.data));
    }
  }
}
ARTS_METHOD_ERROR_CATCH
