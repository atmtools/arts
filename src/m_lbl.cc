#include <absorptionlines.h>
#include <enumsAbsorptionBandSortingOption.h>
#include <lbl.h>
#include <partfun.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <functional>
#include <iterator>
#include <limits>
#include <ranges>
#include <span>
#include <unordered_map>

#include "arts_omp.h"
#include "artstime.h"
#include "configtypes.h"
#include "debug.h"
#include "enumsLineByLineCutoffType.h"
#include "isotopologues.h"
#include "jacobian.h"
#include "lbl_data.h"
#include "lbl_lineshape_linemixing.h"
#include "lineshapemodel.h"
#include "matpack_view.h"
#include "minimize.h"
#include "path_point.h"
#include "quantum_numbers.h"
#include "rtepack.h"
#include "sorting.h"
#include "species.h"
#include "species_tags.h"
#include "xml_io.h"
#include "xml_io_base.h"

LineByLineLineshape toLineshape(const LineShapeTypeOld old_ls,
                                const AbsorptionPopulationTypeOld old_pop) try {
  if (old_ls == LineShapeTypeOld::VP and
      old_pop == AbsorptionPopulationTypeOld::LTE) {
    return LineByLineLineshape::VP_LTE;
  }

  ARTS_USER_ERROR(
      "New code does not support combination of {} and {}", old_ls, old_pop)
}
ARTS_METHOD_ERROR_CATCH

void absorption_bandsFromAbsorbtionLines(
    AbsorptionBands& absorption_bands,
    const ArrayOfArrayOfSpeciesTag& absorption_species,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species) try {
  ARTS_USER_ERROR_IF(absorption_species.size() != abs_lines_per_species.size(),
                     "{} != {}",
                     absorption_species.size(),
                     abs_lines_per_species.size())

  absorption_bands = {};
  absorption_bands.reserve(std::transform_reduce(
      abs_lines_per_species.begin(),
      abs_lines_per_species.end(),
      Size{0},
      std::plus<>{},
      [](const ArrayOfAbsorptionLines& lines) { return lines.size(); }));

  for (auto& abs_lines : abs_lines_per_species) {
    for (auto& old_band : abs_lines) {
      auto& new_band = absorption_bands[old_band.quantumidentity];

      new_band.lineshape =
          toLineshape(old_band.lineshapetype, old_band.population);
      new_band.cutoff = to<LineByLineCutoffType>(toString(old_band.cutoff));
      new_band.cutoff_value = old_band.cutofffreq;
      new_band.lines.reserve(old_band.lines.size());

      for (auto& old_line : old_band.lines) {
        auto& new_line = new_band.emplace_back();

        new_line.a      = old_line.A;
        new_line.f0     = old_line.F0;
        new_line.e0     = old_line.E0;
        new_line.gu     = old_line.gupp;
        new_line.gl     = old_line.glow;
        new_line.z.gu() = old_line.zeeman.gu();
        new_line.z.gl() = old_line.zeeman.gl();
        new_line.qn     = old_line.localquanta;

        new_line.ls.one_by_one = false;
        new_line.ls.T0         = old_band.T0;
        new_line.ls.single_models.reserve(old_band.broadeningspecies.size());

        for (Size i = 0; i < old_band.broadeningspecies.size(); i++) {
          auto& new_line_lsspec = new_line.ls.single_models.emplace_back();

          new_line_lsspec.species = old_band.broadeningspecies[i];
          for (std::string_view strvar :
               {"G0", "D0", "G2", "D2", "ETA", "FVC", "Y", "G", "DV"}) {
            auto old_value =
                old_line.lineshape[i].Get(to<LineShapeVariableOld>(strvar));
            const auto new_var = to<LineShapeModelVariable>(strvar);

            switch (old_value.type) {
              case LineShapeTemperatureModelOld::None: break;
              case LineShapeTemperatureModelOld::T0:
                new_line_lsspec.data.emplace_back(
                    new_var,
                    lbl::temperature::data{LineShapeModelType::T0,
                                           Vector{old_value.X0}});
                break;
              case LineShapeTemperatureModelOld::T1:
                new_line_lsspec.data.emplace_back(
                    new_var,
                    lbl::temperature::data{LineShapeModelType::T1,
                                           Vector{old_value.X0, old_value.X1}});
                break;
              case LineShapeTemperatureModelOld::T2:
                new_line_lsspec.data.emplace_back(
                    new_var,
                    lbl::temperature::data{
                        LineShapeModelType::T2,
                        Vector{old_value.X0, old_value.X1, old_value.X2}});
                break;
              case LineShapeTemperatureModelOld::T3:
                new_line_lsspec.data.emplace_back(
                    new_var,
                    lbl::temperature::data{LineShapeModelType::T3,
                                           Vector{old_value.X0, old_value.X1}});
                break;
              case LineShapeTemperatureModelOld::T4:
                new_line_lsspec.data.emplace_back(
                    new_var,
                    lbl::temperature::data{
                        LineShapeModelType::T4,
                        Vector{old_value.X0, old_value.X1, old_value.X2}});
                break;
              case LineShapeTemperatureModelOld::T5:
                new_line_lsspec.data.emplace_back(
                    new_var,
                    lbl::temperature::data{LineShapeModelType::T5,
                                           Vector{old_value.X0, old_value.X1}});
                break;
              case LineShapeTemperatureModelOld::LM_AER:
                new_line_lsspec.data.emplace_back(
                    new_var,
                    lbl::temperature::data{LineShapeModelType::AER,
                                           Vector{old_value.X0,
                                                  old_value.X1,
                                                  old_value.X2,
                                                  old_value.X3}});
                break;
              case LineShapeTemperatureModelOld::DPL:
                new_line_lsspec.data.emplace_back(
                    new_var,
                    lbl::temperature::data{LineShapeModelType::DPL,
                                           Vector{old_value.X0,
                                                  old_value.X1,
                                                  old_value.X2,
                                                  old_value.X3}});
                break;
              case LineShapeTemperatureModelOld::POLY:
                new_line_lsspec.data.emplace_back(
                    new_var,
                    lbl::temperature::data{LineShapeModelType::POLY,
                                           Vector{old_value.X0,
                                                  old_value.X1,
                                                  old_value.X2,
                                                  old_value.X3}});
                break;
            }
          }
        }
      }
    }
  }

  for (auto& [key, bnd] : absorption_bands) {
    bnd.sort();
    for (auto& line : bnd.lines) {
      line.ls.clear_zeroes();
    }
  }
}
ARTS_METHOD_ERROR_CATCH

std::pair<LineShapeTypeOld, AbsorptionPopulationTypeOld>
toLineshapeAndPolpulation(LineByLineLineshape x) try {
  if (x == LineByLineLineshape::VP_LTE) {
    return {LineShapeTypeOld::VP, AbsorptionPopulationTypeOld::LTE};
  }

  ARTS_USER_ERROR("Old code does not support conversion from {}", x)
}
ARTS_METHOD_ERROR_CATCH

void abs_linesFromArrayOfAbsorptionBands(
    ArrayOfAbsorptionLines& abs_lines,
    const AbsorptionBands& absorption_bands) try {
  abs_lines.resize(0);
  abs_lines.reserve(absorption_bands.size());

  for (auto& [key, band] : absorption_bands) {
    if (band.lines.size() == 0) continue;

    const Size i0 = abs_lines.size();
    AbsorptionLines old_band;

    old_band.quantumidentity = key;
    old_band.cutoff        = to<AbsorptionCutoffTypeOld>(toString(band.cutoff));
    old_band.cutofffreq    = band.cutoff_value;
    const auto [ls, pop]   = toLineshapeAndPolpulation(band.lineshape);
    old_band.lineshapetype = ls;
    old_band.population    = pop;
    old_band.normalization = AbsorptionNormalizationTypeOld::SFS;
    old_band.linemixinglimit = -1;
    old_band.lines.resize(1);

    for (auto& line : band.lines) {
      AbsorptionSingleLine& old_line = old_band.lines.front();

      old_band.broadeningspecies.resize(line.ls.single_models.size());
      for (Size i = 0; i < line.ls.single_models.size(); i++) {
        old_band.broadeningspecies[i] = line.ls.single_models[i].species;
      }
      old_band.bathbroadening =
          old_band.broadeningspecies.back() == SpeciesEnum::Bath;

      old_line.A           = line.a;
      old_line.F0          = line.f0;
      old_line.E0          = line.e0;
      old_line.gupp        = line.gu;
      old_line.glow        = line.gl;
      old_line.zeeman.gu() = line.z.gu();
      old_line.zeeman.gl() = line.z.gl();
      old_line.localquanta = line.qn;
      old_line.I0 =
          -(Math::pow2(Constant::c) / (8 * Constant::pi)) * line.f0 *
          std::expm1(-Constant::h * line.f0 / (Constant::k * line.ls.T0)) *
          line.s(line.ls.T0,
                 PartitionFunctions::Q(line.ls.T0, key.Isotopologue()));

      old_line.lineshape.resize(line.ls.single_models.size());
      for (Size i = 0; i < line.ls.single_models.size(); i++) {
        for (auto& [var, data] : line.ls.single_models[i].data) {
          auto& old_line_ls = old_line.lineshape[i];
          ARTS_USER_ERROR_IF(
              data.X().size() > 4,
              "Old code does not support more than 4 temperature parameters.")
          LineShape::ModelParameters old_model(
              to<LineShapeTemperatureModelOld>(toString(data.Type())),
              data.X());
          old_line_ls.Set(to<LineShapeVariableOld>(toString(var)), old_model);
        }
      }

      if (auto ptr = std::find_if(abs_lines.begin() + i0,
                                  abs_lines.end(),
                                  [&old_band](const AbsorptionLines& l) {
                                    return old_band.Match(l).first;
                                  });
          ptr == abs_lines.end()) {
        abs_lines.push_back(old_band);
      } else {
        ptr->AppendSingleLine(old_line);
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

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

void absorption_bandsSelectFrequency(AbsorptionBands& absorption_bands,
                                     const Numeric& fmin,
                                     const Numeric& fmax,
                                     const Index& by_line) try {
  std::vector<QuantumIdentifier> to_remove;

  for (auto& [key, band] : absorption_bands) {
    if (band.lines.front().f0 > fmax or band.lines.back().f0 < fmin) {
      to_remove.push_back(key);
    }
  }

  for (const auto& key : to_remove) {
    absorption_bands.erase(key);
  }

  if (by_line) {
    for (auto& [key, band] : absorption_bands) {
      auto& lines = band.lines;
      lines.erase(std::remove_if(lines.begin(),
                                 lines.end(),
                                 [fmin, fmax](const lbl::line& l) {
                                   return l.f0 < fmin or l.f0 > fmax;
                                 }),
                  lines.end());
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void absorption_bandsRemoveID(AbsorptionBands& absorption_bands,
                              const QuantumIdentifier& id) try {
  absorption_bands.erase(id);
}
ARTS_METHOD_ERROR_CATCH

void sortedIndexOfBands(ArrayOfIndex& sorted_idxs,
                        const AbsorptionBands& absorption_bands,
                        const String& criteria,
                        const Index& reverse,
                        const Numeric& temperature) try {
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
            [T = temperature, ir = key.Isotopologue()](const lbl::line& l) {
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
    auto [first, last] = std::ranges::equal_range(
        qid_sorter | std::views::drop(i),
        qid_sorter[i],
        [](const auto& p1, const auto& p2) {
          return p1.qid.isotopologue_index < p2.qid.isotopologue_index;
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
  const bool ignore_missing = static_cast<bool>(ignore_missing_);
  absorption_bands.clear();

  const String my_base = complete_basename(basename);

  std::set<SpeciesIsotope> isotopologues;
  for (auto& specs : absorbtion_species) {
    for (auto& spec : specs) {
      if (spec.type == SpeciesTagType::Plain) {
        if (spec.is_joker()) {
          for (auto&& isot : Species::isotopologues(spec.Spec())) {
            if (is_predefined_model(isot)) continue;
            if (isot.joker()) continue;
            isotopologues.insert(isot);
          }
        } else {
          isotopologues.insert(spec.Isotopologue());
        }
      }
    }
  }

  for (auto& isot : isotopologues) {
    String filename = my_base + isot.FullName() + ".xml";
    if (find_xml_file_existence(filename)) {
      AbsorptionBands other;
      xml_read_from_file(filename, other);
      absorption_bands.insert(std::make_move_iterator(other.begin()),
                              std::make_move_iterator(other.end()));
    } else {
      ARTS_USER_ERROR_IF(not ignore_missing, "File {} not found", filename)
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void absorption_bandsReadSplit(AbsorptionBands& absorption_bands,
                               const String& dir) try {
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

  ARTS_USER_ERROR_IF(not error.empty(), "{}", error)

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
  auto create_if_not = [](const std::filesystem::path& path) {
    if (not std::filesystem::exists(path)) {
      std::filesystem::create_directories(path);
    }
    return path;
  };

  const auto p = create_if_not(dir);

  std::unordered_map<SpeciesIsotope, AbsorptionBands> isotopologues_data;
  for (auto& [key, band] : absorption_bands) {
    isotopologues_data[key.Isotopologue()][key] = band;
  }

  for (const auto& [isot, bands] : isotopologues_data) {
    xml_write_to_file(
        (p / var_string(isot, ".xml")).string(), bands, FileType::ascii, 0);
  }
}
ARTS_METHOD_ERROR_CATCH

void absorption_bandsSetZeeman(AbsorptionBands& absorption_bands,
                               const SpeciesIsotope& species,
                               const Numeric& fmin,
                               const Numeric& fmax,
                               const Index& _on) try {
  const bool on = static_cast<bool>(_on);

  for (auto& [key, band] : absorption_bands) {
    if (key.Isotopologue() != species) continue;

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
  const auto n = arts_omp_get_max_threads();
  if (n == 1 or arts_omp_in_parallel() or n > f_grid.size()) {
    lbl::calculate(pm,
                   sv,
                   dpm,
                   dsv,
                   f_grid,
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
    for (Index i = 0; i < n; i++) {
      try {
        lbl::calculate(pm.slice(ompv[i].first, ompv[i].second),
                       sv.slice(ompv[i].first, ompv[i].second),
                       dpm(joker, Range(ompv[i].first, ompv[i].second)),
                       dsv(joker, Range(ompv[i].first, ompv[i].second)),
                       static_cast<const Vector&>(f_grid).slice(ompv[i].first,
                                                                ompv[i].second),
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

    ARTS_USER_ERROR_IF(not error.empty(), "{}", error)
  }
}
ARTS_METHOD_ERROR_CATCH

void absorption_bandsReadHITRAN(AbsorptionBands& absorption_bands,
                                const String& filename,
                                const Vector2& frequency_range,
                                const String& line_strength_option,
                                const Index& compute_zeeman_parameters) try {
  using namespace Quantum::Number;

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

void absorption_bandsLineMixingAdaptation(
    AbsorptionBands& absorption_bands,
    const LinemixingEcsData& ecs_data,
    const AtmPoint& atmospheric_point,
    const AscendingGrid& temperatures,
    const QuantumIdentifier& band_key,
    const Index& rosenkranz_fit_order,
    const Index& polynomial_fit_degree) try {
  ARTS_USER_ERROR_IF(temperatures.empty() or temperatures.front() <= 0.0,
                     "Need a positive temperature grid")

  ARTS_USER_ERROR_IF(rosenkranz_fit_order != 1 and rosenkranz_fit_order != 2,
                     "Only 1 or 2 is supported for the ordered fit")
  ARTS_USER_ERROR_IF(
      polynomial_fit_degree < 1 or
          polynomial_fit_degree > temperatures.size() - 1,
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
                                     ecs_data.data.at(band_key.Isotopologue()),
                                     atmospheric_point,
                                     temperatures);

  band.sort(LineByLineVariable::f0);

  Matrix lbl_str(M, N);
  for (Size i = 0; i < M; i++) {
    const Numeric Q =
        PartitionFunctions::Q(temperatures[i], band_key.Isotopologue());
    for (Size k = 0; k < N; k++) {
      auto& line    = band.lines[k];
      lbl_str(i, k) = line.s(temperatures[i], Q) * Math::pow2(Constant::c) /
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
        lbl_val(i, j, k) = Complex{
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
      auto s = eqv_str(i, j, joker);
      auto v = eqv_val(i, j, joker);
      bubble_sort_by(
          [&](auto I1, auto I2) { return v[I1].real() > v[I2].real(); }, s, v);
    }
  }

  eqv_val -= lbl_val;

  for (Size j = 0; j < K; j++) eqv_str(joker, j, joker) /= lbl_str;

  eqv_val.real() /= Math::pow2(atmospheric_point.pressure);
  eqv_val.imag() /= Math::pow3(atmospheric_point.pressure);

  eqv_str.real() -= 1.0;
  eqv_str.real() /= Math::pow2(atmospheric_point.pressure);
  eqv_str.imag() /= atmospheric_point.pressure;

  using namespace Minimize;
  for (Size j = 0; j < K; j++) {
    for (Size k = 0; k < N; k++) {
      if (rosenkranz_fit_order >= 1) {
        auto [success_y, yfit] = curve_fit<Polynom>(
            temperatures, eqv_str(joker, j, k).imag(), polynomial_fit_degree);
        ARTS_USER_ERROR_IF(
            not success_y, "Cannot fit y for line {} of band {}", k, band_key)
        band.lines[k].ls.single_models[j].data.emplace_back(
            LineShapeModelVariable::Y,
            lbl::temperature::data{LineShapeModelType::POLY, Vector{yfit}});
      }

      if (rosenkranz_fit_order >= 2) {
        auto [success_g, gfit] = curve_fit<Polynom>(
            temperatures, eqv_str(joker, j, k).real(), polynomial_fit_degree);
        ARTS_USER_ERROR_IF(
            not success_g, "Cannot fit g for line {} of band {}", k, band_key)
        band.lines[k].ls.single_models[j].data.emplace_back(
            LineShapeModelVariable::G,
            lbl::temperature::data{LineShapeModelType::POLY, Vector{gfit}});

        auto [success_d, dfit] = curve_fit<Polynom>(
            temperatures, eqv_val(joker, j, k).real(), polynomial_fit_degree);
        ARTS_USER_ERROR_IF(
            not success_d, "Cannot fit dv for line {} of band {}", k, band_key)
        band.lines[k].ls.single_models[j].data.emplace_back(
            LineShapeModelVariable::DV,
            lbl::temperature::data{LineShapeModelType::POLY, Vector{dfit}});
      }
    }
  }

  for (auto& line : band.lines) line.ls.one_by_one = orig_one_by_one;
}
ARTS_METHOD_ERROR_CATCH
