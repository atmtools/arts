#include <absorptionlines.h>
#include <lbl.h>
#include <partfun.h>

#include <exception>
#include <stdexcept>

#include "arts_omp.h"
#include "debug.h"
#include "lbl_data.h"
#include "new_jacobian.h"
#include "quantum_numbers.h"

lbl::Lineshape toLineshape(const LineShape::Type old_ls,
                           const Absorption::PopulationType old_pop) {
  if (old_ls == LineShape::Type::VP and
      old_pop == Absorption::PopulationType::LTE) {
    return lbl::Lineshape::VP_LTE;
  }

  ARTS_USER_ERROR(
      "New code does not support combination of ", old_ls, " and ", old_pop)
}

void absorption_bandsFromAbsorbtionLines(
    AbsorptionBands& absorption_bands,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const Numeric& allowed_linestrength_error) {
  absorption_bands.resize(0);
  absorption_bands.reserve(std::transform_reduce(
      abs_lines_per_species.begin(),
      abs_lines_per_species.end(),
      Size{0},
      std::plus<>{},
      [](const ArrayOfAbsorptionLines& lines) { return lines.size(); }));

  std::string errors;
  for (auto& abs_lines : abs_lines_per_species) {
    for (auto& old_band : abs_lines) {
      auto& [new_key, new_band] = absorption_bands.emplace_back();
      new_key = old_band.quantumidentity;
      new_band.lineshape =
          toLineshape(old_band.lineshapetype, old_band.population);
      new_band.cutoff = lbl::toCutoffTypeOrThrow(toString(old_band.cutoff));
      new_band.cutoff_value = old_band.cutofffreq;
      new_band.lines.reserve(old_band.lines.size());

      for (auto& old_line : old_band.lines) {
        auto& new_line = new_band.emplace_back();

        new_line.a = old_line.A;
        new_line.f0 = old_line.F0;
        new_line.e0 = old_line.E0;
        new_line.gu = old_line.gupp;
        new_line.gl = old_line.glow;
        new_line.z.gu() = old_line.zeeman.gu();
        new_line.z.gl() = old_line.zeeman.gl();
        new_line.qn = old_line.localquanta;

        new_line.ls.one_by_one = false;
        new_line.ls.T0 = old_band.T0;
        new_line.ls.single_models.reserve(old_band.broadeningspecies.size());

        if (allowed_linestrength_error > 0.0) {
          const Numeric new_i0 =
              -new_line.f0 *
              std::expm1(-Constant::h * new_line.f0 /
                         (Constant::k * new_line.ls.T0)) *
              new_line.s(
                  old_band.T0,
                  PartitionFunctions::Q(old_band.T0, new_key.Isotopologue()));
          const Numeric percentage_diff =
              100 * std::abs(1.0 - old_line.I0 / new_i0);
          if (percentage_diff > allowed_linestrength_error)
            errors += var_string(new_key,
                                 " diff: ",
                                 percentage_diff,
                                 "%; line: ",
                                 old_line,
                                 '\n');
        }

        for (Size i = 0; i < old_band.broadeningspecies.size(); i++) {
          auto& new_line_lsspec = new_line.ls.single_models.emplace_back();

          new_line_lsspec.species = old_band.broadeningspecies[i];
          for (std::string_view strvar :
               {"G0", "D0", "G2", "D2", "ETA", "FVC", "Y", "G", "DV"}) {
            auto old_value =
                old_line.lineshape[i].Get(LineShape::toVariableOrThrow(strvar));
            const auto new_var = lbl::line_shape::tovariableOrThrow(strvar);

            switch (old_value.type) {
              case LineShapeTemperatureModel::None:
                break;
              case LineShapeTemperatureModel::T0:
                new_line_lsspec.data.emplace_back(
                    new_var,
                    lbl::temperature::data{lbl::temperature::model_type::T0,
                                           Vector{old_value.X0}});
                break;
              case LineShapeTemperatureModel::T1:
                new_line_lsspec.data.emplace_back(
                    new_var,
                    lbl::temperature::data{lbl::temperature::model_type::T1,
                                           Vector{old_value.X0, old_value.X1}});
                break;
              case LineShapeTemperatureModel::T2:
                new_line_lsspec.data.emplace_back(
                    new_var,
                    lbl::temperature::data{
                        lbl::temperature::model_type::T2,
                        Vector{old_value.X0, old_value.X1, old_value.X2}});
                break;
              case LineShapeTemperatureModel::T3:
                new_line_lsspec.data.emplace_back(
                    new_var,
                    lbl::temperature::data{lbl::temperature::model_type::T3,
                                           Vector{old_value.X0, old_value.X1}});
                break;
              case LineShapeTemperatureModel::T4:
                new_line_lsspec.data.emplace_back(
                    new_var,
                    lbl::temperature::data{
                        lbl::temperature::model_type::T4,
                        Vector{old_value.X0, old_value.X1, old_value.X2}});
                break;
              case LineShapeTemperatureModel::T5:
                new_line_lsspec.data.emplace_back(
                    new_var,
                    lbl::temperature::data{lbl::temperature::model_type::T5,
                                           Vector{old_value.X0, old_value.X1}});
                break;
              case LineShapeTemperatureModel::LM_AER:
                new_line_lsspec.data.emplace_back(
                    new_var,
                    lbl::temperature::data{lbl::temperature::model_type::AER,
                                           Vector{old_value.X0,
                                                  old_value.X1,
                                                  old_value.X2,
                                                  old_value.X3}});
                break;
              case LineShapeTemperatureModel::DPL:
                new_line_lsspec.data.emplace_back(
                    new_var,
                    lbl::temperature::data{lbl::temperature::model_type::DPL,
                                           Vector{old_value.X0,
                                                  old_value.X1,
                                                  old_value.X2,
                                                  old_value.X3}});
                break;
              case LineShapeTemperatureModel::POLY:
                new_line_lsspec.data.emplace_back(
                    new_var,
                    lbl::temperature::data{lbl::temperature::model_type::POLY,
                                           Vector{old_value.X0,
                                                  old_value.X1,
                                                  old_value.X2,
                                                  old_value.X3}});
                break;
              case LineShapeTemperatureModel::FINAL:;
            }
          }
        }
      }
    }
  }

  if (errors.size()) {
    throw std::runtime_error(var_string(errors));
  }

  std::ranges::sort(
      absorption_bands,
      [](const lbl::band_data& lhs, const lbl::band_data& rhs) {
        return lhs.size() > rhs.size();
      },
      &lbl::band::data);
  for (auto& bnd : absorption_bands) bnd.data.sort();
}

std::vector<std::pair<Index, Index>> omp_offset_count(const Index N,
                                                      const Index n) {
  std::vector<std::pair<Index, Index>> result(n, {0, 0});
  const Index dn = N / n;
  result.front().second = dn;

  for (Index i = 1; i < n - 1; i++) {
    result[i].first = result[i - 1].first + dn;
    result[i].second = dn;
  }

  result.back().first = result[n - 2].first + dn;
  result.back().second = N - result.back().first;

  return result;
}

void absorption_bandsSelectFrequency(AbsorptionBands& absorption_bands,
                                     const Numeric& fmin,
                                     const Numeric& fmax) {
  std::vector<Size> to_remove;

  for (Size i = 0; i < absorption_bands.size(); i++) {
    if (absorption_bands[i].data.lines.front().f0 > fmax or
        absorption_bands[i].data.lines.back().f0 < fmin) {
      to_remove.push_back(i);
    }
  }

  for (auto i : to_remove | std::views::reverse) {
    absorption_bands.erase(absorption_bands.begin() + i);
  }
}

void absorption_bandsRemoveID(AbsorptionBands& absorption_bands,
                              const QuantumIdentifier& id) {
  for (Size i = 0; i < absorption_bands.size(); i++) {
    if (id == absorption_bands[i].key) {
      absorption_bands.erase(absorption_bands.begin() + i);
      return;
    }
  }
  ARTS_USER_ERROR("Did not find band of ID: ", id)
}

void absorption_bandsKeepID(AbsorptionBands& absorption_bands,
                            const QuantumIdentifier& id,
                            const Index& line) {
  for (Size i = 0; i < absorption_bands.size(); i++) {
    if (id == absorption_bands[i].key) {
      absorption_bands = {absorption_bands[i]};

      if (line >= 0) {
        ARTS_USER_ERROR_IF(
            static_cast<Size>(line) >= absorption_bands[0].data.lines.size(),
            "Line index out of range: ",
            line)
        absorption_bands[0].data.lines = {absorption_bands[0].data.lines[line]};
      }

      return;
    }
  }
}

void propmat_clearskyAddLines2(PropmatVector& pm,
                               PropmatMatrix& dpm,
                               const Vector& f_grid,
                               const JacobianTargets& jacobian_targets,
                               const AbsorptionBands& absorption_bands,
                               const AtmPoint& atm_point) {
  const auto n = arts_omp_get_max_threads();
  if (n == 1 or arts_omp_in_parallel() or n > f_grid.size()) {
    lbl::calculate(
        pm, dpm, f_grid, jacobian_targets, absorption_bands, atm_point);
  } else {
    const auto ompv = omp_offset_count(f_grid.size(), n);
    std::string error;
#pragma omp parallel for
    for (Index i = 0; i < n; i++) {
      try {
        lbl::calculate(
            pm.slice(ompv[i].first, ompv[i].second),
            dpm(joker,
                Range(ompv[i].first, ompv[i].second)),  // FIXME: IF DERIVS
            f_grid.slice(ompv[i].first, ompv[i].second),
            jacobian_targets,
            absorption_bands,
            atm_point);
      } catch (std::exception& e) {
#pragma omp critical
        if (error.empty()) error = var_string(e.what(), '\n');
      }
    }

    ARTS_USER_ERROR_IF(not error.empty(), error)
  }
}
