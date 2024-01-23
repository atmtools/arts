/**
 * @file m_linemixing.cc
 * @author Richard Larsson
 * @date 2020-06-23
 * 
 * @brief User interface for dealing with pure line mixing calculations
 * 
 * Note: If defined using parameterized form, the normal line-functions
 * approach is faster and more appropriate.  These functions should first
 * compute the relaxation, not simply use the relaxation
 */

#include <matpack.h>

#include <stdexcept>

#include "arts_conversions.h"
#include "atm.h"
#include "debug.h"
#include "hitran_species.h"
#include "isotopologues.h"
#include "jacobian.h"
#include "lbl_lineshape_linemixing.h"
#include "lbl_temperature_model.h"
#include "linemixing.h"
#include "linemixing_hitran.h"

void abs_hitran_relmat_dataReadHitranRelmatDataAndLines(
    HitranRelaxationMatrixData& abs_hitran_relmat_data,
    ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const String& basedir,
    const Numeric& linemixinglimit,
    const Numeric& fmin,
    const Numeric& fmax,
    const Numeric& stot,
    const String& mode) {
  const lm_hitran_2017::ModeOfLineMixing intmode =
      lm_hitran_2017::toModeOfLineMixingOrThrow(mode);

  const SpeciesIsotopologueRatios isotopologue_ratios =
      Hitran::isotopologue_ratios();

  ARTS_USER_ERROR_IF(abs_species.size() not_eq abs_lines_per_species.size(),
                     "Bad size of input species+lines");

  ArrayOfAbsorptionLines lines;
  lm_hitran_2017::read(abs_hitran_relmat_data,
                       lines,
                       isotopologue_ratios,
                       basedir,
                       linemixinglimit,
                       fmin,
                       fmax,
                       stot,
                       intmode);
  std::for_each(lines.begin(), lines.end(), [](auto& band) {
    band.sort_by_frequency();
  });  // Sort so largest frequency is last
  ArrayOfIndex used(lines.size(), false);

  bool emptied = false;
  for (Size i = 0; i < abs_species.size(); i++) {
    for (Size j = 0; j < abs_species[i].size(); j++) {
      if (abs_species[i][j].Spec() not_eq Species::fromShortName("CO2"))
        continue;

      if (not emptied) {
        abs_lines_per_species[i].resize(0);
        emptied = true;
      }

      for (Size k = 0; k < lines.size(); k++) {
        if (used[k]) continue;

        const Numeric lf{abs_species[i][j].lower_freq > 0
                             ? abs_species[i][j].lower_freq
                             : -std::numeric_limits<Numeric>::max()};
        const Numeric uf{abs_species[i][j].upper_freq > 0
                             ? abs_species[i][j].upper_freq
                             : std::numeric_limits<Numeric>::max()};

        // Select lines with correct Isotopologue and one line center within the range
        if ((abs_species[i][j].Isotopologue() == lines[k].Isotopologue() or
             abs_species[i][j].Isotopologue() ==
                 Species::select_joker("CO2")) and
            (lines[k].lines.front().F0 <= uf and
             lines[k].lines.back().F0 >= lf)) {
          used[k] = true;  // The lines should not be copied into other places
          abs_lines_per_species[i].push_back(lines[k]);
        }
      }
    }
  }
}

void propagation_matrixAddHitranLineMixingLines(
    PropmatVector& propagation_matrix,
    const HitranRelaxationMatrixData& abs_hitran_relmat_data,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const Vector& f_grid,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfSpeciesTag& select_abs_species,
    const JacobianTargets& jacobian_targets,
    const AtmPoint& atm_point) {
  ARTS_USER_ERROR_IF(jacobian_targets.target_count(),
                     "Cannot support any Jacobian at this time");
  ARTS_USER_ERROR_IF(abs_species.size() not_eq abs_lines_per_species.size(),
                     "Bad size of input species+lines");

  // vmrs should be [air, water, co2]  FIXME: confirm, because code disagreed
  const Numeric water = atm_point[Species::Species::Water];
  const Numeric co2 = atm_point[Species::Species::CarbonDioxide];
  const Vector vmrs{co2, water, 1.0 - co2 - water};

  for (Size i = 0; i < abs_species.size(); i++) {
    if (select_abs_species.size() and select_abs_species not_eq abs_species[i])
      continue;

    if (abs_lines_per_species[i].size() and
        (abs_lines_per_species[i].front().population ==
             Absorption::PopulationType::ByHITRANFullRelmat or
         abs_lines_per_species[i].front().population ==
             Absorption::PopulationType::ByHITRANRosenkranzRelmat)) {
      // const auto compres = lm_hitran_2017::compute(
      //     abs_hitran_relmat_data, abs_lines_per_species[i], isotopologue_ratios,
      //     atm_point.pressure, atm_point.temperature, vmrs, f_grid);
      // for (Index iv = 0; iv < f_grid.nelem(); iv++) {
      //   propagation_matrix[iv].A() += compres[iv];
      // }
      ARTS_ASSERT(false, "Fix isotopologues")
    }
  }
}

void abs_lines_per_speciesAdaptHitranLineMixing(
    ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const HitranRelaxationMatrixData& abs_hitran_relmat_data,
    const Vector& t_grid,
    const Numeric& pressure,
    const Index& order) {
  for (auto& abs_lines : abs_lines_per_species) {
    for (auto& band : abs_lines) {
      if (band.population == Absorption::PopulationType::ByHITRANFullRelmat or
          band.population ==
              Absorption::PopulationType::ByHITRANRosenkranzRelmat) {
        lm_hitran_2017::hitran_lm_eigenvalue_adaptation(
            band, t_grid, abs_hitran_relmat_data, pressure, order);
      }
    }
  }
}

void propagation_matrixAddOnTheFlyLineMixing(
    PropmatVector& propagation_matrix,
    PropmatMatrix& propagation_matrix_jacobian,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const MapOfErrorCorrectedSuddenData& ecs_data,
    const Vector& f_grid,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfSpeciesTag& select_abs_species,
    const JacobianTargets& jacobian_targets,
    const AtmPoint& atm_point) {
  ARTS_USER_ERROR_IF(abs_species.size() not_eq abs_lines_per_species.size(),
                     "Bad size of input species+lines");

  for (Size i = 0; i < abs_species.size(); i++) {
    if (select_abs_species.size() and select_abs_species not_eq abs_species[i])
      continue;
    for (auto& band : abs_lines_per_species[i]) {
      if (band.OnTheFlyLineMixing() and band.DoLineMixing(atm_point.pressure)) {
        // vmrs should be for the line
        const Vector line_shape_vmr = band.BroadeningSpeciesVMR(atm_point);
        const Numeric this_vmr = atm_point[abs_species[i].Species()] *
                                 atm_point[band.Isotopologue()];
        const auto [abs, dabs, error] = Absorption::LineMixing::ecs_absorption(
            atm_point.temperature,
            0,
            atm_point.pressure,
            this_vmr,
            line_shape_vmr,
            ecs_data[band.quantumidentity],
            f_grid,
            Zeeman::Polarization::None,
            band,
            jacobian_targets);
        for (Index iv = 0; iv < f_grid.nelem(); iv++) {
          propagation_matrix[iv].A() += abs[iv].real();
        }

        for (auto& atm : jacobian_targets.atm()) {
          const auto j = atm.target_pos;
          if (atm.type == abs_species[i].Species()) {
            for (Index iv = 0; iv < f_grid.nelem(); iv++) {
              propagation_matrix_jacobian[j][iv].A() += abs[iv].real();
            }
          } else {
            for (Index iv = 0; iv < f_grid.nelem(); iv++) {
              propagation_matrix_jacobian[j][iv].A() += dabs[j][iv].real();
            }
          }
        }

        // FIXME: Add line parameters!
      }
    }
  }
}

void abs_linesAdaptOnTheFlyLineMixing(
    ArrayOfAbsorptionLines& abs_lines,
    const MapOfErrorCorrectedSuddenData& ecs_data,
    const Vector& t_grid,
    const Numeric& pressure,
    const Index& order,
    const Index& robust,
    const Index& rosenkranz_adaptation) {
  for (auto& band : abs_lines) {
    if (band.population ==
            Absorption::PopulationType::ByRovibLinearDipoleLineMixing or
        band.population == Absorption::PopulationType::ByMakarovFullRelmat) {
      Absorption::LineMixing::ecs_eigenvalue_adaptation(
          band,
          t_grid,
          ecs_data[band.quantumidentity],
          pressure,
          order,
          robust,
          rosenkranz_adaptation);
    }
  }
}

void abs_lines_per_speciesAdaptOnTheFlyLineMixing(
    ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const MapOfErrorCorrectedSuddenData& ecs_data,
    const Vector& t_grid,
    const Numeric& pressure,
    const Index& order,
    const Index& robust,
    const Index& rosenkranz_adaptation) {
  for (auto& abs_lines : abs_lines_per_species) {
    abs_linesAdaptOnTheFlyLineMixing(abs_lines,
                                     ecs_data,
                                     t_grid,
                                     pressure,
                                     order,
                                     robust,
                                     rosenkranz_adaptation);
  }
}

void ecs_dataInit(MapOfErrorCorrectedSuddenData& ecs_data) {
  ecs_data.resize(0);
}

void ecs_dataAddSpeciesData(MapOfErrorCorrectedSuddenData& ecs_data,
                            const AtmField& atm_field,
                            const QuantumIdentifier& qid,
                            const String& species,
                            const String& scaling_type,
                            const Vector& scaling,
                            const String& beta_type,
                            const Vector& beta,
                            const String& lambda_type,
                            const Vector& lambda,
                            const String& collisional_distance_type,
                            const Vector& collisional_distance) {
  ARTS_ASSERT(false, "Fix isotopologues")

  const Species::Species spec = Species::fromShortName(species);
  ARTS_USER_ERROR_IF(not good_enum(spec), "Invalid species: ", species)
  auto& data = ecs_data[qid][spec];
  data.scaling = LineShapeModelParameters(
      LineShape::toTemperatureModelOrThrow(scaling_type), scaling);
  data.beta = LineShapeModelParameters(
      LineShape::toTemperatureModelOrThrow(beta_type), beta);
  data.lambda = LineShapeModelParameters(
      LineShape::toTemperatureModelOrThrow(lambda_type), lambda);
  data.collisional_distance = LineShapeModelParameters(
      LineShape::toTemperatureModelOrThrow(collisional_distance_type),
      collisional_distance);
  //data.mass = Species::mean_mass(spec, isotopologue_ratios);
  ARTS_ASSERT(false, "Fix isotopologues")
}

void ecs_dataAddMeanAir(MapOfErrorCorrectedSuddenData& ecs_data,
                        const Vector& vmrs,
                        const ArrayOfSpeciesTag& specs) {
  ARTS_USER_ERROR_IF(specs.size() not_eq static_cast<Size>(vmrs.nelem()),
                     "Bad sizes of specs and vmrs\nspecs: [",
                     specs,
                     "]\nvmrs: [",
                     vmrs,
                     "]")

  static_assert(LineShapeModelParameters::N == 4);
  const auto scale = [](LineShapeModelParameters& mp, Numeric x) {
    mp.X0 *= x;
    mp.X1 *= x;
    mp.X2 *= x;
    mp.X3 *= x;
  };
  const auto scale_and_add_or_throw = [](LineShapeModelParameters& mp1,
                                         const LineShapeModelParameters& mp2,
                                         Numeric x) {
    ARTS_USER_ERROR_IF(mp1.type not_eq mp2.type,
                       "Can only scale and add same type\nmp1: ",
                       mp1,
                       "\nmp2: ",
                       mp2)
    mp1.X0 += mp2.X0 * x;
    mp1.X1 += mp2.X1 * x;
    mp1.X2 += mp2.X2 * x;
    mp1.X3 += mp2.X3 * x;
  };

  for (auto& ecs : ecs_data) {
    SpeciesErrorCorrectedSuddenData airdata;

    bool found = false;
    Numeric sumvmr = 0;
    for (Size i = 0; i < specs.size(); i++) {
      ARTS_USER_ERROR_IF(not specs[i].Isotopologue().joker(),
                         "Can only have joker species, finds: [",
                         specs,
                         ']')

      auto spec = specs[i].Spec();

      if (spec == Species::Species::Bath) continue;

      auto& data = ecs[spec];

      // This means we never had the data so we should skip
      if (std::isinf(data.mass)) {
        continue;
      }

      const Numeric vmr = vmrs[i];
      if (not found) {
        airdata.scaling = data.scaling;
        airdata.beta = data.beta;
        airdata.lambda = data.lambda;
        airdata.collisional_distance = data.collisional_distance;
        airdata.mass = data.mass * vmr;
        scale(airdata.scaling, vmr);
        scale(airdata.beta, vmr);
        scale(airdata.lambda, vmr);
        scale(airdata.collisional_distance, vmr);

        found = true;
      } else {
        scale_and_add_or_throw(airdata.scaling, data.scaling, vmr);
        scale_and_add_or_throw(airdata.beta, data.beta, vmr);
        scale_and_add_or_throw(airdata.lambda, data.lambda, vmr);
        scale_and_add_or_throw(
            airdata.collisional_distance, data.collisional_distance, vmr);
        airdata.mass += data.mass * vmr;
      }

      sumvmr += vmrs[i];
    }

    if (sumvmr not_eq 0 and sumvmr not_eq 1) {
      scale(airdata.scaling, 1.0 / sumvmr);
      scale(airdata.beta, 1.0 / sumvmr);
      scale(airdata.lambda, 1.0 / sumvmr);
      scale(airdata.collisional_distance, 1.0 / sumvmr);
      airdata.mass /= sumvmr;
    }

    ecs[Species::Species::Bath] = airdata;
  }
}

void ecs_dataAddMeanAirNEWNEW(LinemixingEcsData& ecs_data,
                              const Vector& vmrs,
                              const ArrayOfSpecies& specs) {
  ARTS_USER_ERROR_IF(static_cast<Size>(vmrs.size()) != specs.size(),
                     "Mismatch dimension of vmrs and specs")
  ARTS_USER_ERROR_IF(
      std::abs(sum(vmrs) - 1.) > 1e-4, "Bad vmrs [sum far from 1]: ", vmrs)

  auto set = [](const lbl::temperature::data& data1,
                const lbl::temperature::data& data2,
                const Numeric vmr,
                bool first) -> lbl::temperature::data {
    auto v = data1.X();
    v *= vmr;

    if (first) {
      return {data1.Type(), v};
    }

    v += data2.X();
    ARTS_USER_ERROR_IF(data1.Type() != data2.Type(),
                       "Type error in species temperature type for ",
                       data1.Type(),
                       " and ",
                       data2.Type())
    return {data1.Type(), v};
  };

  for (auto& [isot, data] : ecs_data) {
    auto& airdata = data[Species::Species::Bath];

    bool first = true;
    for (Index i = 0; i < vmrs.size(); i++) {
      const auto spec = specs[i];
      const auto vmr = vmrs[i];

      try {
        auto& specdata = data.at(spec);
        airdata.scaling = set(specdata.scaling, airdata.scaling, vmr, first);
        airdata.beta = set(specdata.beta, airdata.beta, vmr, first);
        airdata.lambda = set(specdata.lambda, airdata.lambda, vmr, first);
        airdata.collisional_distance = set(specdata.collisional_distance,
                                           airdata.collisional_distance,
                                           vmr,
                                           first);
        first = false;
      } catch (std::out_of_range&) {
        ARTS_USER_ERROR(
            "Missing species ", spec, " in ecs_data of isotopologue ", isot)
      } catch (std::exception& e) {
        ARTS_USER_ERROR("Error for species ",
                        spec,
                        " in ecs_data of isotopologue ",
                        isot,
                        ":\n",
                        e.what())
      }
    }
  }
}

void ecs_dataAddMakarov2020(MapOfErrorCorrectedSuddenData& ecs_data,
                            const AtmField& isotopologue_ratios) {
  ARTS_ASSERT(false, "Fix isotopologues")
  // The band is ignored
  auto& ecs = ecs_data[QuantumIdentifier("O2-66")];

  // All species have the same effect, so just copy the values but change the mass (allow new mass for Air)

  ecs[Species::Species::Oxygen].scaling =
      LineShapeModelParameters(LineShapeTemperatureModel::T0, 1.0, 0, 0, 0);
  ecs[Species::Species::Oxygen].collisional_distance = LineShapeModelParameters(
      LineShapeTemperatureModel::T0, Conversion::angstrom2meter(0.61), 0, 0, 0);
  ecs[Species::Species::Oxygen].lambda =
      LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.39, 0, 0, 0);
  ecs[Species::Species::Oxygen].beta =
      LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.567, 0, 0, 0);
  // ecs[Species::Species::Oxygen].mass = 31.9898;

  ecs[Species::Species::Nitrogen].scaling =
      LineShapeModelParameters(LineShapeTemperatureModel::T0, 1.0, 0, 0, 0);
  ecs[Species::Species::Nitrogen].collisional_distance =
      LineShapeModelParameters(LineShapeTemperatureModel::T0,
                               Conversion::angstrom2meter(0.61),
                               0,
                               0,
                               0);
  ecs[Species::Species::Nitrogen].lambda =
      LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.39, 0, 0, 0);
  ecs[Species::Species::Nitrogen].beta =
      LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.567, 0, 0, 0);
  // ecs[Species::Species::Nitrogen].mass = 28.8504;
}

void ecs_dataInitNEWNEW(LinemixingEcsData& ecs_data) { ecs_data.clear(); }

void ecs_dataAddMakarov2020NEWNEW(LinemixingEcsData& ecs_data) {
  using enum lbl::temperature::model_type;
  using data = lbl::temperature::data;

  auto& ecs = ecs_data[SpeciesIsotopeRecord("O2-66")];

  // All species have the same effect, so just copy the values but change the mass (allow new mass for Air)
  auto& oxy = ecs[Species::Species::Oxygen];
  oxy.scaling = data(T0, {1.0});
  oxy.collisional_distance = data(T0, {Conversion::angstrom2meter(0.61)});
  oxy.lambda = data(T0, {0.39});
  oxy.beta = data(T0, {0.567});

  auto& nit = ecs[Species::Species::Nitrogen];
  nit.scaling = data(T0, {1.0});
  nit.collisional_distance = data(T0, {Conversion::angstrom2meter(0.61)});
  nit.lambda = data(T0, {0.39});
  nit.beta = data(T0, {0.567});
}

void ecs_dataAddRodrigues1997(MapOfErrorCorrectedSuddenData& ecs_data,
                              const AtmField& isotopologue_ratios) {
  ARTS_ASSERT(false, "Fix isotopologues")
  for (const auto* key : {"CO2-626", "CO2-628", "CO2-636"}) {
    auto& ecs = ecs_data[QuantumIdentifier(key)];

    ecs[Species::Species::Nitrogen].scaling =
        LineShapeModelParameters(LineShapeTemperatureModel::T1,
                                 Conversion::kaycm_per_atm2hz_per_pa(0.0180),
                                 0.85,
                                 0,
                                 0);
    ecs[Species::Species::Nitrogen].lambda = LineShapeModelParameters(
        LineShapeTemperatureModel::T1, 0.81, 0.0152, 0, 0);
    ecs[Species::Species::Nitrogen].beta =
        LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.008, 0, 0, 0);
    ecs[Species::Species::Nitrogen].collisional_distance =
        LineShapeModelParameters(LineShapeTemperatureModel::T0,
                                 Conversion::angstrom2meter(2.2),
                                 0,
                                 0,
                                 0);
    // ecs[Species::Species::Nitrogen].mass =
    //     Species::mean_mass(Species::Species::Nitrogen, isotopologue_ratios);

    ecs[Species::Species::Oxygen].scaling =
        LineShapeModelParameters(LineShapeTemperatureModel::T1,
                                 Conversion::kaycm_per_atm2hz_per_pa(0.0168),
                                 0.5,
                                 0,
                                 0);
    ecs[Species::Species::Oxygen].lambda = LineShapeModelParameters(
        LineShapeTemperatureModel::T1, 0.82, -0.091, 0, 0);
    ecs[Species::Species::Oxygen].beta =
        LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.007, 0, 0, 0);
    ecs[Species::Species::Oxygen].collisional_distance =
        LineShapeModelParameters(LineShapeTemperatureModel::T0,
                                 Conversion::angstrom2meter(2.4),
                                 0,
                                 0,
                                 0);
    // ecs[Species::Species::Oxygen].mass =
    //     Species::mean_mass(Species::Species::Oxygen, isotopologue_ratios);
  }
}

void ecs_dataAddRodrigues1997NEWNEW(LinemixingEcsData& ecs_data) {
  using enum lbl::temperature::model_type;
  using data = lbl::temperature::data;

  for (const std::string_view key : {"CO2-626", "CO2-628", "CO2-636"}) {
    auto& ecs = ecs_data[SpeciesIsotopeRecord(key)];

    ecs[Species::Species::Nitrogen].scaling =
        data(T1, {Conversion::kaycm_per_atm2hz_per_pa(0.0180), 0.85});
    ecs[Species::Species::Nitrogen].lambda = data(T1, {.81, 0.0152});
    ecs[Species::Species::Nitrogen].beta = data(T0, {.008});
    ecs[Species::Species::Nitrogen].collisional_distance =
        data(T0, {Conversion::angstrom2meter(2.2)});

    ecs[Species::Species::Oxygen].scaling =
        data(T1, {Conversion::kaycm_per_atm2hz_per_pa(0.0168), 0.5});
    ecs[Species::Species::Oxygen].lambda = data(T1, {.82, -0.091});
    ecs[Species::Species::Oxygen].beta = data(T0, {.007});
    ecs[Species::Species::Oxygen].collisional_distance =
        data(T0, {Conversion::angstrom2meter(2.4)});
  }
}

void ecs_dataAddTran2006(MapOfErrorCorrectedSuddenData& ecs_data,
                         const AtmField& isotopologue_ratios) {
  ARTS_ASSERT(false, "Fix isotopologues")
  auto& ecs = ecs_data[QuantumIdentifier("O2-66")];

  ecs[Species::Species::Oxygen].scaling =
      LineShapeModelParameters(LineShapeTemperatureModel::T1,
                               Conversion::kaycm_per_atm2hz_per_pa(0.0275),
                               1.01,
                               0,
                               0);
  ecs[Species::Species::Oxygen].collisional_distance = LineShapeModelParameters(
      LineShapeTemperatureModel::T0, Conversion::angstrom2meter(1.05), 0, 0, 0);
  ecs[Species::Species::Oxygen].lambda =
      LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.935, 0, 0, 0);
  ecs[Species::Species::Oxygen].beta =
      LineShapeModelParameters(LineShapeTemperatureModel::T0, 0, 0, 0, 0);
  // ecs[Species::Species::Oxygen].mass =
  //     Species::mean_mass(Species::Species::Oxygen, isotopologue_ratios);

  ecs[Species::Species::Nitrogen].scaling =
      LineShapeModelParameters(LineShapeTemperatureModel::T1,
                               Conversion::kaycm_per_atm2hz_per_pa(0.0285),
                               1.03,
                               0,
                               0);
  ecs[Species::Species::Nitrogen].collisional_distance =
      LineShapeModelParameters(LineShapeTemperatureModel::T0,
                               Conversion::angstrom2meter(1.00),
                               0,
                               0,
                               0);
  ecs[Species::Species::Nitrogen].lambda =
      LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.95, 0, 0, 0);
  ecs[Species::Species::Nitrogen].beta =
      LineShapeModelParameters(LineShapeTemperatureModel::T0, 0, 0, 0, 0);
  // ecs[Species::Species::Nitrogen].mass =
  //     Species::mean_mass(Species::Species::Nitrogen, isotopologue_ratios);
}

void ecs_dataAddTran2011(MapOfErrorCorrectedSuddenData& ecs_data,
                         const AtmField& isotopologue_ratios) {
  ARTS_ASSERT(false, "Fix isotopologues")
  for (const auto* key : {"CO2-626", "CO2-628", "CO2-636"}) {
    auto& ecs = ecs_data[QuantumIdentifier(key)];

    ecs[Species::Species::CarbonDioxide].scaling =
        LineShapeModelParameters(LineShapeTemperatureModel::T0,
                                 Conversion::kaycm_per_atm2hz_per_pa(0.019),
                                 0,
                                 0,
                                 0);
    ecs[Species::Species::CarbonDioxide].lambda =
        LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.61, 0, 0, 0);
    ecs[Species::Species::CarbonDioxide].beta =
        LineShapeModelParameters(LineShapeTemperatureModel::T0, 0.052, 0, 0, 0);
    ecs[Species::Species::CarbonDioxide].collisional_distance =
        LineShapeModelParameters(LineShapeTemperatureModel::T0,
                                 Conversion::angstrom2meter(5.5),
                                 0,
                                 0,
                                 0);
    // ecs[Species::Species::CarbonDioxide].mass = Species::mean_mass(
    //     Species::Species::CarbonDioxide, isotopologue_ratios);
  }
}

void ecs_dataAddTran2011NEWNEW(LinemixingEcsData& ecs_data) {
  using enum lbl::temperature::model_type;
  using data = lbl::temperature::data;

  for (const std::string_view key : {"CO2-626", "CO2-628", "CO2-636"}) {
    auto& ecs = ecs_data[SpeciesIsotopeRecord(key)];

    ecs[Species::Species::CarbonDioxide].scaling =
        data(T0, {Conversion::kaycm_per_atm2hz_per_pa(0.019)});
    ecs[Species::Species::CarbonDioxide].lambda = data(T0, {0.61});
    ecs[Species::Species::CarbonDioxide].beta = data(T0, {0.052});
    ecs[Species::Species::CarbonDioxide].collisional_distance =
        data(T0, {Conversion::angstrom2meter(5.5)});
  }
}
