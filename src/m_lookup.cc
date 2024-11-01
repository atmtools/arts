#include <jacobian.h>
#include <lookup_map.h>

#include <algorithm>
#include <ranges>

#include "mh_checks.h"

void absorption_lookup_tableInit(
    AbsorptionLookupTables& absorption_lookup_table) {
  absorption_lookup_table.clear();
}

template <bool calc>
std::conditional_t<calc, Vector, void> _propagation_matrixAddLookup(
    PropmatVector& propagation_matrix [[maybe_unused]],
    PropmatMatrix& propagation_matrix_jacobian [[maybe_unused]],
    const AscendingGrid& frequency_grid,
    const JacobianTargets& jacobian_targets [[maybe_unused]],
    const SpeciesEnum& propagation_matrix_select_species,
    const AbsorptionLookupTables& absorption_lookup_table,
    const AtmPoint& atmospheric_point,
    const Index& no_negative_absorption,
    const Index& p_interp_order,
    const Index& t_interp_order,
    const Index& water_interp_order,
    const Index& f_interp_order,
    const Numeric& extpolfac) {
  if constexpr (calc) {
    Vector absorption(frequency_grid.size(), 0.0);

    if (propagation_matrix_select_species == "Bath"_spec) {
      for (auto& [spec, data] : absorption_lookup_table) {
        data.absorption(absorption,
                        spec,
                        p_interp_order,
                        t_interp_order,
                        water_interp_order,
                        f_interp_order,
                        atmospheric_point,
                        frequency_grid,
                        extpolfac);
      }
    } else {
      absorption_lookup_table.at(propagation_matrix_select_species)
          .absorption(absorption,
                      propagation_matrix_select_species,
                      p_interp_order,
                      t_interp_order,
                      water_interp_order,
                      f_interp_order,
                      atmospheric_point,
                      frequency_grid,
                      extpolfac);
    }

    return absorption;
  } else {
    const Vector absorption = _propagation_matrixAddLookup<not calc>(
        propagation_matrix,
        propagation_matrix_jacobian,
        frequency_grid,
        jacobian_targets,
        propagation_matrix_select_species,
        absorption_lookup_table,
        atmospheric_point,
        no_negative_absorption,
        p_interp_order,
        t_interp_order,
        water_interp_order,
        f_interp_order,
        extpolfac);

    for (Index i = 0; i < frequency_grid.size(); i++) {
      if (no_negative_absorption == 0 or absorption[i] > 0.0) {
        propagation_matrix[i].A() += absorption[i];
      }
    }

    if (jacobian_targets.atm().size()) {
      Vector d_absorption;

      for (auto& jacobian_target : jacobian_targets.atm()) {
        ARTS_USER_ERROR_IF(
            not std::isnormal(jacobian_target.d),
            "The target {} is not good, it lacks a perturbation value.",
            jacobian_target);

        if (jacobian_target.is_wind()) {
          const AscendingGrid frequency_grid2(
              frequency_grid.begin(),
              frequency_grid.end(),
              [d = jacobian_target.d](Numeric x) { return x + d; });
          d_absorption = _propagation_matrixAddLookup<not calc>(
              propagation_matrix,
              propagation_matrix_jacobian,
              frequency_grid2,
              {},
              propagation_matrix_select_species,
              absorption_lookup_table,
              atmospheric_point,
              no_negative_absorption,
              p_interp_order,
              t_interp_order,
              water_interp_order,
              f_interp_order,
              extpolfac);

        } else {
          AtmPoint atm_point               = atmospheric_point;
          atm_point[jacobian_target.type] += jacobian_target.d;

          d_absorption = _propagation_matrixAddLookup<not calc>(
              propagation_matrix,
              propagation_matrix_jacobian,
              frequency_grid,
              {},
              propagation_matrix_select_species,
              absorption_lookup_table,
              atm_point,
              no_negative_absorption,
              p_interp_order,
              t_interp_order,
              water_interp_order,
              f_interp_order,
              extpolfac);
        }

        const Numeric d_inv = 1.0 / jacobian_target.d;
        for (Index i = 0; i < frequency_grid.size(); i++) {
          if (no_negative_absorption == 0 or d_absorption[i] > 0.0) {
            propagation_matrix_jacobian(jacobian_target.target_pos, i).A() =
                (d_absorption[i] - absorption[i]) * d_inv;
          }
        }
      }
    }
  }
}

void propagation_matrixAddLookup(
    PropmatVector& propagation_matrix,
    PropmatMatrix& propagation_matrix_jacobian,
    const AscendingGrid& frequency_grid,
    const JacobianTargets& jacobian_targets,
    const SpeciesEnum& propagation_matrix_select_species,
    const AbsorptionLookupTables& absorption_lookup_table,
    const AtmPoint& atmospheric_point,
    const Index& no_negative_absorption,
    const Index& p_interp_order,
    const Index& t_interp_order,
    const Index& water_interp_order,
    const Index& f_interp_order,
    const Numeric& extpolfac) try {
  _propagation_matrixAddLookup<false>(propagation_matrix,
                                      propagation_matrix_jacobian,
                                      frequency_grid,
                                      jacobian_targets,
                                      propagation_matrix_select_species,
                                      absorption_lookup_table,
                                      atmospheric_point,
                                      no_negative_absorption,
                                      p_interp_order,
                                      t_interp_order,
                                      water_interp_order,
                                      f_interp_order,
                                      extpolfac);
}
ARTS_METHOD_ERROR_CATCH

void absorption_lookup_tablePrecompute(
    AbsorptionLookupTables& absorption_lookup_table,
    const ArrayOfAtmPoint& ray_path_atmospheric_point,
    const AscendingGrid& frequency_grid,
    const AbsorptionBands& absorption_bands,
    const LinemixingEcsData& ecs_data,
    const AscendingGrid& temperature_perturbation,
    const AscendingGrid& water_perturbation,
    const SpeciesEnum& select_species) {
  absorption_lookup_table[select_species] = AbsorptionLookupTable(
      select_species,
      ray_path_atmospheric_point,
      std::make_shared<const AscendingGrid>(frequency_grid),
      absorption_bands,
      ecs_data,
      std::make_shared<const AscendingGrid>(temperature_perturbation),
      std::make_shared<const AscendingGrid>(water_perturbation));
}

void absorption_lookup_tablePrecomputeAll(
    AbsorptionLookupTables& absorption_lookup_table,
    const ArrayOfAtmPoint& ray_path_atmospheric_point,
    const AscendingGrid& frequency_grid,
    const AbsorptionBands& absorption_bands,
    const LinemixingEcsData& ecs_data,
    const AscendingGrid& temperature_perturbation,
    const AscendingGrid& water_perturbation,
    const ArrayOfSpeciesEnum& water_affected_species) {
  const AscendingGrid empty{};

  ArrayOfSpeciesEnum lut_species;
  const auto species_not_in_lut =
      std::views::transform(
          [](const auto& pair) { return pair.first.Species(); }) |
      std::views::filter([&lut_species](const SpeciesEnum& s) {
        return lut_species.end() == std::ranges::find(lut_species, s);
      });

  std::ranges::copy(absorption_bands | species_not_in_lut,
                    std::back_inserter(lut_species));

  std::ranges::sort(lut_species);

  for (auto spec : water_affected_species) {
    ARTS_USER_ERROR_IF(
        not std::ranges::binary_search(lut_species, spec),
        R"(Missing a species in the absorption bands that is marked as affected by water vapor.
  Species:                                       {}
  Absorption band species:
    {:B,}
  All species marked as affected by water vapor:
    {:B,}
)",
        spec,
        lut_species,
        water_affected_species)
  }

  for (SpeciesEnum spec : lut_species) {
    const bool do_water_perturb =
        std::ranges::any_of(water_affected_species, Cmp::eq(spec));

    absorption_lookup_tablePrecompute(
        absorption_lookup_table,
        ray_path_atmospheric_point,
        frequency_grid,
        absorption_bands,
        ecs_data,
        temperature_perturbation,
        do_water_perturb ? water_perturbation : empty,
        spec);
  }

  //! Make the newly added LUT entries share grids with eachother
  if (lut_species.size() > 1) {
    const std::shared_ptr<const AscendingGrid> fs =
        absorption_lookup_table[lut_species.front()].f_grid;
    const std::shared_ptr<const DescendingGrid> ps =
        absorption_lookup_table[lut_species.front()].log_p_grid;
    const std::shared_ptr<const AscendingGrid> ts =
        absorption_lookup_table[lut_species.front()].t_pert;
    const std::shared_ptr<const AscendingGrid> ws =
        water_affected_species.empty()
            ? nullptr
            : absorption_lookup_table[water_affected_species.front()].w_pert;

    for (SpeciesEnum spec : lut_species | std::ranges::views::drop(1)) {
      const bool do_water_perturb =
          std::ranges::any_of(water_affected_species, Cmp::eq(spec));

      absorption_lookup_table[spec].f_grid     = fs;
      absorption_lookup_table[spec].log_p_grid = ps;
      absorption_lookup_table[spec].t_pert     = ts;
      if (do_water_perturb) absorption_lookup_table[spec].w_pert = ws;
    }
  }
}

void absorption_lookup_tableFromProfiles(
    AbsorptionLookupTables& absorption_lookup_table,
    const AscendingGrid& frequency_grid,
    const AbsorptionBands& absorption_bands,
    const LinemixingEcsData& ecs_data,
    const DescendingGrid& pressure_profile,
    const Vector& temperature_profile,
    const SpeciesEnumVectors& vmr_profiles,
    const AscendingGrid& temperature_perturbation,
    const AscendingGrid& water_perturbation,
    const ArrayOfSpeciesEnum& water_affected_species,
    const String& isoratio_option) {
  absorption_lookup_tableInit(absorption_lookup_table);

  ArrayOfAtmPoint ray_path_atmospheric_point(
      pressure_profile.size(), AtmPoint{to<IsoRatioOption>(isoratio_option)});

  ARTS_USER_ERROR_IF(not same_shape(pressure_profile, temperature_profile),
                     "Pressure and temperature profiles must agree in size.");

  for (Index i = 0; i < pressure_profile.size(); i++) {
    ray_path_atmospheric_point[i].pressure    = pressure_profile[i];
    ray_path_atmospheric_point[i].temperature = temperature_profile[i];
  }

  for (auto& [spec, prof] : vmr_profiles) {
    ARTS_USER_ERROR_IF(
        not same_shape(pressure_profile, prof),
        "Pressure and VMR profiles must agree in size, fails for species {}",
        spec);

    for (Index i = 0; i < prof.size(); i++) {
      ray_path_atmospheric_point[i][spec] = prof[i];
    }
  }

  absorption_lookup_tablePrecomputeAll(absorption_lookup_table,
                                       ray_path_atmospheric_point,
                                       frequency_grid,
                                       absorption_bands,
                                       ecs_data,
                                       temperature_perturbation,
                                       water_perturbation,
                                       water_affected_species);
}
