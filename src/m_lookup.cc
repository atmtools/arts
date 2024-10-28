#include <jacobian.h>
#include <lookup_map.h>

#include <algorithm>
#include <set>

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
    const Numeric& extpolfac) {
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
      std::make_shared<AscendingGrid>(frequency_grid),
      std::make_shared<ArrayOfAtmPoint>(ray_path_atmospheric_point),
      absorption_bands,
      ecs_data,
      temperature_perturbation,
      water_perturbation);
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
  std::set<SpeciesEnum> species_set;

  for (auto& [qid, _] : absorption_bands) {
    if (const SpeciesEnum s = qid.Species(); not species_set.contains(s)) {
      species_set.insert(s);
      const AscendingGrid local_water_perturbation =
          std::ranges::any_of(water_affected_species, Cmp::eq(s))
              ? water_perturbation
              : AscendingGrid{};
      absorption_lookup_tablePrecompute(absorption_lookup_table,
                                        ray_path_atmospheric_point,
                                        frequency_grid,
                                        absorption_bands,
                                        ecs_data,
                                        temperature_perturbation,
                                        local_water_perturbation,
                                        s);
    }
  }
}

void absorption_lookup_tableInit(
    AbsorptionLookupTables& absorption_lookup_table) {
  absorption_lookup_table = {};
}