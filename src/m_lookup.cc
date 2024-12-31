#include <jacobian.h>
#include <lookup_map.h>

#include <algorithm>
#include <ranges>
#include <set>

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

    for (Size i = 0; i < frequency_grid.size(); i++) {
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
        for (Size i = 0; i < frequency_grid.size(); i++) {
          if (no_negative_absorption == 0 or d_absorption[i] > 0.0) {
            propagation_matrix_jacobian[jacobian_target.target_pos, i].A() =
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

void ray_path_atmospheric_pointExtendInPressure(
    ArrayOfAtmPoint& ray_path_atmospheric_point,
    const Numeric& extended_max_pressure,
    const Numeric& extended_min_pressure,
    const String& extrapolation_option) try {
  lookup::extend_atmosphere(
      ray_path_atmospheric_point,
      to<InterpolationExtrapolation>(extrapolation_option),
      extended_max_pressure,
      extended_min_pressure);
}
ARTS_METHOD_ERROR_CATCH

void absorption_lookup_tablePrecompute(
    AbsorptionLookupTables& absorption_lookup_table,
    const ArrayOfAtmPoint& ray_path_atmospheric_point,
    const AscendingGrid& frequency_grid,
    const AbsorptionBands& absorption_bands,
    const LinemixingEcsData& ecs_data,
    const SpeciesEnum& select_species,
    const AscendingGrid& temperature_perturbation,
    const AscendingGrid& water_perturbation) {
  absorption_lookup_table[select_species] = {
      select_species,
      ray_path_atmospheric_point,
      std::make_shared<const AscendingGrid>(frequency_grid),
      absorption_bands,
      ecs_data,
      temperature_perturbation.empty()
          ? nullptr
          : std::make_shared<const AscendingGrid>(temperature_perturbation),
      water_perturbation.empty()
          ? nullptr
          : std::make_shared<const AscendingGrid>(water_perturbation)};
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

  const auto f = std::make_shared<const AscendingGrid>(frequency_grid);
  const auto t =
      std::make_shared<const AscendingGrid>(temperature_perturbation);
  const auto w =
      water_affected_species.empty()
          ? nullptr
          : std::make_shared<const AscendingGrid>(water_perturbation);

  for (SpeciesEnum s : lut_species) {
    const bool do_water_perturb =
        std::ranges::any_of(water_affected_species, Cmp::eq(s));

    if (do_water_perturb) {
      absorption_lookup_table[s] = {
          s, ray_path_atmospheric_point, f, absorption_bands, ecs_data, t, w};
    } else {
      absorption_lookup_table[s] = {
          s, ray_path_atmospheric_point, f, absorption_bands, ecs_data, t};
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

  ARTS_USER_ERROR_IF(not same_shape<1>(pressure_profile.vec(), temperature_profile),
                     "Pressure and temperature profiles must agree in size.");

  for (Size i = 0; i < pressure_profile.size(); i++) {
    ray_path_atmospheric_point[i].pressure    = pressure_profile[i];
    ray_path_atmospheric_point[i].temperature = temperature_profile[i];
  }

  for (auto& [spec, prof] : vmr_profiles) {
    ARTS_USER_ERROR_IF(
        not same_shape<1>(pressure_profile.vec(), prof),
        "Pressure and VMR profiles must agree in size, fails for species {}",
        spec);

    for (Size i = 0; i < prof.size(); i++) {
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

void absorption_lookup_tableSimpleWide(
    AbsorptionLookupTables& absorption_lookup_table,
    const AscendingGrid& frequency_grid,
    const AbsorptionBands& absorption_bands,
    const LinemixingEcsData& ecs_data,
    const ArrayOfSpeciesEnum& water_affected_species,
    const Vector2& pressure_range,
    const Vector2& temperature_range,
    const Vector2& water_vmr_range,
    const String& isoratio_option,
    const Numeric& vmr_value,
    const Index& atmospheric_steps,
    const Index& temperature_perturbation_steps,
    const Index& water_vmr_perturbation_steps) {
  ARTS_USER_ERROR_IF(pressure_range[1] <= pressure_range[0],
                     "Pressure range must be increasing, got {:B,}",
                     pressure_range);
  ARTS_USER_ERROR_IF(temperature_range[1] <= temperature_range[0],
                     "Temperature range must be increasing, got {:B,}",
                     temperature_range);
  ARTS_USER_ERROR_IF(water_vmr_range[1] <= water_vmr_range[0],
                     "Water VMR range must be increasing, got {:B,}",
                     water_vmr_range);
  ARTS_USER_ERROR_IF(
      vmr_value <= 0.0, "VMR value must be positive, got {}", vmr_value);
  ARTS_USER_ERROR_IF(atmospheric_steps < 2,
                     "Must have at least two atmospheric levels, got {}",
                     atmospheric_steps);
  ARTS_USER_ERROR_IF(
      temperature_perturbation_steps < 0,
      "Must have zero or positive temperature perturbations steps, got {}",
      temperature_perturbation_steps);
  ARTS_USER_ERROR_IF(
      water_vmr_perturbation_steps < 0,
      "Must have zero or positive water VMR perturbations steps, got {}",
      water_vmr_perturbation_steps);

  const Vector temperature_profile(atmospheric_steps, mean(temperature_range));
  const DescendingGrid pressure_profile =
      nlogspace(pressure_range[1], pressure_range[0], atmospheric_steps);

  const AscendingGrid water_perturbation = [vmr_value,
                                            water_vmr_range,
                                            water_vmr_perturbation_steps]() {
    Vector out = nlinspace(
        water_vmr_range[0], water_vmr_range[1], water_vmr_perturbation_steps);
    out /= vmr_value;
    return out;
  }();
  const AscendingGrid temperature_perturbation =
      nlinspace(temperature_range[0] - temperature_profile[0],
                temperature_range[1] - temperature_profile[0],
                temperature_perturbation_steps);

  const SpeciesEnumVectors vmr_profiles =
      [species = lbl::species_in_bands(absorption_bands),
       vmr_value,
       atmospheric_steps]() {
        SpeciesEnumVectors out;
        for (SpeciesEnum s : species) {
          out[s] = Vector(atmospheric_steps, vmr_value);
        }
        return out;
      }();

  absorption_lookup_tableFromProfiles(absorption_lookup_table,
                                      frequency_grid,
                                      absorption_bands,
                                      ecs_data,
                                      pressure_profile,
                                      temperature_profile,
                                      vmr_profiles,
                                      temperature_perturbation,
                                      water_perturbation,
                                      water_affected_species,
                                      isoratio_option);
}
