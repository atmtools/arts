#include <workspace.h>

void abs_ecs_dataAddMeanAir(LinemixingEcsData& abs_ecs_data,
                            const Vector& vmrs,
                            const ArrayOfSpeciesEnum& specs) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(static_cast<Size>(vmrs.size()) != specs.size(),
                     "Mismatch dimension of vmrs and specs")
  ARTS_USER_ERROR_IF(
      std::abs(sum(vmrs) - 1.) > 1e-4, "Bad vmrs [sum far from 1]: {}", vmrs)

  auto set = [](const lbl::temperature::data& data1,
                const lbl::temperature::data& data2,
                const Numeric vmr,
                bool first) -> lbl::temperature::data {
    auto v  = data1.X();
    v      *= vmr;

    if (first) {
      return {data1.Type(), v};
    }

    v += data2.X();
    ARTS_USER_ERROR_IF(data1.Type() != data2.Type(),
                       "Type error in species temperature type for {} and {}",
                       data1.Type(),
                       data2.Type())
    return {data1.Type(), v};
  };

  for (auto& [isot, data] : abs_ecs_data) {
    auto& airdata = data[SpeciesEnum::Bath];

    bool first = true;
    for (Size i = 0; i < vmrs.size(); i++) {
      const auto spec = specs[i];
      const auto vmr  = vmrs[i];

      try {
        auto& specdata  = data.at(spec);
        airdata.scaling = set(specdata.scaling, airdata.scaling, vmr, first);
        airdata.beta    = set(specdata.beta, airdata.beta, vmr, first);
        airdata.lambda  = set(specdata.lambda, airdata.lambda, vmr, first);
        airdata.collisional_distance = set(specdata.collisional_distance,
                                           airdata.collisional_distance,
                                           vmr,
                                           first);
        first                        = false;
      } catch (std::out_of_range&) {
        ARTS_USER_ERROR(
            "Missing species {} in abs_ecs_data of isotopologue {}", spec, isot)
      } catch (std::exception& e) {
        ARTS_USER_ERROR(
            "Error for species {}"
            " in abs_ecs_data of isotopologue {}:\n{}",
            spec,
            isot,
            e.what())
      }
    }
  }
}

void abs_ecs_dataInit(LinemixingEcsData& abs_ecs_data) {
  ARTS_TIME_REPORT
  abs_ecs_data.clear();
}

void abs_ecs_dataAddMakarov2020(LinemixingEcsData& abs_ecs_data) {
  ARTS_TIME_REPORT

  using enum LineShapeModelType;
  using data = lbl::temperature::data;

  auto& ecs = abs_ecs_data["O2-66"_isot];

  // All species have the same effect, so just copy the values but change the mass (allow new mass for Air)
  auto& oxy                = ecs[SpeciesEnum::Oxygen];
  oxy.scaling              = data(T0, {1.0});
  oxy.collisional_distance = data(T0, {Conversion::angstrom2meter(0.61)});
  oxy.lambda               = data(T0, {0.39});
  oxy.beta                 = data(T0, {0.567});

  auto& nit                = ecs[SpeciesEnum::Nitrogen];
  nit.scaling              = data(T0, {1.0});
  nit.collisional_distance = data(T0, {Conversion::angstrom2meter(0.61)});
  nit.lambda               = data(T0, {0.39});
  nit.beta                 = data(T0, {0.567});
}

void abs_ecs_dataAddRodrigues1997(LinemixingEcsData& abs_ecs_data) {
  ARTS_TIME_REPORT

  using enum LineShapeModelType;
  using data = lbl::temperature::data;

  for (const auto isot : {"CO2-626"_isot, "CO2-628"_isot, "CO2-636"_isot}) {
    auto& ecs = abs_ecs_data[isot];

    ecs[SpeciesEnum::Nitrogen].scaling =
        data(T1, {Conversion::kaycm_per_atm2hz_per_pa(0.0180), 0.85});
    ecs[SpeciesEnum::Nitrogen].lambda = data(T1, {.81, 0.0152});
    ecs[SpeciesEnum::Nitrogen].beta   = data(T0, {.008});
    ecs[SpeciesEnum::Nitrogen].collisional_distance =
        data(T0, {Conversion::angstrom2meter(2.2)});

    ecs[SpeciesEnum::Oxygen].scaling =
        data(T1, {Conversion::kaycm_per_atm2hz_per_pa(0.0168), 0.5});
    ecs[SpeciesEnum::Oxygen].lambda = data(T1, {.82, -0.091});
    ecs[SpeciesEnum::Oxygen].beta   = data(T0, {.007});
    ecs[SpeciesEnum::Oxygen].collisional_distance =
        data(T0, {Conversion::angstrom2meter(2.4)});
  }
}

void abs_ecs_dataAddTran2011(LinemixingEcsData& abs_ecs_data) {
  ARTS_TIME_REPORT

  using enum LineShapeModelType;
  using data = lbl::temperature::data;

  for (const std::string_view key : {"CO2-626", "CO2-628", "CO2-636"}) {
    auto& ecs = abs_ecs_data[SpeciesIsotope(key)];

    ecs[SpeciesEnum::CarbonDioxide].scaling =
        data(T0, {Conversion::kaycm_per_atm2hz_per_pa(0.019)});
    ecs[SpeciesEnum::CarbonDioxide].lambda = data(T0, {0.61});
    ecs[SpeciesEnum::CarbonDioxide].beta   = data(T0, {0.052});
    ecs[SpeciesEnum::CarbonDioxide].collisional_distance =
        data(T0, {Conversion::angstrom2meter(5.5)});
  }
}

void abs_ecs_dataAddBoulet1999(LinemixingEcsData& abs_ecs_data) {
  ARTS_TIME_REPORT

  using enum LineShapeModelType;
  using data = lbl::temperature::data;

  // NH3-4111 (14NH3) ECS-EP parameters for H2 and He broadening
  // Based on Boulet et al. (1999) and Neshyba & Gamache (1999)
  // for NH3 in a Jovian atmosphere.
  auto& ecs = abs_ecs_data["NH3-4111"_isot];

  // H2 broadening parameters
  auto& h2                = ecs[SpeciesEnum::Hydrogen];
  h2.scaling              = data(T1, {Conversion::kaycm_per_atm2hz_per_pa(0.040), 0.73});
  h2.lambda               = data(T0, {0.65});
  h2.beta                 = data(T0, {0.006});
  h2.collisional_distance = data(T0, {Conversion::angstrom2meter(2.3)});

  // He broadening parameters
  auto& he                = ecs[SpeciesEnum::Helium];
  he.scaling              = data(T1, {Conversion::kaycm_per_atm2hz_per_pa(0.018), 0.55});
  he.lambda               = data(T0, {0.58});
  he.beta                 = data(T0, {0.003});
  he.collisional_distance = data(T0, {Conversion::angstrom2meter(1.8)});

  // NH3 self-broadening parameters (approximate, based on H2 values)
  // Needed because HITRAN catalog lines carry an NH3 self-broadening
  // model.  In a Jovian atmosphere NH3 VMR is ~300 ppm so this
  // contribution is negligible, but the entry must exist.
  auto& nh3                = ecs[SpeciesEnum::Ammonia];
  nh3.scaling              = data(T1, {Conversion::kaycm_per_atm2hz_per_pa(0.040), 0.73});
  nh3.lambda               = data(T0, {0.65});
  nh3.beta                 = data(T0, {0.006});
  nh3.collisional_distance = data(T0, {Conversion::angstrom2meter(2.3)});
}

void abs_ecs_dataAddPH3Preliminary(LinemixingEcsData& abs_ecs_data) {
  ARTS_TIME_REPORT

  using enum LineShapeModelType;
  using data = lbl::temperature::data;

  // PH3-1111 (31PH3) ECS-EP parameters — preliminary/placeholder
  // Should be refined against laboratory measurements.
  auto& ecs = abs_ecs_data["PH3-1111"_isot];

  // H2 broadening parameters (approximate, based on similarity to NH3)
  auto& h2                = ecs[SpeciesEnum::Hydrogen];
  h2.scaling              = data(T1, {Conversion::kaycm_per_atm2hz_per_pa(0.035), 0.70});
  h2.lambda               = data(T0, {0.60});
  h2.beta                 = data(T0, {0.005});
  h2.collisional_distance = data(T0, {Conversion::angstrom2meter(2.5)});

  // He broadening parameters (approximate)
  auto& he                = ecs[SpeciesEnum::Helium];
  he.scaling              = data(T1, {Conversion::kaycm_per_atm2hz_per_pa(0.015), 0.50});
  he.lambda               = data(T0, {0.55});
  he.beta                 = data(T0, {0.003});
  he.collisional_distance = data(T0, {Conversion::angstrom2meter(2.0)});

  // PH3 self-broadening parameters (approximate, based on H2 values)
  auto& ph3                = ecs[SpeciesEnum::Phosphine];
  ph3.scaling              = data(T1, {Conversion::kaycm_per_atm2hz_per_pa(0.035), 0.70});
  ph3.lambda               = data(T0, {0.60});
  ph3.beta                 = data(T0, {0.005});
  ph3.collisional_distance = data(T0, {Conversion::angstrom2meter(2.5)});
}

void abs_ecs_dataAddPieroni1999(LinemixingEcsData& abs_ecs_data) {
  ARTS_TIME_REPORT

  using enum LineShapeModelType;
  using data = lbl::temperature::data;

  // CH4-211 (12CH4) ECS-EP parameters for H2 and He broadening
  // Based on:
  //   Pieroni et al. (1999) — ECS-EP line mixing in CH4 nu3
  //   Pine & Gabard (2003)  — Q-branch line mixing data
  // Parameters adapted for a Jovian (H2/He) atmosphere.
  auto& ecs = abs_ecs_data["CH4-211"_isot];

  // H2 broadening parameters
  auto& h2                = ecs[SpeciesEnum::Hydrogen];
  h2.scaling              = data(T1, {Conversion::kaycm_per_atm2hz_per_pa(0.060), 0.75});
  h2.lambda               = data(T0, {0.70});
  h2.beta                 = data(T0, {0.008});
  h2.collisional_distance = data(T0, {Conversion::angstrom2meter(2.4)});

  // He broadening parameters
  auto& he                = ecs[SpeciesEnum::Helium];
  he.scaling              = data(T1, {Conversion::kaycm_per_atm2hz_per_pa(0.025), 0.56});
  he.lambda               = data(T0, {0.60});
  he.beta                 = data(T0, {0.004});
  he.collisional_distance = data(T0, {Conversion::angstrom2meter(1.9)});

  // CH4 self-broadening parameters (approximate, based on H2 values)
  // Needed because HITRAN catalog lines carry a CH4 self-broadening model.
  // In a Jovian atmosphere CH4 VMR is ~0.2% so this contribution is small.
  auto& ch4                = ecs[SpeciesEnum::Methane];
  ch4.scaling              = data(T1, {Conversion::kaycm_per_atm2hz_per_pa(0.060), 0.75});
  ch4.lambda               = data(T0, {0.70});
  ch4.beta                 = data(T0, {0.008});
  ch4.collisional_distance = data(T0, {Conversion::angstrom2meter(2.4)});
}
