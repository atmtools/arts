#include <workspace.h>

void ecs_dataAddMeanAir(LinemixingEcsData& ecs_data,
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

  for (auto& [isot, data] : ecs_data) {
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
            "Missing species {} in ecs_data of isotopologue {}", spec, isot)
      } catch (std::exception& e) {
        ARTS_USER_ERROR(
            "Error for species {}"
            " in ecs_data of isotopologue {}:\n{}",
            spec,
            isot,
            e.what())
      }
    }
  }
}

void ecs_dataInit(LinemixingEcsData& ecs_data) {
  ARTS_TIME_REPORT
  ecs_data.clear();
}

void ecs_dataAddMakarov2020(LinemixingEcsData& ecs_data) {
  ARTS_TIME_REPORT

  using enum LineShapeModelType;
  using data = lbl::temperature::data;

  auto& ecs = ecs_data["O2-66"_isot];

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

void ecs_dataAddRodrigues1997(LinemixingEcsData& ecs_data) {
  ARTS_TIME_REPORT

  using enum LineShapeModelType;
  using data = lbl::temperature::data;

  for (const auto isot : {"CO2-626"_isot, "CO2-628"_isot, "CO2-636"_isot}) {
    auto& ecs = ecs_data[isot];

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

void ecs_dataAddTran2011(LinemixingEcsData& ecs_data) {
  ARTS_TIME_REPORT

  using enum LineShapeModelType;
  using data = lbl::temperature::data;

  for (const std::string_view key : {"CO2-626", "CO2-628", "CO2-636"}) {
    auto& ecs = ecs_data[SpeciesIsotope(key)];

    ecs[SpeciesEnum::CarbonDioxide].scaling =
        data(T0, {Conversion::kaycm_per_atm2hz_per_pa(0.019)});
    ecs[SpeciesEnum::CarbonDioxide].lambda = data(T0, {0.61});
    ecs[SpeciesEnum::CarbonDioxide].beta   = data(T0, {0.052});
    ecs[SpeciesEnum::CarbonDioxide].collisional_distance =
        data(T0, {Conversion::angstrom2meter(5.5)});
  }
}
