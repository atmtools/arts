#include <arts_conversions.h>
#include <debug.h>
#include <lbl_lineshape_linemixing.h>
#include <time_report.h>

#include <algorithm>
#include <cmath>

void abs_ecs_dataAddMeanAir(LinemixingEcsData& abs_ecs_data, const Vector& vmrs, const ArrayOfSpeciesEnum& specs) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(vmrs.size() != specs.size() or vmrs.empty(), "Expected equally sized, nonempty vmrs and species");
  for (Size i = 0; i < vmrs.size(); ++i) {
    ARTS_USER_ERROR_IF(not std::isfinite(vmrs[i]) or vmrs[i] < 0, "Invalid VMR for species {}: {}", specs[i], vmrs[i]);
    ARTS_USER_ERROR_IF(specs[i] == SpeciesEnum::Bath, "The bath species cannot be an input to mean-air averaging");
  }
  const Numeric total = sum(vmrs);
  ARTS_USER_ERROR_IF(not std::isfinite(total) or std::abs(total - 1) > 1e-4, "Bad vmrs [sum far from 1]: {}", vmrs);

  // Keep the legacy coefficient-averaging approximation for catalogues with
  // bath widths. Exact partner mixing requires separate diagonal widths too.
  const auto add = [](const lbl::temperature::data& source,
                      const lbl::temperature::data& accumulated,
                      Numeric                       weight,
                      bool                          first) -> lbl::temperature::data {
    auto coefficients = source.X();
    ARTS_USER_ERROR_IF(
        coefficients.empty() or not std::ranges::all_of(coefficients, [](Numeric x) { return std::isfinite(x); }),
        "Mean-air averaging requires finite, nonempty temperature-model coefficients");
    if (not first) {
      ARTS_USER_ERROR_IF(source.Type() != accumulated.Type() or coefficients.size() != accumulated.X().size(),
                         "Incompatible temperature models in mean-air averaging: {} and {}",
                         source,
                         accumulated);
    }
    coefficients *= weight;
    if (not first) coefficients += accumulated.X();
    ARTS_USER_ERROR_IF(not std::ranges::all_of(coefficients, [](Numeric x) { return std::isfinite(x); }),
                       "Non-finite mean-air temperature-model coefficients");
    return {source.Type(), coefficients};
  };

  // A failure for any isotopologue must not leave partial bath data behind.
  auto updated = abs_ecs_data;
  for (auto& [isot, data] : updated) {
    LinemixingSingleEcsData air;
    bool                    first = true;
    for (Size i = 0; i < vmrs.size(); ++i) {
      if (vmrs[i] == 0) continue;
      const auto source = data.find(specs[i]);
      ARTS_USER_ERROR_IF(source == data.end(), "Missing species {} in abs_ecs_data of isotopologue {}", specs[i], isot);
      const auto&   partner    = source->second;
      const Numeric weight     = vmrs[i] / total;
      air.scaling              = add(partner.scaling, air.scaling, weight, first);
      air.beta                 = add(partner.beta, air.beta, weight, first);
      air.lambda               = add(partner.lambda, air.lambda, weight, first);
      air.collisional_distance = add(partner.collisional_distance, air.collisional_distance, weight, first);
      first                    = false;
    }
    data[SpeciesEnum::Bath] = std::move(air);
  }
  abs_ecs_data.swap(updated);
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

  // The fitted ECS parameters are shared; each partner still has its own mass
  // and pressure-broadening coefficients in the relaxation matrix.
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

    ecs[SpeciesEnum::Nitrogen].scaling              = data(T1, {Conversion::kaycm_per_atm2hz_per_pa(0.0180), 0.85});
    ecs[SpeciesEnum::Nitrogen].lambda               = data(T1, {.81, 0.0152});
    ecs[SpeciesEnum::Nitrogen].beta                 = data(T0, {.008});
    ecs[SpeciesEnum::Nitrogen].collisional_distance = data(T0, {Conversion::angstrom2meter(2.2)});

    ecs[SpeciesEnum::Oxygen].scaling              = data(T1, {Conversion::kaycm_per_atm2hz_per_pa(0.0168), 0.5});
    ecs[SpeciesEnum::Oxygen].lambda               = data(T1, {.82, -0.091});
    ecs[SpeciesEnum::Oxygen].beta                 = data(T0, {.007});
    ecs[SpeciesEnum::Oxygen].collisional_distance = data(T0, {Conversion::angstrom2meter(2.4)});
  }
}

void abs_ecs_dataAddTran2011(LinemixingEcsData& abs_ecs_data) {
  ARTS_TIME_REPORT

  using enum LineShapeModelType;
  using data = lbl::temperature::data;

  for (const auto key : {"CO2-626"_isot, "CO2-628"_isot, "CO2-636"_isot}) {
    auto& ecs = abs_ecs_data[key];

    ecs[SpeciesEnum::CarbonDioxide].scaling              = data(T0, {Conversion::kaycm_per_atm2hz_per_pa(0.019)});
    ecs[SpeciesEnum::CarbonDioxide].lambda               = data(T0, {0.61});
    ecs[SpeciesEnum::CarbonDioxide].beta                 = data(T0, {0.052});
    ecs[SpeciesEnum::CarbonDioxide].collisional_distance = data(T0, {Conversion::angstrom2meter(5.5)});
  }
}
