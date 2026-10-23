#include <arts_constants.h>
#include <physics_funcs.h>
#include <workspace.h>
#include <workspace_agenda_creator.h>

namespace {
Numeric water_mass_density(const AtmPoint& atm_point) {
  constexpr SpeciesEnum water = "H2O"_spec;

  if (not atm_point.has(water))
    arts::user_error(
        "Cannot derive water mass density: the atmospheric point has no H2O VMR. "
        "Set H2O or provide water_mass_density explicitly.");

  return atm_point.number_density(water) * atm_point.mean_mass(water) * Constant::m_u;
}

constexpr matpack::cdata_t<SpeciesEnum, 6> microwave_refractive_species = {
    "N2"_spec, "O2"_spec, "CO2"_spec, "H2"_spec, "He"_spec, "H2O"_spec};

constexpr matpack::cdata_t<Numeric, 6> microwave_reference_refractivities = {
    293.81e-6, 266.95e-6, 495.16e-6, 135.77e-6, 34.51e-6, 5338.89e-6};
}  // namespace

void single_dispersionAddGasMicrowavesEarth(Numeric&        single_dispersion,
                                            const AtmPoint& atm_point,
                                            const Numeric&  k1,
                                            const Numeric&  k2,
                                            const Numeric&  k3) try {
  ARTS_TIME_REPORT

  if (atm_point.temperature <= 0.0)
    arts::user_error("Microwave gas refractivity requires a positive atmospheric temperature, got {} K",
                     atm_point.temperature);
  if (atm_point.pressure < 0.0)
    arts::user_error("Microwave gas refractivity requires a non-negative atmospheric pressure, got {} Pa",
                     atm_point.pressure);

  constexpr SpeciesEnum water          = "H2O"_spec;
  const Numeric         water_pressure = atm_point.has(water) ? atm_point.pressure * atm_point[water] : 0.0;

  single_dispersion +=
      (k1 * (atm_point.pressure - water_pressure) + (k2 + k3 / atm_point.temperature) * water_pressure) /
      atm_point.temperature;
}
ARTS_METHOD_ERROR_CATCH

void single_dispersionAddGasMicrowavesGeneral(Numeric& single_dispersion, const AtmPoint& atm_point) try {
  ARTS_TIME_REPORT

  if (atm_point.temperature <= 0.0)
    arts::user_error("Microwave gas refractivity requires a positive atmospheric temperature, got {} K",
                     atm_point.temperature);
  if (atm_point.pressure < 0.0)
    arts::user_error("Microwave gas refractivity requires a non-negative atmospheric pressure, got {} Pa",
                     atm_point.pressure);

  Numeric refractivity = 0.0;
  Numeric vmr_sum      = 0.0;
  for (Size i = 0; i < microwave_refractive_species.size(); ++i) {
    const auto species = microwave_refractive_species[i];
    if (not atm_point.has(species)) continue;

    const Numeric vmr  = atm_point[species];
    vmr_sum           += vmr;
    refractivity      += microwave_reference_refractivities[i] * vmr;
  }

  if (vmr_sum == 0.0) return;

  constexpr Numeric reference_pressure     = Conversion::torr2pa(760.0);
  constexpr Numeric reference_temperature  = 273.15;
  single_dispersion                       += refractivity / vmr_sum * (atm_point.pressure / reference_pressure) *
                                             (reference_temperature / atm_point.temperature);
}
ARTS_METHOD_ERROR_CATCH

void single_dispersionAddWaterVisibleNIRHarvey98(Numeric&        single_dispersion,
                                                 const Numeric&  freq,
                                                 const AtmPoint& atm_point,
                                                 const Numeric&  water_mass_density_,
                                                 const Index&    check_validity) try {
  ARTS_TIME_REPORT

  const Numeric density = water_mass_density_ < 0.0 ? water_mass_density(atm_point) : water_mass_density_;
  single_dispersion +=
      refractive_index_water_visible_nir_harvey98(freq, atm_point.temperature, density, check_validity != 0) - 1.0;
}
ARTS_METHOD_ERROR_CATCH

void single_propmat_agendaSetWaterVisibleNIRHarvey98(Agenda&        single_propmat_agenda,
                                                     const Numeric& water_mass_density_,
                                                     const Index&   check_validity) {
  ARTS_TIME_REPORT

  AgendaCreator creator("single_propmat_agenda");
  creator.add("single_propmatInit");
  creator.add("single_dispersionAddWaterVisibleNIRHarvey98",
              SetWsv{"water_mass_density", water_mass_density_},
              SetWsv{"check_validity", check_validity});
  single_propmat_agenda = std::move(creator).finalize(true);
}

void single_propmat_agendaSetGasMicrowavesEarth(Agenda&        single_propmat_agenda,
                                                const Numeric& k1,
                                                const Numeric& k2,
                                                const Numeric& k3) {
  ARTS_TIME_REPORT

  AgendaCreator creator("single_propmat_agenda");
  creator.add("single_propmatInit");
  creator.add("single_dispersionAddGasMicrowavesEarth", SetWsv{"k1", k1}, SetWsv{"k2", k2}, SetWsv{"k3", k3});
  single_propmat_agenda = std::move(creator).finalize(true);
}

void single_propmat_agendaSetGasMicrowavesGeneral(Agenda& single_propmat_agenda) {
  ARTS_TIME_REPORT

  AgendaCreator creator("single_propmat_agenda");
  creator.add("single_propmatInit");
  creator.add("single_dispersionAddGasMicrowavesGeneral");
  single_propmat_agenda = std::move(creator).finalize(true);
}
