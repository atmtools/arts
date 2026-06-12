#include "workspace_agenda_creator.h"

#include <auto_wsa.h>
#include <auto_wsa_options.h>
#include <workspace.h>

#include <string>
#include <unordered_map>

SetWsv::SetWsv(std::string n) : name(std::move(n)) {
  if (auto ind = name.find('='); ind not_eq name.npos) {
    *this = SetWsv(name.substr(0, ind), name.substr(ind + 1));
  }
}

AgendaCreator& AgendaCreator::add(const std::string& name,
                                  std::vector<SetWsv>&& v) {
  std::vector<std::string> args{};
  std::unordered_map<std::string, std::string> kwargs{};

  for (auto& wsv : v) {
    if (wsv.wsv) {
      a.add(Method{named_input_prefix + wsv.name, wsv.wsv.value()});
      kwargs[wsv.name] = named_input_prefix + wsv.name;
    } else if (wsv.other) {
      kwargs[wsv.name] = wsv.other.value();
    } else {
      args.push_back(wsv.name);
    }
  }

  a.add(Method{name, args, kwargs});

  return *this;
}

AgendaCreator& AgendaCreator::ignore(const std::string& name) {
  a.add(Method{"Ignore", {name}, {}});

  return *this;
}

Agenda AgendaCreator::finalize(bool fix) && {
  Agenda ag{std::move(a)};
  ag.finalize(fix);
  return ag;
};

Agenda get_spectral_propmat_scat_agenda(const std::string_view option) {
  AgendaCreator agenda("spectral_propmat_scat_agenda");

  using enum spectral_propmat_scat_agendaPredefined;
  switch (to<spectral_propmat_scat_agendaPredefined>(option)) {
    case AirSimple:
      agenda.add("spectral_propmat_scatInit");
      agenda.add("spectral_propmat_scatAirSimple");
  }

  return std::move(agenda).finalize(false);
}

Agenda get_spectral_propmat_scat_spectral_agenda(
    const std::string_view option) {
  AgendaCreator agenda("spectral_propmat_scat_spectral_agenda");

  using enum spectral_propmat_scat_spectral_agendaPredefined;
  switch (to<spectral_propmat_scat_spectral_agendaPredefined>(option)) {
    case FromSpeciesTRO:
      agenda.add("spectral_propmat_scatSpectralInit");
      agenda.add("spectral_propmat_scatAddSpectralScatteringSpeciesTRO");
  }

  return std::move(agenda).finalize(false);
}

Agenda get_spectral_propmat_agenda(const std::string_view option) {
  AgendaCreator agenda("spectral_propmat_agenda");

  using enum spectral_propmat_agendaPredefined;
  switch (to<spectral_propmat_agendaPredefined>(option)) {
    case Empty: agenda.add("spectral_propmatInit");
  }

  return std::move(agenda).finalize(true);
}

Agenda get_spectral_rad_observer_agenda(const std::string_view option) {
  AgendaCreator agenda("spectral_rad_observer_agenda");

  using enum spectral_rad_observer_agendaPredefined;
  switch (to<spectral_rad_observer_agendaPredefined>(option)) {
    case Emission:
      agenda.add("ray_path_observer_agendaExecute");
      agenda.add("spectral_radClearskyEmission");
      agenda.add("spectral_rad_jacAddSensorJacobianPerturbations");
      break;
    case EmissionAdaptiveHalfsteps:
      agenda.add("ray_path_observer_agendaExecute");
      agenda.add("ray_pointBackground");
      agenda.add("spectral_rad_bkgAgendasAtEndOfPath");
      agenda.add("atm_pathFromPath");
      agenda.add("freq_grid_pathFromPath");
      agenda.add("spectral_propmat_pathFromPath");
      agenda.add("spectral_propmat_pathAddAdaptiveHalfPath");
      agenda.add("spectral_radSetToBackground");
      agenda.add("spectral_radSinglePathEmissionFrequencyLoop");
      agenda.add("spectral_rad_jacAddSensorJacobianPerturbations");
      break;
    case EmissionNoSensor:
      agenda.add("ray_path_observer_agendaExecute");
      agenda.add("spectral_radClearskyEmission");
      break;
  }

  return std::move(agenda).finalize(false);
}

Agenda get_spectral_rad_space_agenda(const std::string_view option) {
  AgendaCreator agenda("spectral_rad_space_agenda");

  using enum spectral_rad_space_agendaPredefined;
  switch (to<spectral_rad_space_agendaPredefined>(option)) {
    case UniformCosmicBackground:
      agenda.add("spectral_radUniformCosmicBackground");
      agenda.add("spectral_rad_jacEmpty");
      break;
    case SunOrCosmicBackground:
      agenda.add("spectral_radSunsOrCosmicBackground");
      agenda.add("spectral_rad_jacEmpty");
      break;
    case Transmission: agenda.add("spectral_radDefaultTransmission"); break;
  }

  return std::move(agenda).finalize(true);
}

Agenda get_spectral_rad_surface_agenda(const std::string_view option) {
  AgendaCreator agenda("spectral_rad_surface_agenda");

  using enum spectral_rad_surface_agendaPredefined;
  switch (to<spectral_rad_surface_agendaPredefined>(option)) {
    case Blackbody:    agenda.add("spectral_radSurfaceBlackbody"); break;
    case Transmission: agenda.add("spectral_radDefaultTransmission"); break;
    case SurfaceReflectance:
      agenda.add("spectral_radSurfaceReflectance");
      break;
  }

  return std::move(agenda).finalize(true);
}

Agenda get_inversion_iterate_agenda(const std::string_view option) {
  AgendaCreator agenda("inversion_iterate_agenda");

  using enum inversion_iterate_agendaPredefined;
  switch (to<inversion_iterate_agendaPredefined>(option)) {
    case Full:
      agenda.add("UpdateModelStates");
      agenda.add("measurement_inversion_agendaExecute");
      break;
  }

  return std::move(agenda).finalize(true);
}

Agenda get_measurement_inversion_agenda(const std::string_view option) {
  AgendaCreator agenda("measurement_inversion_agenda");

  using enum measurement_inversion_agendaPredefined;
  switch (to<measurement_inversion_agendaPredefined>(option)) {
    case LowMemory:
      agenda.add("measurement_vec_errorFromModelState");
      agenda.add("jac_targetsConditionalClear");
      agenda.add("measurement_vecFromSensor", SetWsv("kernel", "Low Memory"s));
      agenda.add("measurement_jacTransformations");
      agenda.add("measurement_vecConditionalAddError");
      agenda.add("measurement_vec_fitFromMeasurement");
      break;
    case HighPerformance:
      agenda.add("measurement_vec_errorFromModelState");
      agenda.add("jac_targetsConditionalClear");
      agenda.add("measurement_vecFromSensor",
                 SetWsv("kernel", "High Performance"s));
      agenda.add("measurement_jacTransformations");
      agenda.add("measurement_vecConditionalAddError");
      agenda.add("measurement_vec_fitFromMeasurement");
      break;
  }

  return std::move(agenda).finalize(true);
}

Agenda get_spectral_surf_refl_agenda(const std::string_view option) {
  AgendaCreator agenda("spectral_surf_refl_agenda");

  using enum spectral_surf_refl_agendaPredefined;
  switch (to<spectral_surf_refl_agendaPredefined>(option)) {
    case FlatScalar: agenda.add("spectral_surf_reflFlatScalar"); break;
    case FlatRealFresnel:
      agenda.add("spectral_surf_reflFlatRealFresnel");
      break;
  }

  return std::move(agenda).finalize(true);
}

Agenda get_disort_settings_downwelling_wrapper_agenda(
    const std::string_view option) {
  AgendaCreator agenda("disort_settings_downwelling_wrapper_agenda");

  using enum disort_settings_downwelling_wrapper_agendaPredefined;
  switch (to<disort_settings_downwelling_wrapper_agendaPredefined>(option)) {
    case Standard:
      agenda.add("disort_settings_agendaExecute");
      agenda.add("disort_settingsDownwellingObserver");
      break;
  }

  return std::move(agenda).finalize(true);
}

Agenda get_single_rad_space_agenda(const std::string_view option) {
  AgendaCreator agenda("single_rad_space_agenda");

  using enum single_rad_space_agendaPredefined;
  switch (to<single_rad_space_agendaPredefined>(option)) {
    case WrapGrid:
      agenda.add("freq_gridFromSingleFrequency");
      agenda.add("spectral_rad_space_agendaExecute");
      agenda.add("single_radFromVector");
      break;
  }

  return std::move(agenda).finalize(true);
}

Agenda get_single_rad_surface_agenda(const std::string_view option) {
  AgendaCreator agenda("single_rad_surface_agenda");

  using enum single_rad_surface_agendaPredefined;
  switch (to<single_rad_surface_agendaPredefined>(option)) {
    case WrapGrid:
      agenda.add("freq_gridFromSingleFrequency");
      agenda.add("spectral_rad_surface_agendaExecute");
      agenda.add("single_radFromVector");
      break;
  }

  return std::move(agenda).finalize(true);
}

Agenda get_ray_point_back_propagation_agenda(const std::string_view option) {
  AgendaCreator agenda("ray_point_back_propagation_agenda");

  using enum ray_point_back_propagation_agendaPredefined;
  switch (to<ray_point_back_propagation_agendaPredefined>(option)) {
    case GeometricStepwise:  agenda.add("ray_pointPastGeometric"); break;
    case RefractiveStepwise: agenda.add("ray_pointPastRefractive"); break;
  }

  return std::move(agenda).finalize(true);
}

Agenda get_spectral_propmat_and_atm_path_agenda(const std::string_view option) {
  AgendaCreator agenda("spectral_propmat_and_atm_path_agenda");

  using enum spectral_propmat_and_atm_path_agendaPredefined;
  switch (to<spectral_propmat_and_atm_path_agendaPredefined>(option)) {
    case Default:
      agenda.add("atm_pathFromPath");
      agenda.add("freq_grid_pathFromPath");
      agenda.add("spectral_propmat_pathFromPath");
      agenda.add("Touch", "ray_path");
      break;
    case AdaptiveHalfPath:
      agenda.add("atm_pathFromPath");
      agenda.add("freq_grid_pathFromPath");
      agenda.add("spectral_propmat_pathFromPath");
      agenda.add("spectral_propmat_pathAddAdaptiveHalfPath");
      break;
    case Profile2Path:
      agenda.add("Profile2Path");
      agenda.add("Touch", "ray_path");
      break;
  }

  return std::move(agenda).finalize(true);
}

Agenda get_ray_path_observer_agenda(const std::string_view option) {
  AgendaCreator agenda("ray_path_observer_agenda");

  using enum ray_path_observer_agendaPredefined;
  switch (to<ray_path_observer_agendaPredefined>(option)) {
    case GeometricProfile:
      agenda.add("ray_pathInit",
                 SetWsv{"pos", "obs_pos"},
                 SetWsv{"los", "obs_los"},
                 SetWsv("as_sensor", Index{1}));
      agenda.add("ray_pathSetGeometricExtremes",
                 SetWsv("surf_search_accuracy", Numeric{0.1}),
                 SetWsv("surf_safe_search", Index{1}));
      agenda.add("ray_pathAddGeometricAltitudeGridCrossings");
      break;
    case GeometricDefault: {
      Agenda out;
      auto& wsm_def = internal_workspace_methods()
                          .at("ray_path_observer_agendaSetGeometric")
                          .gin_value;
      ray_path_observer_agendaSetGeometric(out,
                                           wsm_def[0]->get<String>(),
                                           wsm_def[1]->get<Numeric>(),
                                           wsm_def[2]->get<Numeric>(),
                                           wsm_def[3]->get<AtmKey>(),
                                           wsm_def[4]->get<Index>(),
                                           wsm_def[5]->get<Index>(),
                                           wsm_def[6]->get<Index>(),
                                           wsm_def[7]->get<Index>(),
                                           wsm_def[8]->get<Index>(),
                                           wsm_def[9]->get<Index>(),
                                           wsm_def[10]->get<Index>());
      return out;
    } break;
  }

  return std::move(agenda).finalize(true);
}
