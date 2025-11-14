#include "workspace_agenda_creator.h"

#include <auto_wsa.h>
#include <auto_wsa_options.h>

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

Agenda get_propagation_matrix_scattering_agenda(const std::string_view option) {
  AgendaCreator agenda("propagation_matrix_scattering_agenda");

  using enum propagation_matrix_scattering_agendaPredefined;
  switch (to<propagation_matrix_scattering_agendaPredefined>(option)) {
    case AirSimple:
      agenda.add("propagation_matrix_scatteringInit");
      agenda.add("propagation_matrix_scatteringAirSimple");
  }

  return std::move(agenda).finalize(false);
}

Agenda get_propagation_matrix_scattering_spectral_agenda(
    const std::string_view option) {
  AgendaCreator agenda("propagation_matrix_scattering_spectral_agenda");

  using enum propagation_matrix_scattering_spectral_agendaPredefined;
  switch (to<propagation_matrix_scattering_spectral_agendaPredefined>(option)) {
    case FromSpeciesTRO:
      agenda.add("propagation_matrix_scatteringSpectralInit");
      agenda.add(
          "propagation_matrix_scatteringAddSpectralScatteringSpeciesTRO");
  }

  return std::move(agenda).finalize(false);
}

Agenda get_propagation_matrix_agenda(const std::string_view option) {
  AgendaCreator agenda("propagation_matrix_agenda");

  using enum propagation_matrix_agendaPredefined;
  switch (to<propagation_matrix_agendaPredefined>(option)) {
    case Empty: agenda.add("propagation_matrixInit");
  }

  return std::move(agenda).finalize(true);
}

Agenda get_spectral_radiance_observer_agenda(const std::string_view option) {
  AgendaCreator agenda("spectral_radiance_observer_agenda");

  using enum spectral_radiance_observer_agendaPredefined;
  switch (to<spectral_radiance_observer_agendaPredefined>(option)) {
    case Emission:
      agenda.add("ray_path_observer_agendaExecute");
      agenda.add("spectral_radianceClearskyEmission");
      agenda.add("spectral_radiance_jacobianAddSensorJacobianPerturbations");
      break;
    case EmissionNoSensor:
      agenda.add("ray_path_observer_agendaExecute");
      agenda.add("spectral_radianceClearskyEmission");
      break;
  }

  return std::move(agenda).finalize(false);
}

Agenda get_spectral_radiance_space_agenda(const std::string_view option) {
  AgendaCreator agenda("spectral_radiance_space_agenda");

  using enum spectral_radiance_space_agendaPredefined;
  switch (to<spectral_radiance_space_agendaPredefined>(option)) {
    case UniformCosmicBackground:
      agenda.add("spectral_radianceUniformCosmicBackground");
      agenda.add("spectral_radiance_jacobianEmpty");
      break;
    case SunOrCosmicBackground:
      agenda.add("spectral_radianceSunsOrCosmicBackground");
      agenda.add("spectral_radiance_jacobianEmpty");
      break;
    case Transmission:
      agenda.add("spectral_radianceDefaultTransmission");
      break;
  }

  return std::move(agenda).finalize(true);
}

Agenda get_spectral_radiance_surface_agenda(const std::string_view option) {
  AgendaCreator agenda("spectral_radiance_surface_agenda");

  using enum spectral_radiance_surface_agendaPredefined;
  switch (to<spectral_radiance_surface_agendaPredefined>(option)) {
    case Blackbody: agenda.add("spectral_radianceSurfaceBlackbody"); break;
    case Transmission:
      agenda.add("spectral_radianceDefaultTransmission");
      break;
    case SurfaceReflectance:
      agenda.add("spectral_radianceSurfaceReflectance");
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
    case Standard:
      agenda.add("measurement_vector_errorFromModelState");
      agenda.add("jacobian_targetsConditionalClear");
      agenda.add("measurement_vectorFromSensor");
      agenda.add("measurement_jacobianTransformations");
      agenda.add("measurement_vectorConditionalAddError");
      agenda.add("measurement_vector_fittedFromMeasurement");
      break;
  }

  return std::move(agenda).finalize(true);
}

Agenda get_surface_reflectance_agenda(const std::string_view option) {
  AgendaCreator agenda("surface_reflectance_agenda");

  using enum surface_reflectance_agendaPredefined;
  switch (to<surface_reflectance_agendaPredefined>(option)) {
    case FlatScalar: agenda.add("surface_reflectanceFlatScalar"); break;
    case FlatRealFresnel:
      agenda.add("surface_reflectanceFlatRealFresnel");
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
      agenda.add("spectral_radiance_space_agendaExecute");
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
      agenda.add("spectral_radiance_surface_agendaExecute");
      agenda.add("single_radFromVector");
      break;
  }

  return std::move(agenda).finalize(true);
}

Agenda get_ray_path_point_back_propagation_agenda(
    const std::string_view option) {
  AgendaCreator agenda("ray_path_point_back_propagation_agenda");

  using enum ray_path_point_back_propagation_agendaPredefined;
  switch (to<ray_path_point_back_propagation_agendaPredefined>(option)) {
    case GeometricStepwise:  agenda.add("ray_path_pointPastGeometric"); break;
    case RefractiveStepwise: agenda.add("ray_path_pointPastRefractive"); break;
  }

  return std::move(agenda).finalize(true);
}
