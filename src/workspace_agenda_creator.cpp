#include "workspace_agenda_creator.h"

#include <arts_options.h>

#include <string>
#include <unordered_map>

SetWsv::SetWsv(std::string n) : name(std::move(n)) {
  if (auto ind = name.find('='); ind not_eq name.npos) {
    *this = SetWsv(name.substr(0, ind), name.substr(ind + 1));
  }
}

AgendaCreator& AgendaCreator::set(const std::string& name,
                                  WorkspaceGroup auto v) {
  a.add(Method{name, std::move(v)});
  return *this;
}

AgendaCreator& AgendaCreator::add(const std::string& name,
                                  std::vector<SetWsv>&& v) {
  std::vector<std::string> args{};
  std::unordered_map<std::string, std::string> kwargs{};

  for (auto& wsv : v) {
    if (wsv.wsv) {
      a.add(Method{"@" + wsv.name, wsv.wsv.value()});
      kwargs[wsv.name] = "@" + wsv.name;
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

Agenda AgendaCreator::finalize() && {
  Agenda ag{std::move(a)};
  ag.finalize(true);
  return ag;
};

Agenda get_iy_main_agenda(const std::string& option) {
  AgendaCreator agenda("iy_main_agenda");

  using enum Options::iy_main_agendaDefaultOptions;
  switch (Options::toiy_main_agendaDefaultOptionsOrThrow(option)) {
    case GeometricEmission:
      agenda.add("ppathGeometric");
      agenda.add("ppathClampAltitude");
      agenda.add("ppvar_atmFromPath");
      agenda.add("ppvar_fFromPath");
      agenda.add("ppvar_rtprop_agendaExecute");
      agenda.add("background_transmittanceFromBack");
      agenda.add("rte_background_agendaExecute");
      agenda.add("spectral_radiance_pathCalcEmission");
      agenda.add("iyCopyPath");
      agenda.add("diy_dxTransform");
      agenda.add("iyUnitConversion");
      agenda.add("iy_auxFromVars");
      agenda.set("geo_pos", Vector{});
      break;
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_iy_loop_freqs_agenda(const std::string& option) {
  AgendaCreator agenda("iy_loop_freqs_agenda");

  using enum Options::iy_loop_freqs_agendaDefaultOptions;
  switch (Options::toiy_loop_freqs_agendaDefaultOptionsOrThrow(option)) {
    case Emission:
      agenda.add("ppathCalc");
      agenda.add("iyEmissionStandard");
      break;
    case Transmission:
      agenda.add("ppathCalc");
      agenda.add("iyTransmissionStandard");
      break;
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_iy_space_agenda(const std::string& option) {
  AgendaCreator agenda("iy_space_agenda");

  using enum Options::iy_space_agendaDefaultOptions;
  switch (Options::toiy_space_agendaDefaultOptionsOrThrow(option)) {
    case CosmicBackground:
      agenda.add("MatrixCBR", "iy", "f_grid");
      break;
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_iy_surface_agenda(const std::string& option) {
  AgendaCreator agenda("iy_surface_agenda");

  using enum Options::iy_surface_agendaDefaultOptions;
  switch (Options::toiy_surface_agendaDefaultOptionsOrThrow(option)) {
    case UseSurfaceRtprop:
      agenda.add("SurfaceDummy");
      agenda.add("iySurfaceRtpropAgenda");
      break;
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_iy_cloudbox_agenda(const std::string& option) {
  AgendaCreator agenda("iy_cloudbox_agenda");

  using enum Options::iy_cloudbox_agendaDefaultOptions;
  switch (Options::toiy_cloudbox_agendaDefaultOptionsOrThrow(option)) {
    case LinInterpField:
      agenda.add("iyInterpCloudboxField");
      break;
    case QuarticInterpField:
      agenda.add("iyInterpCloudboxField", SetWsv{"za_interp_order", Index{4}});
      break;
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_refr_index_air_agenda(const std::string& option) {
  AgendaCreator agenda("refr_index_air_agenda");

  using enum Options::refr_index_air_agendaDefaultOptions;
  switch (Options::torefr_index_air_agendaDefaultOptionsOrThrow(option)) {
    case NoRefrac:
      agenda.set("refr_index_air", Numeric{1.0});
      agenda.set("refr_index_air_group", Numeric{1.0});
      break;
    case GasMicrowavesEarth:
      agenda.set("refr_index_air", Numeric{1.0});
      agenda.set("refr_index_air_group", Numeric{1.0});
      agenda.add("refr_index_airMicrowavesEarth");
      break;
    case GasInfraredEarth:
      agenda.set("refr_index_air", Numeric{1.0});
      agenda.set("refr_index_air_group", Numeric{1.0});
      agenda.add("refr_index_airInfraredEarth");
      break;
    case GasMicrowavesGeneral:
      agenda.set("refr_index_air", Numeric{1.0});
      agenda.set("refr_index_air_group", Numeric{1.0});
      agenda.add("refr_index_airMicrowavesGeneral");
      break;
    case FreeElectrons:
      agenda.set("refr_index_air", Numeric{1.0});
      agenda.set("refr_index_air_group", Numeric{1.0});
      agenda.add("refr_index_airFreeElectrons");
      break;
    case GasMicrowavesGeneralAndElectrons:
      agenda.set("refr_index_air", Numeric{1.0});
      agenda.set("refr_index_air_group", Numeric{1.0});
      agenda.add("refr_index_airMicrowavesGeneral");
      agenda.add("refr_index_airFreeElectrons");
      break;
    case GasMicrowavesEarthAndElectrons:
      agenda.set("refr_index_air", Numeric{1.0});
      agenda.set("refr_index_air_group", Numeric{1.0});
      agenda.add("refr_index_airMicrowavesEarth");
      agenda.add("refr_index_airFreeElectrons");
      break;
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_gas_scattering_agenda(const std::string& option) {
  AgendaCreator agenda("gas_scattering_agenda");

  using enum Options::gas_scattering_agendaDefaultOptions;
  switch (Options::togas_scattering_agendaDefaultOptionsOrThrow(option)) {
    case Dummy:
      break;
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_surface_rtprop_agenda(const std::string& option) {
  AgendaCreator agenda("surface_rtprop_agenda");

  using enum Options::surface_rtprop_agendaDefaultOptions;
  switch (Options::tosurface_rtprop_agendaDefaultOptionsOrThrow(option)) {
    case Blackbody_SurfTFromt_surface:
      agenda.add("InterpSurfaceFieldToPosition");
      agenda.add("surfaceBlackbody");
      break;
    case Blackbody_SurfTFromt_field:
      agenda.add("surface_pointFromAtm");
      agenda.add("surfaceBlackbody");
      break;
    case Specular_NoPol_ReflFix_SurfTFromt_surface:
      agenda.add("specular_losCalc");
      agenda.add("InterpSurfaceFieldToPosition");
      agenda.add("surfaceFlatScalarReflectivity");
      break;
    case Specular_NoPol_ReflFix_SurfTFromt_field:
      agenda.add("specular_losCalc");
      agenda.add("surface_pointFromAtm");
      agenda.add("surfaceFlatScalarReflectivity");
      break;
    case Specular_WithPol_ReflFix_SurfTFromt_surface:
      agenda.add("specular_losCalc");
      agenda.add("InterpSurfaceFieldToPosition");
      agenda.add("surfaceFlatReflectivity");
      break;
    case lambertian_ReflFix_SurfTFromt_surface:
      agenda.add("InterpSurfaceFieldToPosition");
      agenda.add("specular_losCalc");
      agenda.add("surfaceLambertianSimple");
      break;
    case lambertian_ReflFix_SurfTFromt_field:
      agenda.add("surface_pointFromAtm");
      agenda.add("specular_losCalc");
      agenda.add("surfaceLambertianSimple");
      break;
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_dobatch_calc_agenda(const std::string& option) {
  AgendaCreator agenda("dobatch_calc_agenda");

  using enum Options::dobatch_calc_agendaDefaultOptions;
  switch (Options::todobatch_calc_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_ybatch_calc_agenda(const std::string& option) {
  AgendaCreator agenda("ybatch_calc_agenda");

  using enum Options::ybatch_calc_agendaDefaultOptions;
  switch (Options::toybatch_calc_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_test_agenda(const std::string& option) {
  AgendaCreator agenda("test_agenda");

  using enum Options::test_agendaDefaultOptions;
  switch (Options::totest_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_spt_calc_agenda(const std::string& option) {
  AgendaCreator agenda("spt_calc_agenda");

  using enum Options::spt_calc_agendaDefaultOptions;
  switch (Options::tospt_calc_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_sensor_response_agenda(const std::string& option) {
  AgendaCreator agenda("sensor_response_agenda");

  using enum Options::sensor_response_agendaDefaultOptions;
  switch (Options::tosensor_response_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_propmat_clearsky_agenda(const std::string& option) {
  AgendaCreator agenda("propmat_clearsky_agenda");

  using enum Options::propmat_clearsky_agendaDefaultOptions;
  switch (Options::topropmat_clearsky_agendaDefaultOptionsOrThrow(option)) {
    case Empty:
      agenda.add("propmat_clearskyInit");
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_pha_mat_spt_agenda(const std::string& option) {
  AgendaCreator agenda("pha_mat_spt_agenda");

  using enum Options::pha_mat_spt_agendaDefaultOptions;
  switch (Options::topha_mat_spt_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_met_profile_calc_agenda(const std::string& option) {
  AgendaCreator agenda("met_profile_calc_agenda");

  using enum Options::met_profile_calc_agendaDefaultOptions;
  switch (Options::tomet_profile_calc_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_iy_radar_agenda(const std::string& option) {
  AgendaCreator agenda("iy_radar_agenda");

  using enum Options::iy_radar_agendaDefaultOptions;
  switch (Options::toiy_radar_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_iy_independent_beam_approx_agenda(const std::string& option) {
  AgendaCreator agenda("iy_independent_beam_approx_agenda");

  using enum Options::iy_independent_beam_approx_agendaDefaultOptions;
  switch (Options::toiy_independent_beam_approx_agendaDefaultOptionsOrThrow(
      option)) {
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_inversion_iterate_agenda(const std::string& option) {
  AgendaCreator agenda("inversion_iterate_agenda");

  using enum Options::inversion_iterate_agendaDefaultOptions;
  switch (Options::toinversion_iterate_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_forloop_agenda(const std::string& option) {
  AgendaCreator agenda("forloop_agenda");

  using enum Options::forloop_agendaDefaultOptions;
  switch (Options::toforloop_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_doit_scat_field_agenda(const std::string& option) {
  AgendaCreator agenda("doit_scat_field_agenda");

  using enum Options::doit_scat_field_agendaDefaultOptions;
  switch (Options::todoit_scat_field_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_doit_rte_agenda(const std::string& option) {
  AgendaCreator agenda("doit_rte_agenda");

  using enum Options::doit_rte_agendaDefaultOptions;
  switch (Options::todoit_rte_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_doit_mono_agenda(const std::string& option) {
  AgendaCreator agenda("doit_mono_agenda");

  using enum Options::doit_mono_agendaDefaultOptions;
  switch (Options::todoit_mono_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_doit_conv_test_agenda(const std::string& option) {
  AgendaCreator agenda("doit_conv_test_agenda");

  using enum Options::doit_conv_test_agendaDefaultOptions;
  switch (Options::todoit_conv_test_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_ppvar_rtprop_agenda(const std::string& option) {
  AgendaCreator agenda("ppvar_rtprop_agenda");

  using enum Options::ppvar_rtprop_agendaDefaultOptions;
  switch (Options::toppvar_rtprop_agendaDefaultOptionsOrThrow(option)) {
    case Propmat:
      agenda.add("ppvar_propmatCalc");
      agenda.add("spectral_radiance_path_sourceFromPropmat");
      agenda.add("ppvar_tramatCalc");
      agenda.add("ppvar_cumtramatForward");
      break;
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_rte_background_agenda(const std::string& option) {
  AgendaCreator agenda("rte_background_agenda");

  using enum Options::rte_background_agendaDefaultOptions;
  switch (Options::torte_background_agendaDefaultOptionsOrThrow(option)) {
    case ByPath:
      agenda.add("iyBackground");
      agenda.add("background_radFromMatrix", "iy_mat=iy");
      break;
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_spectral_radiance_background_space_agenda(const std::string& option) {
  AgendaCreator agenda("spectral_radiance_background_surface_agenda");

  using enum Options::spectral_radiance_background_space_agendaDefaultOptions;
  switch (Options::tospectral_radiance_background_space_agendaDefaultOptionsOrThrow(option)) {
    case UniformCosmicBackground:
      agenda.add("spectral_radiance_backgroundUniformCosmicBackground");
      break;
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}

Agenda get_spectral_radiance_background_surface_agenda(const std::string& option) {
  AgendaCreator agenda("spectral_radiance_background_surface_agenda");

  using enum Options::spectral_radiance_background_surface_agendaDefaultOptions;
  switch (Options::tospectral_radiance_background_surface_agendaDefaultOptionsOrThrow(option)) {
    case Blackbody:
      agenda.add("spectral_radiance_backgroundSurfaceBlackbody");
      break;
    case FINAL:
      break;
  }

  return std::move(agenda).finalize();
}
