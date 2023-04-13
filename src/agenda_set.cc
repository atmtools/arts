#include "agenda_set.h"

#include <default_gins.h>

#include "arts_options.h"
#include "transmissionmatrix.h"

namespace AgendaManip {
std::ostream& operator<<(std::ostream& os, const AgendaMethodVariable& x) {
  return os << x.name << "@" << x.ws_pos;
}

Array<AgendaMethodVariable> sorted_mdrecord(Workspace& ws,
                                            const String& method_name) {
  static const Index any_pos = global_data::WsvGroupMap.at("Any");

  auto& method =
      global_data::md_data_raw.at(global_data::MdRawMap.at(method_name));

  Array<AgendaMethodVariable> var_order;

  for (auto& output : method.Out()) {
    auto& var = var_order.emplace_back();
    var.name = ws.wsv_data_ptr->operator[](output).Name();
    var.group = ws.wsv_data_ptr->operator[](output).Group();
    var.ws_pos = output;
    var.input = false;
  }

  for (Index i = 0; i < method.GOut().nelem(); i++) {
    auto& var = var_order.emplace_back();
    var.name = method.GOut()[i];
    var.group = method.GOutType()[i];
    var.g_pos = i;
    var.input = false;
    var.any = var.group == any_pos;
  }

  for (auto& input : method.In()) {
    // Skip inout
    if (std::any_of(method.Out().begin(),
                    method.Out().end(),
                    [input](auto& output) { return input == output; }))
      continue;

    auto& var = var_order.emplace_back();
    var.name = ws.wsv_data_ptr->operator[](input).Name();
    var.group = ws.wsv_data_ptr->operator[](input).Group();
    var.ws_pos = input;
    var.input = true;
  }

  for (Index i = 0; i < method.GIn().nelem(); i++) {
    auto& var = var_order.emplace_back();
    var.name = method.GIn()[i];
    var.group = method.GInType()[i];
    var.g_pos = i;
    var.input = true;
    var.any = var.group == any_pos;

    // Create defaults (on the workspace)
    if (method.GInDefault()[i] not_eq NODEF) {
      var.ws_pos =
          create_workspace_gin_default_internal(ws, method_name, var.name);
    }
  }

  return var_order;
}

std::pair<ArrayOfIndex, ArrayOfIndex> split_io(
    Array<AgendaMethodVariable>& var_order) {
  ArrayOfIndex in, out;
  for (auto& var : var_order) {
    if (var.input)
      in.push_back(var.ws_pos);
    else
      out.push_back(var.ws_pos);
  }
  return {in, out};
}

SetWsv::SetWsv(const std::string& x, const std::string& y)
    : test(opt::NameOnly), str(x + "=" + y) {}

SetWsv::SetWsv(std::string_view x) : test(opt::NameOnly), str(x) {}

MethodVariable::MethodVariable(Workspace& ws,
                               const Array<AgendaMethodVariable>& list,
                               const SetWsv& wsv) {
  using enum SetWsv::opt;
  switch (wsv.test) {
    case NameOnly: {
      const std::string_view expr = wsv.str;
      auto equal_sign = expr.find('=');
      if (equal_sign == expr.npos) {
        positional = true;
        ws_pos = ws.WsvMap_ptr->at(expr);
      } else {
        positional = false;
        auto rhs = expr;
        rhs.remove_prefix(equal_sign);

        auto lhs = expr;
        lhs.remove_suffix(rhs.size());

        while (lhs.size() and (std::isspace(lhs.front()) or lhs.front() == '='))
          lhs.remove_prefix(1);
        while (lhs.size() and (std::isspace(lhs.back()) or lhs.back() == '='))
          lhs.remove_suffix(1);
        while (rhs.size() and (std::isspace(rhs.front()) or rhs.front() == '='))
          rhs.remove_prefix(1);
        while (rhs.size() and (std::isspace(rhs.back()) or rhs.back() == '='))
          rhs.remove_suffix(1);

        method_position(list, lhs);
        ws_pos = ws.WsvMap_ptr->at(rhs);
      }
    } break;
    case ValueOnly:
      positional = true;
      wsv_position(ws, wsv.val);
      break;
    case NameAndValue:
      positional = false;
      method_position(list, wsv.str);
      wsv_position(ws, wsv.val);
      break;
  }
}

void MethodVariable::add_del(Workspace& ws, Agenda& a) const {
  if (new_value)
    a.push_back(MRecord(
        global_data::MdMap.at(var_string(
            "Delete_sg_",
            global_data::wsv_groups.at(ws.wsv_data_ptr->at(ws_pos).Group()))),
        {},
        {ws_pos},
        {},
        Agenda{ws}));
}

void MethodVariable::add_set(Workspace& ws, Agenda& a) const {
  if (new_value)
    a.push_back(MRecord(
        global_data::MdMap.at(var_string(
            global_data::wsv_groups.at(ws.wsv_data_ptr->at(ws_pos).Group())
                .name,
            "Set")),
        {ws_pos},
        {},
        ws.wsv_data_ptr->at(ws_pos).default_value(),
        Agenda{ws}));
}

void MethodVariable::method_position(const Array<AgendaMethodVariable>& list,
                                     const std::string_view rhs) {
  auto ptr = std::find_if(
      list.begin(), list.end(), [rhs](auto& x) { return x.name == rhs; });
  ARTS_ASSERT(
      ptr not_eq list.end(), "Wrongly named parameter selection: ", rhs);
  method_pos = std::distance(list.begin(), ptr);
}

void MethodVariable::wsv_position(Workspace& ws, const TokVal& value) {
  new_value = true;
  ws_pos = ws.add_wsv(WsvRecord(
      var_string("::wsv", ws.nelem()).c_str(), "method value", value.type(), value));
}

AgendaCreator::AgendaCreator(Workspace& workspace, const char* name)
    : ws(workspace), agenda(ws) {
  agenda.set_name(name);
}

//! Check the agenda and move it out of here
Agenda AgendaCreator::finalize() {
  auto& agrecord = global_data::agenda_data.at(global_data::AgendaMap.at(agenda.name()));

  ArrayOfIndex input, output;
  for (auto& method: agenda.Methods()) {
    for (auto& wsv: method.In()) input.push_back(wsv);
    for (auto& wsv: method.Out()) output.push_back(wsv);
  }
  std::sort(input.begin(), input.end());
  std::sort(output.begin(), output.end());

  //! Add Ignore(WSV)
  for (auto& wsv: agrecord.In()) {
    auto test = [wsv](auto& x){return x == wsv;};
    if (std::none_of(input.begin(), input.end(), test) and std::none_of(output.begin(), output.end(), test)) {
      auto& group = global_data::wsv_groups.at(ws.wsv_data_ptr -> at(wsv).Group()).name;
      auto pos = global_data::MdMap.at(var_string("Ignore_sg_", group));
      agenda.push_back(MRecord(pos, {}, {wsv}, {}, Agenda(ws)));
    }
  }

  //! Add Touch(WSV)
  for (auto& wsv: agrecord.Out()) {
    auto test = [wsv](auto& x){return x == wsv;};
    if (std::none_of(output.begin(), output.end(), test)) {
      auto& group = global_data::wsv_groups.at(ws.wsv_data_ptr -> at(wsv).Group()).name;
      auto pos = global_data::MdMap.at(var_string("Touch_sg_", group));
      agenda.push_back(MRecord(pos, {wsv}, {}, {}, Agenda(ws)));
    }
  }

  //! Finally check that the agenda is OK (it should be a developer error if this fails!)
  agenda.check(ws, Verbosity{});

  return agenda;
}

void AgendaCreator::set(const std::string_view var, const TokVal& value) {
  auto pos = global_data::MdMap.at(var_string(value.type(), "Set"));
  auto ws_pos = ws.WsvMap_ptr -> at(var);
  if (value.holdsAgenda())
    agenda.push_back(MRecord(pos, {ws_pos}, {}, {}, Agenda{value}));
  else
    agenda.push_back(MRecord(pos, {ws_pos}, {}, value, Agenda{ws}));
}

void AgendaCreator::ignore(const std::string_view var) { add("Ignore", var); }

Agenda get_iy_main_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "iy_main_agenda");

  using enum Options::iy_main_agendaDefaultOptions;
  switch (Options::toiy_main_agendaDefaultOptionsOrThrow(option)) {
    case Emission:
      agenda.add("ppathCalc");
      agenda.add("iyEmissionStandard");
      agenda.set("geo_pos", Vector{});
      break;
    case EmissionPlaneParallel:
      agenda.add("ppathPlaneParallel");
      agenda.add("iyEmissionStandard");
      agenda.set("geo_pos", Vector{});
      break;
    case Clearsky:
      agenda.add("ppathCalc");
      agenda.add("iyClearsky");
      agenda.set("geo_pos", Vector{});
      break;
    case Transmission:
      agenda.add("ppathCalc", SetWsv{"cloudbox_on", Index{0}});
      agenda.add("iyTransmissionStandard");
      agenda.set("geo_pos", Vector{});
      break;
    case TransmissionUnitUnpolIntensity:
      agenda.add(
          "MatrixUnitIntensity", "iy_transmitter", "stokes_dim", "f_grid");
      agenda.add("ppathCalc", SetWsv{"cloudbox_on", Index{0}});
      agenda.add("iyTransmissionStandard");
      agenda.set("geo_pos", Vector{});
      break;
    case TransmissionUnitPolIntensity:
      agenda.add("iy_transmitterSinglePol");
      agenda.add("ppathCalc", SetWsv{"cloudbox_on", Index{0}});
      agenda.add("iyTransmissionStandard");
      agenda.set("geo_pos", Vector{});
      break;
    case Freqloop:
      agenda.add("iyLoopFrequencies");
      agenda.set("geo_pos", Vector{});
      agenda.ignore("diy_dx");
      break;
    case ScattMC:
      agenda.add("iyMC");
      agenda.set("geo_pos", Vector{});
      agenda.ignore("diy_dx");
      break;
    case FINAL:
      break;
  }

  return agenda.finalize();
}

Agenda get_iy_loop_freqs_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "iy_loop_freqs_agenda");

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

  return agenda.finalize();
}

Agenda get_iy_space_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "iy_space_agenda");

  using enum Options::iy_space_agendaDefaultOptions;
  switch (Options::toiy_space_agendaDefaultOptionsOrThrow(option)) {
    case CosmicBackground:
      agenda.add("MatrixCBR", "iy", "stokes_dim", "f_grid");
      break;
    case FINAL:
      break;
  }

  return agenda.finalize();
}

Agenda get_iy_surface_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "iy_surface_agenda");

  using enum Options::iy_surface_agendaDefaultOptions;
  switch (Options::toiy_surface_agendaDefaultOptionsOrThrow(option)) {
    case UseSurfaceRtprop:
      agenda.add("SurfaceDummy");
      agenda.add("iySurfaceRtpropAgenda");
      break;
    case FINAL:
      break;
  }

  return agenda.finalize();
}

Agenda get_iy_cloudbox_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "iy_cloudbox_agenda");

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

  return agenda.finalize();
}

Agenda get_ppath_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "ppath_agenda");

  using enum Options::ppath_agendaDefaultOptions;
  switch (Options::toppath_agendaDefaultOptionsOrThrow(option)) {
    case FollowSensorLosPath:
      agenda.add("ppathStepByStep");
      break;
    case PlaneParallel:
      agenda.add("ppathPlaneParallel");
      break;
    case TransmitterReceiverPath:
      agenda.add("rte_losGeometricFromRtePosToRtePos2");
      agenda.add("ppathFromRtePos2");
      break;
    case FINAL:
      break;
  }

  return agenda.finalize();
}

Agenda get_ppath_step_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "ppath_step_agenda");

  using enum Options::ppath_step_agendaDefaultOptions;
  switch (Options::toppath_step_agendaDefaultOptionsOrThrow(option)) {
    case GeometricPath:
      agenda.add("ppath_stepGeometric");
      break;
    case RefractedPath:
      agenda.add("ppath_stepRefractionBasic");
      break;
    case FINAL:
      break;
  }

  return agenda.finalize();
}

Agenda get_refr_index_air_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "refr_index_air_agenda");

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

  return agenda.finalize();
}

Agenda get_water_p_eq_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "water_p_eq_agenda");

  using enum Options::water_p_eq_agendaDefaultOptions;
  switch (Options::towater_p_eq_agendaDefaultOptionsOrThrow(option)) {
    case MK05:
      agenda.add("water_p_eq_fieldMK05");
      break;
    case FINAL:
      break;
  }

  return agenda.finalize();
}

Agenda get_gas_scattering_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "gas_scattering_agenda");

  using enum Options::gas_scattering_agendaDefaultOptions;
  switch (Options::togas_scattering_agendaDefaultOptionsOrThrow(option)) {
    case Dummy:
      break;
    case FINAL:
      break;
  }

  return agenda.finalize();
}

Agenda get_surface_rtprop_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "surface_rtprop_agenda");

  using enum Options::surface_rtprop_agendaDefaultOptions;
  switch (Options::tosurface_rtprop_agendaDefaultOptionsOrThrow(option)) {
    case Blackbody_SurfTFromt_surface:
      agenda.add("InterpSurfaceFieldToPosition",
                 "out=surface_skin_t",
                 "field=t_surface");
      agenda.add("surfaceBlackbody");
      break;
    case Blackbody_SurfTFromt_field:
      agenda.add(
          "InterpAtmFieldToPosition", "out=surface_skin_t", "field=t_field");
      agenda.add("surfaceBlackbody");
      break;
    case Specular_NoPol_ReflFix_SurfTFromt_surface:
      agenda.add("specular_losCalcOld");
      agenda.add("InterpSurfaceFieldToPosition",
                 "out=surface_skin_t",
                 "field=t_surface");
      agenda.add("surfaceFlatScalarReflectivity");
      break;
    case Specular_NoPol_ReflFix_SurfTFromt_field:
      agenda.add("specular_losCalcOld");
      agenda.add(
          "InterpAtmFieldToPosition", "out=surface_skin_t", "field=t_field");
      agenda.add("surfaceFlatScalarReflectivity");
      break;
    case Specular_WithPol_ReflFix_SurfTFromt_surface:
      agenda.add("specular_losCalcOld");
      agenda.add("InterpSurfaceFieldToPosition",
                 "out=surface_skin_t",
                 "field=t_surface");
      agenda.add("surfaceFlatReflectivity");
      break;
    case lambertian_ReflFix_SurfTFromt_surface:
      agenda.add("InterpSurfaceFieldToPosition",
                 "out=surface_skin_t",
                 "field=t_surface");
      agenda.add("specular_losCalcOld");
      agenda.add("surfaceLambertianSimple");
      break;
    case lambertian_ReflFix_SurfTFromt_field:
      agenda.add(
          "InterpAtmFieldToPosition", "out=surface_skin_t", "field=t_field");
      agenda.add("specular_losCalcOld");
      agenda.add("surfaceLambertianSimple");
      break;
    case FINAL:
      break;
  }

  return agenda.finalize();
}

Agenda get_g0_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "g0_agenda");

  using enum Options::g0_agendaDefaultOptions;
  switch (Options::tog0_agendaDefaultOptionsOrThrow(option)) {
    case Earth:
      agenda.add("g0Earth");
      break;
    case Io:
      agenda.add("g0Io");
      break;
    case Jupiter:
      agenda.add("g0Jupiter");
      break;
    case Mars:
      agenda.add("g0Mars");
      break;
    case Venus:
      agenda.add("g0Venus");
      break;
    case FINAL:
      break;
  }

  return agenda.finalize();
}

Agenda get_dobatch_calc_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "dobatch_calc_agenda");

  using enum Options::dobatch_calc_agendaDefaultOptions;
  switch (Options::todobatch_calc_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }
  
  return agenda.finalize();
}

Agenda get_ybatch_calc_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "ybatch_calc_agenda");

  using enum Options::ybatch_calc_agendaDefaultOptions;
  switch (Options::toybatch_calc_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }
  
  return agenda.finalize();
}

Agenda get_test_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "test_agenda");

  using enum Options::test_agendaDefaultOptions;
  switch (Options::totest_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }
  
  return agenda.finalize();
}

Agenda get_spt_calc_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "spt_calc_agenda");

  using enum Options::spt_calc_agendaDefaultOptions;
  switch (Options::tospt_calc_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }
  
  return agenda.finalize();
}

Agenda get_sensor_response_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "sensor_response_agenda");

  using enum Options::sensor_response_agendaDefaultOptions;
  switch (Options::tosensor_response_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }
  
  return agenda.finalize();
}

Agenda get_propmat_clearsky_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "propmat_clearsky_agenda");

  using enum Options::propmat_clearsky_agendaDefaultOptions;
  switch (Options::topropmat_clearsky_agendaDefaultOptionsOrThrow(option)) {
    case Empty:
      agenda.add("propmat_clearskyInit");
    case FINAL:
      break;
  }
  
  return agenda.finalize();
}

Agenda get_pha_mat_spt_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "pha_mat_spt_agenda");

  using enum Options::pha_mat_spt_agendaDefaultOptions;
  switch (Options::topha_mat_spt_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }
  
  return agenda.finalize();
}

Agenda get_met_profile_calc_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "met_profile_calc_agenda");

  using enum Options::met_profile_calc_agendaDefaultOptions;
  switch (Options::tomet_profile_calc_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }
  
  return agenda.finalize();
}

Agenda get_main_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "main_agenda");

  using enum Options::main_agendaDefaultOptions;
  switch (Options::tomain_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }
  
  return agenda.finalize();
}

Agenda get_jacobian_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "jacobian_agenda");

  using enum Options::jacobian_agendaDefaultOptions;
  switch (Options::tojacobian_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }
  
  return agenda.finalize();
}

Agenda get_iy_radar_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "iy_radar_agenda");

  using enum Options::iy_radar_agendaDefaultOptions;
  switch (Options::toiy_radar_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }
  
  return agenda.finalize();
}

Agenda get_iy_independent_beam_approx_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "iy_independent_beam_approx_agenda");

  using enum Options::iy_independent_beam_approx_agendaDefaultOptions;
  switch (Options::toiy_independent_beam_approx_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }
  
  return agenda.finalize();
}

Agenda get_inversion_iterate_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "inversion_iterate_agenda");

  using enum Options::inversion_iterate_agendaDefaultOptions;
  switch (Options::toinversion_iterate_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }
  
  return agenda.finalize();
}

Agenda get_forloop_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "forloop_agenda");

  using enum Options::forloop_agendaDefaultOptions;
  switch (Options::toforloop_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }
  
  return agenda.finalize();
}

Agenda get_doit_scat_field_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "doit_scat_field_agenda");

  using enum Options::doit_scat_field_agendaDefaultOptions;
  switch (Options::todoit_scat_field_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }
  
  return agenda.finalize();
}

Agenda get_doit_rte_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "doit_rte_agenda");

  using enum Options::doit_rte_agendaDefaultOptions;
  switch (Options::todoit_rte_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }
  
  return agenda.finalize();
}

Agenda get_doit_mono_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "doit_mono_agenda");

  using enum Options::doit_mono_agendaDefaultOptions;
  switch (Options::todoit_mono_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }
  
  return agenda.finalize();
}

Agenda get_doit_conv_test_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "doit_conv_test_agenda");

  using enum Options::doit_conv_test_agendaDefaultOptions;
  switch (Options::todoit_conv_test_agendaDefaultOptionsOrThrow(option)) {
    case FINAL:
      break;
  }
  
  return agenda.finalize();
}

Agenda get_ppvar_rtprop_agenda(Workspace &ws, const String &option) {
  AgendaCreator agenda(ws, "ppvar_rtprop_agenda");

  using enum Options::ppvar_rtprop_agendaDefaultOptions;
  switch (Options::toppvar_rtprop_agendaDefaultOptionsOrThrow(option)) {
  case Propmat:
      agenda.add("ppvar_propmatCalc");
      agenda.add("ppvar_srcFromPropmat");
      break;
  case FINAL:
      break;
  }

  return agenda.finalize();
}
}  // namespace AgendaManip
