#include "agenda_set.h"

#include <default_gins.h>

#include "arts_options.h"

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
  static std::size_t n = 0;
  ws_pos = ws.add_wsv(WsvRecord(
      var_string("::wsv", n++).c_str(), "method value", value.type(), value));
}

AgendaCreator::AgendaCreator(Workspace& workspace, const char* name)
    : ws(workspace), agenda(ws) {
  agenda.set_name(name);
}

//! Check the agenda and move it out of here
Agenda AgendaCreator::finalize() {
  agenda.check(ws, Verbosity{});
  return std::move(agenda);
}

Agenda get_iy_main_agenda(Workspace& ws, const String& option) {
  AgendaCreator agenda(ws, "iy_main_agenda");

  using enum Options::iy_main_agendaDefaultOptions;
  switch (Options::toiy_main_agendaDefaultOptionsOrThrow(option)) {
    case Emission:
      agenda.add("ppathCalc");
      agenda.add("iyEmissionStandard");
      agenda.add("VectorSet", "geo_pos", Vector{});
      break;
    case Clearsky:
      agenda.add("ppathCalc");
      agenda.add("iyClearsky");
      agenda.add("VectorSet", "geo_pos", Vector{});
      break;
    case Transmission:
      agenda.add("Ignore", "iy_unit");
      agenda.add("Ignore", "iy_id");
      agenda.add("ppathCalc", SetWsv{"cloudbox_on", Index{0}});
      agenda.add("iyTransmissionStandard");
      agenda.add("VectorSet", "geo_pos", Vector{});
      break;
    case TransmissionUnitUnpolIntensity:
      agenda.add("Ignore", "iy_unit");
      agenda.add("Ignore", "iy_id");
      agenda.add(
          "MatrixUnitIntensity", "iy_transmitter", "stokes_dim", "f_grid");
      agenda.add("ppathCalc", SetWsv{"cloudbox_on", Index{0}});
      agenda.add("iyTransmissionStandard");
      agenda.add("VectorSet", "geo_pos", Vector{});
      break;
    case TransmissionUnitPolIntensity:
      agenda.add("Ignore", "iy_unit");
      agenda.add("Ignore", "iy_id");
      agenda.add("iy_transmitterSinglePol");
      agenda.add("ppathCalc", SetWsv{"cloudbox_on", Index{0}});
      agenda.add("iyTransmissionStandard");
      agenda.add("VectorSet", "geo_pos", Vector{});
      break;
    case Freqloop:
      agenda.add("Ignore", "diy_dx");
      agenda.add("Ignore", "iy_id");
      agenda.add("Ignore", "iy_unit");
      agenda.add("Ignore", "nlte_field");
      agenda.add("Ignore", "cloudbox_on");
      agenda.add("Ignore", "jacobian_do");
      agenda.add("iyLoopFrequencies");
      agenda.add("Touch", "ppath");
      agenda.add("VectorSet", "geo_pos", Vector{});
      break;
    case ScattMC:
      agenda.add("Ignore", "rte_pos2");
      agenda.add("Ignore", "diy_dx");
      agenda.add("Ignore", "iy_id");
      agenda.add("Ignore", "nlte_field");
      agenda.add("iyMC");
      agenda.add("Touch", "ppath");
      agenda.add("VectorSet", "geo_pos", Vector{});
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
      agenda.add("Ignore", "iy_id");
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
      agenda.add("Ignore", "rtp_pos");
      agenda.add("Ignore", "rtp_los");
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
      agenda.add("Ignore", "rte_pos2");
      agenda.add("ppathStepByStep");
      break;
    case PlaneParallel:
      agenda.add("Ignore", "ppath_lraytrace");
      agenda.add("Ignore", "rte_pos2");
      agenda.add("Ignore", "t_field");
      agenda.add("Ignore", "vmr_field");
      agenda.add("Ignore", "f_grid");
      agenda.add("ppathPlaneParallel");
      break;
    case TransmitterReceiverPath:
      agenda.add("Ignore", "cloudbox_on");
      agenda.add("Ignore", "ppath_inside_cloudbox_do");
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
      agenda.add("Ignore", "f_grid");
      agenda.add("Ignore", "ppath_lraytrace");
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
      agenda.add("Ignore", "f_grid");
      agenda.add("Ignore", "rtp_pressure");
      agenda.add("Ignore", "rtp_temperature");
      agenda.add("Ignore", "rtp_vmr");
      agenda.add("NumericSet", "refr_index_air", Numeric{1.0});
      agenda.add("NumericSet", "refr_index_air_group", Numeric{1.0});
      break;
    case GasMicrowavesEarth:
      agenda.add("Ignore", "f_grid");
      agenda.add("NumericSet", "refr_index_air", Numeric{1.0});
      agenda.add("NumericSet", "refr_index_air_group", Numeric{1.0});
      agenda.add("refr_index_airMicrowavesEarth");
      break;
    case GasInfraredEarth:
      agenda.add("Ignore", "f_grid");
      agenda.add("Ignore", "rtp_vmr");
      agenda.add("NumericSet", "refr_index_air", Numeric{1.0});
      agenda.add("NumericSet", "refr_index_air_group", Numeric{1.0});
      agenda.add("refr_index_airInfraredEarth");
      break;
    case GasMicrowavesGeneral:
      agenda.add("Ignore", "f_grid");
      agenda.add("NumericSet", "refr_index_air", Numeric{1.0});
      agenda.add("NumericSet", "refr_index_air_group", Numeric{1.0});
      agenda.add("refr_index_airMicrowavesGeneral");
      break;
    case FreeElectrons:
      agenda.add("Ignore", "rtp_pressure");
      agenda.add("Ignore", "rtp_temperature");
      agenda.add("NumericSet", "refr_index_air", Numeric{1.0});
      agenda.add("NumericSet", "refr_index_air_group", Numeric{1.0});
      agenda.add("refr_index_airFreeElectrons");
      break;
    case GasMicrowavesGeneralAndElectrons:
      agenda.add("NumericSet", "refr_index_air", Numeric{1.0});
      agenda.add("NumericSet", "refr_index_air_group", Numeric{1.0});
      agenda.add("refr_index_airMicrowavesGeneral");
      agenda.add("refr_index_airFreeElectrons");
      break;
    case GasMicrowavesEarthAndElectrons:
      agenda.add("NumericSet", "refr_index_air", Numeric{1.0});
      agenda.add("NumericSet", "refr_index_air_group", Numeric{1.0});
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
      agenda.add("Touch", "gas_scattering_coef");
      agenda.add("Touch", "gas_scattering_mat");
      agenda.add("Touch", "gas_scattering_fct_legendre");
      agenda.add("Ignore", "f_grid");
      agenda.add("Ignore", "rtp_pressure");
      agenda.add("Ignore", "rtp_temperature");
      agenda.add("Ignore", "rtp_vmr");
      agenda.add("Ignore", "gas_scattering_los_in");
      agenda.add("Ignore", "gas_scattering_los_out");
      agenda.add("Ignore", "gas_scattering_output_type");
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
      agenda.add("specular_losCalc");
      agenda.add("InterpSurfaceFieldToPosition",
                 "out=surface_skin_t",
                 "field=t_surface");
      agenda.add("surfaceFlatScalarReflectivity");
      break;
    case Specular_NoPol_ReflFix_SurfTFromt_field:
      agenda.add("specular_losCalc");
      agenda.add(
          "InterpAtmFieldToPosition", "out=surface_skin_t", "field=t_field");
      agenda.add("surfaceFlatScalarReflectivity");
      break;
    case Specular_WithPol_ReflFix_SurfTFromt_surface:
      agenda.add("specular_losCalc");
      agenda.add("InterpSurfaceFieldToPosition",
                 "out=surface_skin_t",
                 "field=t_surface");
      agenda.add("surfaceFlatReflectivity");
      break;
    case lambertian_ReflFix_SurfTFromt_surface:
      agenda.add("InterpSurfaceFieldToPosition",
                 "out=surface_skin_t",
                 "field=t_surface");
      agenda.add("specular_losCalc");
      agenda.add("surfaceLambertianSimple");
      break;
    case lambertian_ReflFix_SurfTFromt_field:
      agenda.add(
          "InterpAtmFieldToPosition", "out=surface_skin_t", "field=t_field");
      agenda.add("specular_losCalc");
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
      agenda.add("Ignore", "lon");
      agenda.add("g0Earth");
      break;
    case Io:
      agenda.add("Ignore", "lon");
      agenda.add("Ignore", "lat");
      agenda.add("g0Io");
      break;
    case Jupiter:
      agenda.add("Ignore", "lon");
      agenda.add("Ignore", "lat");
      agenda.add("g0Jupiter");
      break;
    case Mars:
      agenda.add("Ignore", "lon");
      agenda.add("Ignore", "lat");
      agenda.add("g0Mars");
      break;
    case Venus:
      agenda.add("Ignore", "lon");
      agenda.add("Ignore", "lat");
      agenda.add("g0Venus");
      break;
    case FINAL:
      break;
  }

  return agenda.finalize();
}
}  // namespace AgendaManip