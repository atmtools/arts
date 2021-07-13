#include <autoarts.h>

namespace ARTS::Agenda {
  Workspace& iy_main_agenda_emission(Workspace& ws) {
    using namespace Agenda::Method;
    using namespace Agenda::Define;
    iy_main_agenda(ws, ppathCalc(ws), iyEmissionStandard(ws));
    return ws;
  }
  
  Workspace& iy_main_agenda_transmission(Workspace& ws) {
    using namespace Agenda::Method;
    using namespace Agenda::Define;
    using namespace Var;
    iy_main_agenda(ws, Ignore(ws, iy_unit(ws)), Ignore(ws, iy_id(ws)),
                   ppathCalc(ws), iyTransmissionStandard(ws));
    return ws;
  }
  
  Workspace& iy_space_agenda_cosmic_background(Workspace& ws) {
    using namespace Agenda::Method;
    using namespace Agenda::Define;
    using namespace Var;
    iy_space_agenda(ws, Ignore(ws, rtp_pos(ws)), Ignore(ws, rtp_los(ws)),
                    MatrixCBR(ws, iy(ws), f_grid(ws)));
    return ws;
  }
  
  Workspace& iy_surface_agenda_use_surface_property(Workspace& ws) {
    using namespace Agenda::Method;
    using namespace Agenda::Define;
    iy_surface_agenda(ws, SurfaceDummy(ws), iySurfaceRtpropAgenda(ws));
    return ws;
  }
  
  Workspace& ppath_agenda_follow_sensor_los(Workspace& ws) {
    using namespace Agenda::Method;
    using namespace Agenda::Define;
    using namespace Var;
    ppath_agenda(ws, Ignore(ws, rte_pos2(ws)), ppathStepByStep(ws));
    return ws;
  }
  
  Workspace& ppath_agenda_plane_parallel(Workspace& ws) {
    using namespace Agenda::Method;
    using namespace Agenda::Define;
    using namespace Var;
    ppath_agenda(ws, Ignore(ws, ppath_lraytrace(ws)), Ignore(ws, rte_pos2(ws)),
                 Ignore(ws, t_field(ws)), Ignore(ws, vmr_field(ws)),
                 Ignore(ws, f_grid(ws)), ppathPlaneParallel(ws));
    return ws;
  }
  
  Workspace& ppath_step_agenda_geometric_path(Workspace& ws) {
    using namespace Agenda::Method;
    using namespace Agenda::Define;
    using namespace Var;
    ppath_step_agenda(ws, Ignore(ws, ppath_lraytrace(ws)), Ignore(ws, f_grid(ws)),
                      ppath_stepGeometric(ws));
    return ws;
  }
  
  Workspace& ppath_step_agenda_refracted_path(Workspace& ws) {
    using namespace Agenda::Method;
    using namespace Agenda::Define;
    ppath_step_agenda(ws, ppath_stepRefractionBasic(ws));
    return ws;
  }
  
  Workspace& propmat_clearsky_agenda_xsec_agenda(Workspace& ws) {
    using namespace Agenda::Method;
    using namespace Agenda::Define;
    using namespace Var;
    propmat_clearsky_agenda(ws, Ignore(ws, rtp_mag(ws)), Ignore(ws, rtp_los(ws)),
                            propmat_clearskyInit(ws),
                            propmat_clearskyAddXsecAgenda(ws),
                            propmat_clearskyAddLines(ws)
                           );
    return ws;
  }
  
  Workspace& propmat_clearsky_agenda_on_the_fly_zeeman(Workspace& ws) {
    using namespace Agenda::Method;
    using namespace Agenda::Define;
    propmat_clearsky_agenda(ws, propmat_clearskyInit(ws),
                            propmat_clearskyAddXsecAgenda(ws),
                            propmat_clearskyAddZeeman(ws),
                            propmat_clearskyAddLines(ws));
    return ws;
  }
  
  Workspace& abs_xsec_agenda_standard(Workspace& ws) {
    using namespace Agenda::Method;
    using namespace Agenda::Define;
    abs_xsec_agenda(ws, abs_xsec_per_speciesInit(ws),
                    abs_xsec_per_speciesAddConts(ws));
    return ws;
  }
  
  Workspace& abs_xsec_agenda_standard_with_cia(Workspace& ws) {
    using namespace Agenda::Method;
    using namespace Agenda::Define;
    abs_xsec_agenda(
      ws, abs_xsec_per_speciesInit(ws),
                    abs_xsec_per_speciesAddConts(ws), abs_xsec_per_speciesAddCIA(ws));
    return ws;
  }
  
  Workspace& surface_rtprop_agenda_blackbody_from_surface(Workspace& ws) {
    using namespace Agenda::Method;
    using namespace Agenda::Define;
    using namespace Var;
    surface_rtprop_agenda(
      ws, InterpSurfaceFieldToPosition(ws, surface_skin_t(ws), t_surface(ws)),
                          surfaceBlackbody(ws));
    return ws;
  }
  
  Workspace& surface_rtprop_agenda_blackbody_from_atmosphere(Workspace& ws) {
    using namespace Agenda::Method;
    using namespace Agenda::Define;
    using namespace Var;
    surface_rtprop_agenda(
      ws, InterpAtmFieldToPosition(ws, surface_skin_t(ws), t_field(ws)),
                          surfaceBlackbody(ws));
    return ws;
  }
  
  Workspace& geo_pos_agenda_empty(Workspace& ws) {
    using namespace Agenda::Method;
    using namespace Agenda::Define;
    using namespace Var;
    geo_pos_agenda(ws, Ignore(ws, ppath(ws)),
                   VectorSet(ws, geo_pos(ws), VectorCreate(ws, {}, "Default")));
    return ws;
  }
  
  Workspace& water_p_eq_agenda_default(Workspace& ws) {
    using namespace Agenda::Method;
    using namespace Agenda::Define;
    water_p_eq_agenda(ws, water_p_eq_fieldMK05(ws));
    return ws;
  }
}  // namespace ARTS::Agenda

namespace ARTS::Continua {
  Workspace& init(Workspace& ws) { Method::abs_cont_descriptionInit(ws); return ws; }
  
  Workspace& addH2OPWR98(Workspace& ws) { Method::abs_cont_descriptionAppend(ws, String{"H2O-PWR98"}, String{"Rosenkranz"}); return ws; }
  
  Workspace& addO2PWR98(Workspace& ws) { Method::abs_cont_descriptionAppend(ws, String{"O2-PWR98"}, String{"Rosenkranz"}); return ws; }
}  // ARTS::Continua

int main() try {
  using namespace ARTS;
  
  auto ws = init(0, 0, 0);
  
  ARTS::Agenda::iy_main_agenda_emission(ws);
  
  ARTS::Agenda::iy_space_agenda_cosmic_background(ws);
  
  ARTS::Agenda::iy_surface_agenda_use_surface_property(ws);
  
  ARTS::Agenda::ppath_agenda_follow_sensor_los(ws);
  
  ARTS::Agenda::ppath_step_agenda_geometric_path(ws);
  
  ARTS::Agenda::propmat_clearsky_agenda_xsec_agenda(ws);
  
  ARTS::Agenda::abs_xsec_agenda_standard(ws);
  
  ARTS::Agenda::surface_rtprop_agenda_blackbody_from_surface(ws);
  
  ARTS::Agenda::geo_pos_agenda_empty(ws);
  
  ARTS::Agenda::water_p_eq_agenda_default(ws);
  
  Method::jacobianOff(ws);
  Method::nlteOff(ws);
  Var::iy_unit(ws) = "PlanckBT";
  Method::Touch(ws, Var::iy_aux_vars(ws));
  Method::Touch(ws, Var::surface_props_names(ws));
  
  ARTS::Continua::init(ws);
  ARTS::Continua::addH2OPWR98(ws);
  ARTS::Continua::addO2PWR98(ws);
  
  Method::abs_speciesSet(ws, ArrayOfString{"H2O-PWR98", "O2-PWR98"});
  
  Method::isotopologue_ratiosInitFromBuiltin(ws);
  Method::VectorNLogSpace(ws, Var::p_grid(ws).value(), 51, 1e+05, 1e-4);
  
  Method::AtmosphereSet1D(ws);
  Var::lat_true(ws) = Var::lat_grid(ws);
  Var::lon_true(ws) = Var::lon_grid(ws);
  Method::Touch(ws, Var::wind_u_field(ws));
  Method::Touch(ws, Var::wind_v_field(ws));
  Method::Touch(ws, Var::wind_w_field(ws));
  Method::Touch(ws, Var::mag_u_field(ws));
  Method::Touch(ws, Var::mag_v_field(ws));
  Method::Touch(ws, Var::mag_w_field(ws));
  Method::Touch(ws, Var::nlte_field(ws));
  Method::Touch(ws, Var::surface_props_data(ws));
  
  Var::rte_alonglos_v(ws) = 0;
  Var::p_hse(ws) = 1e5;
  Var::t_field(ws) = Tensor3(51, 1, 1, 250.0);
  Var::vmr_field(ws) = Tensor4(2, 51, 1, 1, 1e-2);
  Var::z_field(ws) = Tensor3(51, 1, 1, 0);
  for (Index i=0; i<51; i++) Var::z_field(ws).value()(i, 0, 0) = 2e3 * Numeric(i);
  Method::Touch(ws, Var::nlte_field(ws));
  
  Method::refellipsoidVenus(ws, String{"Sphere"});
  Method::z_surfaceConstantAltitude(ws);
  Var::t_surface(ws) = Matrix(1, 1, 250.0);
  
  Method::Touch(ws, Var::abs_lines(ws));
  Method::abs_lines_per_speciesCreateFromLines(ws);
  
  Var::abs_f_interp_order(ws) = 1;
  Var::stokes_dim(ws) = 1;
  Var::ppath_lraytrace(ws) = 1e3;
  Var::ppath_lmax(ws) = 1e3;
  Var::rt_integration_option(ws) = "default";
  
  Method::VectorNLinSpace(ws, Var::f_grid(ws).value(), 10001, 22e9 - 500e6, 22e9 + 500e6);
  
  Var::sensor_pos(ws) = Matrix(1, 1, 100);
  Var::sensor_los(ws) = Matrix(1, 1, 75);
  Method::Touch(ws, Var::transmitter_pos(ws));
  Method::sensorOff(ws);
  Method::cloudboxOff(ws);
  
  Method::atmgeom_checkedCalc(ws);
  Method::atmfields_checkedCalc(ws);
  Method::cloudbox_checkedCalc(ws);
  Method::sensor_checkedCalc(ws);
  Method::propmat_clearsky_agenda_checkedCalc(ws);
  Method::abs_xsec_agenda_checkedCalc(ws);
  Method::lbl_checkedCalc(ws);
  
  Method::yCalc(ws);
  for (auto& n: Var::y(ws).value()) std::cout << n << ',';
  std::cout  << '\n';
  return EXIT_SUCCESS;
} catch(const std::exception& e) {
  std::ostringstream os;
  os << "EXITING WITH ERROR:\n" << e.what() << '\n';
  std::cerr << os.str();
  return EXIT_FAILURE;
}
