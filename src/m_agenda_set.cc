#include <default_gins.h>

#include "agenda_class.h"
#include "agenda_set.h"
#include "global_data.h"

using namespace AgendaManip;

void iy_main_agendaSet(Workspace& ws,
                       Agenda& out,
                       const String& option,
                       const Verbosity&) {
  out = get_iy_main_agenda(ws, option);
}

void iy_main_agendaSetByPart(Workspace &ws, Agenda &out,
                             const String &rte_option,
                             const String &propagation_properties_option,
                             const String &background_option,
                             const String &ppath_option, const Verbosity &) {
  out = get_iy_main_agenda(ws, rte_option, propagation_properties_option,
                           background_option, ppath_option);
}

void iy_loop_freqs_agendaSet(Workspace& ws,
                             Agenda& out,
                             const String& option,
                             const Verbosity&) {
  out = get_iy_loop_freqs_agenda(ws, option);
}

void iy_space_agendaSet(Workspace& ws,
                        Agenda& out,
                        const String& option,
                        const Verbosity&) {
  out = get_iy_space_agenda(ws, option);
}

void iy_surface_agendaSet(Workspace& ws,
                          Agenda& out,
                          const String& option,
                          const Verbosity&) {
  out = get_iy_surface_agenda(ws, option);
}

void iy_cloudbox_agendaSet(Workspace& ws,
                           Agenda& out,
                           const String& option,
                           const Verbosity&) {
  out = get_iy_cloudbox_agenda(ws, option);
}

void ppath_agendaSet(Workspace& ws,
                     Agenda& out,
                     const String& option,
                     const Verbosity&) {
  out = get_ppath_agenda(ws, option);
}

void ppath_step_agendaSet(Workspace& ws,
                          Agenda& out,
                          const String& option,
                          const Verbosity&) {
  out = get_ppath_step_agenda(ws, option);
}

void refr_index_air_agendaSet(Workspace& ws,
                              Agenda& out,
                              const String& option,
                              const Verbosity&) {
  out = get_refr_index_air_agenda(ws, option);
}

void water_p_eq_agendaSet(Workspace& ws,
                          Agenda& out,
                          const String& option,
                          const Verbosity&) {
  out = get_water_p_eq_agenda(ws, option);
}

void gas_scattering_agendaSet(Workspace& ws,
                              Agenda& out,
                              const String& option,
                              const Verbosity&) {
  out = get_gas_scattering_agenda(ws, option);
}

void surface_rtprop_agendaSet(Workspace& ws,
                              Agenda& out,
                              const String& option,
                              const Verbosity&) {
  out = get_surface_rtprop_agenda(ws, option);
}

void g0_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_g0_agenda(ws, option);
}

void dobatch_calc_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_dobatch_calc_agenda(ws, option);
}

void ybatch_calc_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_ybatch_calc_agenda(ws, option);
}

void test_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_test_agenda(ws, option);
}

void spt_calc_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_spt_calc_agenda(ws, option);
}

void sensor_response_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_sensor_response_agenda(ws, option);
}

void propmat_clearsky_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_propmat_clearsky_agenda(ws, option);
}

void pha_mat_spt_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_pha_mat_spt_agenda(ws, option);
}

void met_profile_calc_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_met_profile_calc_agenda(ws, option);
}

void main_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_main_agenda(ws, option);
}

void jacobian_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_jacobian_agenda(ws, option);
}

void iy_radar_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_iy_radar_agenda(ws, option);
}

void iy_independent_beam_approx_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_iy_independent_beam_approx_agenda(ws, option);
}

void inversion_iterate_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_inversion_iterate_agenda(ws, option);
}

void forloop_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_forloop_agenda(ws, option);
}

void doit_scat_field_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_doit_scat_field_agenda(ws, option);
}

void doit_rte_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_doit_rte_agenda(ws, option);
}

void doit_mono_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_doit_mono_agenda(ws, option);
}

void doit_conv_test_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_doit_conv_test_agenda(ws, option);
}

void ppvar_rtprop_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_ppvar_rtprop_agenda(ws, option);
}

void rte_background_agendaSet(Workspace& ws,
                  Agenda& out,
                  const String& option,
                  const Verbosity&) {
  out = get_rte_background_agenda(ws, option);
}
