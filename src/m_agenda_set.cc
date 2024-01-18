#include <workspace.h>

void iy_main_agendaSet(Agenda& out, const String& option) {
  out = get_iy_main_agenda(option);
}

void iy_loop_freqs_agendaSet(Agenda& out, const String& option) {
  out = get_iy_loop_freqs_agenda(option);
}

void iy_space_agendaSet(Agenda& out, const String& option) {
  out = get_iy_space_agenda(option);
}

void iy_surface_agendaSet(Agenda& out, const String& option) {
  out = get_iy_surface_agenda(option);
}

void iy_cloudbox_agendaSet(Agenda& out, const String& option) {
  out = get_iy_cloudbox_agenda(option);
}

void refr_index_air_agendaSet(Agenda& out, const String& option) {
  out = get_refr_index_air_agenda(option);
}

void gas_scattering_agendaSet(Agenda& out, const String& option) {
  out = get_gas_scattering_agenda(option);
}

void surface_rtprop_agendaSet(Agenda& out, const String& option) {
  out = get_surface_rtprop_agenda(option);
}

void dobatch_calc_agendaSet(Agenda& out, const String& option) {
  out = get_dobatch_calc_agenda(option);
}

void ybatch_calc_agendaSet(Agenda& out, const String& option) {
  out = get_ybatch_calc_agenda(option);
}

void test_agendaSet(Agenda& out, const String& option) {
  out = get_test_agenda(option);
}

void spt_calc_agendaSet(Agenda& out, const String& option) {
  out = get_spt_calc_agenda(option);
}

void sensor_response_agendaSet(Agenda& out, const String& option) {
  out = get_sensor_response_agenda(option);
}

void propmat_clearsky_agendaSet(Agenda& out, const String& option) {
  out = get_propmat_clearsky_agenda(option);
}

void pha_mat_spt_agendaSet(Agenda& out, const String& option) {
  out = get_pha_mat_spt_agenda(option);
}

void met_profile_calc_agendaSet(Agenda& out, const String& option) {
  out = get_met_profile_calc_agenda(option);
}

void iy_radar_agendaSet(Agenda& out, const String& option) {
  out = get_iy_radar_agenda(option);
}

void iy_independent_beam_approx_agendaSet(Agenda& out, const String& option) {
  out = get_iy_independent_beam_approx_agenda(option);
}

void inversion_iterate_agendaSet(Agenda& out, const String& option) {
  out = get_inversion_iterate_agenda(option);
}

void forloop_agendaSet(Agenda& out, const String& option) {
  out = get_forloop_agenda(option);
}

void doit_scat_field_agendaSet(Agenda& out, const String& option) {
  out = get_doit_scat_field_agenda(option);
}

void doit_rte_agendaSet(Agenda& out, const String& option) {
  out = get_doit_rte_agenda(option);
}

void doit_mono_agendaSet(Agenda& out, const String& option) {
  out = get_doit_mono_agenda(option);
}

void doit_conv_test_agendaSet(Agenda& out, const String& option) {
  out = get_doit_conv_test_agenda(option);
}

void ppvar_rtprop_agendaSet(Agenda& out, const String& option) {
  out = get_ppvar_rtprop_agenda(option);
}

void rte_background_agendaSet(Agenda& out, const String& option) {
  out = get_rte_background_agenda(option);
}

void spectral_radiance_background_surface_agendaSet(Agenda& out, const String& option) {
  out = get_spectral_radiance_background_surface_agenda(option);
}

void spectral_radiance_background_space_agendaSet(Agenda& out, const String& option) {
  out = get_spectral_radiance_background_space_agenda(option);
}
