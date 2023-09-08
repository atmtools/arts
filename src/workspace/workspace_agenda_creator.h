#pragma once

#include <optional>
#include "workspace_method_class.h"

#include "workspace_agenda_class.h"

struct SetWsv {
  std::string name;
  std::optional<std::string> other{std::nullopt};
  std::optional<Wsv> wsv{std::nullopt};

  SetWsv(std::string n);
  SetWsv(const char * const n) : SetWsv(std::string{n}) {}
  SetWsv(std::string n, WorkspaceGroup auto v) : name(std::move(n)), wsv(std::move(v)) {}
  SetWsv(std::string n, std::string v) : name(std::move(n)), other(std::move(v)) {}
};

class AgendaCreator {
  Agenda a;

public:
  AgendaCreator(std::string name) : a(std::move(name)) {}
  AgendaCreator(AgendaCreator&&) = delete;
  AgendaCreator(const AgendaCreator&) = delete;

  AgendaCreator& set(const std::string& name, WorkspaceGroup auto v);

  AgendaCreator& add(const std::string& name, std::vector<SetWsv>&& v);

  AgendaCreator& ignore(const std::string& name);
  
  template<std::convertible_to<SetWsv> ... T>
  AgendaCreator& add(const std::string& name, T&& ... v) {
    return add(name, std::vector<SetWsv>{std::forward<T>(v)...});
  }

  Agenda finalize() &&;
};

Agenda get_iy_main_agenda(const std::string& option);
Agenda get_iy_loop_freqs_agenda(const std::string& option);
Agenda get_iy_space_agenda(const std::string& option);
Agenda get_iy_surface_agenda(const std::string& option);
Agenda get_iy_cloudbox_agenda(const std::string& option);
Agenda get_ppath_agenda(const std::string& option);
Agenda get_ppath_step_agenda(const std::string& option);
Agenda get_refr_index_air_agenda(const std::string& option);
Agenda get_water_p_eq_agenda(const std::string& option);
Agenda get_gas_scattering_agenda(const std::string& option);
Agenda get_surface_rtprop_agenda(const std::string& option);
Agenda get_g0_agenda(const std::string& option);
Agenda get_dobatch_calc_agenda(const std::string& option);
Agenda get_ybatch_calc_agenda(const std::string& option);
Agenda get_test_agenda(const std::string& option);
Agenda get_spt_calc_agenda(const std::string& option);
Agenda get_sensor_response_agenda(const std::string& option);
Agenda get_propmat_clearsky_agenda(const std::string& option);
Agenda get_pha_mat_spt_agenda(const std::string& option);
Agenda get_met_profile_calc_agenda(const std::string& option);
Agenda get_main_agenda(const std::string& option);
Agenda get_jacobian_agenda(const std::string& option);
Agenda get_iy_radar_agenda(const std::string& option);
Agenda get_iy_independent_beam_approx_agenda(const std::string& option);
Agenda get_inversion_iterate_agenda(const std::string& option);
Agenda get_forloop_agenda(const std::string& option);
Agenda get_doit_scat_field_agenda(const std::string& option);
Agenda get_doit_rte_agenda(const std::string& option);
Agenda get_doit_mono_agenda(const std::string& option);
Agenda get_doit_conv_test_agenda(const std::string& option);
Agenda get_ppvar_rtprop_agenda(const std::string& option);
Agenda get_rte_background_agenda(const std::string& option);
