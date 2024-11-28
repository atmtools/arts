#include <workspace.h>

void propagation_matrix_agendaSet(Agenda& out, const String& option) {
  out = get_propagation_matrix_agenda(option);
}

void propagation_matrix_scattering_agendaSet(Agenda& out,
                                             const String& option) {
  out = get_propagation_matrix_scattering_agenda(option);
}

void propagation_matrix_scattering_spectral_agendaSet(Agenda& out,
                                                      const String& option) {
  out = get_propagation_matrix_scattering_spectral_agenda(option);
}

void spectral_radiance_surface_agendaSet(Agenda& out, const String& option) {
  out = get_spectral_radiance_surface_agenda(option);
}

void spectral_radiance_space_agendaSet(Agenda& out, const String& option) {
  out = get_spectral_radiance_space_agenda(option);
}

void spectral_radiance_observer_agendaSet(Agenda& out, const String& option) {
  out = get_spectral_radiance_observer_agenda(option);
}

void ray_path_observer_agendaSet(Agenda& out, const String& option) {
  out = get_ray_path_observer_agenda(option);
}

void disort_settings_agendaSet(Agenda& out, const String& option) {
  out = get_disort_settings_agenda(option);
}
