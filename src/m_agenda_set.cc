#include <workspace.h>

void propagation_matrix_agendaSet(Agenda& out, const String& option) {
  out = get_propagation_matrix_agenda(option);
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

void propagation_path_observer_agendaSet(Agenda& out, const String& option) {
  out = get_propagation_path_observer_agenda(option);
}
