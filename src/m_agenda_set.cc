#include <workspace.h>

void propagation_matrix_agendaSet(Agenda& out, const String& option) {
  out = get_propagation_matrix_agenda(option);
}

void spectral_radiance_background_surface_agendaSet(Agenda& out, const String& option) {
  out = get_spectral_radiance_background_surface_agenda(option);
}

void spectral_radiance_background_space_agendaSet(Agenda& out, const String& option) {
  out = get_spectral_radiance_background_space_agenda(option);
}
