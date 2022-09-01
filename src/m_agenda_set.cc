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
