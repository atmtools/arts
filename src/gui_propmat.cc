#include <atomic>
#include <chrono>
#include <cstdlib>
#include <functional>
#include <future>
#include <iterator>
#include <mutex>
#include <stdexcept>
#include <thread>

#include "agenda_class.h"
#include "artstime.h"
#include "auto_md.h"
#include "debug.h"
#include "energylevelmap.h"
#include "jacobian.h"
#include "messages.h"
#include "propagationmatrix.h"
#include "species_tags.h"
#include "transmissionmatrix.h"
#include "workspace_ng.h"

#ifdef ARTS_GUI_ENABLED
#include <gui/propmat.h>

namespace PropmatClearskyAgendaGUI {
void compute(Workspace& ws,
             ARTSGUI::PropmatClearsky::ComputeValues& v,
             const Agenda& propmat_clearsky_agenda) {
  propmat_clearsky_agendaExecute(ws,
                                 v.pm,
                                 v.sv,
                                 v.aopm,
                                 v.aosv,
                                 v.jacobian_quantities,
                                 v.select_abs_species,
                                 v.f_grid,
                                 v.rtp_mag,
                                 v.rtp_los,
                                 v.rtp_pressure,
                                 v.rtp_temperature,
                                 v.rtp_nlte,
                                 v.rtp_vmr,
                                 propmat_clearsky_agenda);
}

bool run(ARTSGUI::PropmatClearsky::ResultsArray& ret,
         ARTSGUI::PropmatClearsky::Control& ctrl,
         Workspace&& ws,
         const Agenda&& propmat_clearsky_agenda,
         ArrayOfRetrievalQuantity& jacobian_quantities,
         ArrayOfSpeciesTag& select_abs_species,
         Vector& f_grid,
         Vector& rtp_mag,
         Vector& rtp_los,
         Numeric& rtp_pressure,
         Numeric& rtp_temperature,
         EnergyLevelMap& rtp_nlte,
         Vector& rtp_vmr,
         Numeric& transmission_distance) {
  for (auto& v : ret) v.ok.store(false);
  ARTSGUI::PropmatClearsky::ComputeValues v;
  bool error = false;

  while (not error) {
    std::this_thread::sleep_for(10ms);

    try {
      if (ctrl.exit.load()) return true;

      if (ctrl.run.load()) {
        std::lock_guard allow_copy{ctrl.copy};

        v.jacobian_quantities = jacobian_quantities;
        v.select_abs_species = select_abs_species;
        v.f_grid = f_grid;
        v.rtp_mag = rtp_mag;
        v.rtp_los = rtp_los;
        v.rtp_pressure = rtp_pressure;
        v.rtp_temperature = rtp_temperature;
        v.rtp_nlte = rtp_nlte;
        v.rtp_vmr = rtp_vmr;
        v.transmission_distance = transmission_distance;
      } else {
        continue;
      }

      compute(ws, v, propmat_clearsky_agenda);

      v.tm = TransmissionMatrix(v.pm.NumberOfFrequencies(), v.pm.StokesDimensions());
      v.aotm = ArrayOfTransmissionMatrix(v.jacobian_quantities.nelem(), v.tm);
      auto local_aotm = v.aotm;
      stepwise_transmission(v.tm,
                            v.aotm,
                            local_aotm,
                            v.pm,
                            v.pm,
                            v.aopm,
                            v.aopm,
                            v.transmission_distance,
                            0,
                            0,
                            -1);

      // Lock after the compute to copy values
      std::lock_guard allow_copy{ctrl.copy};

      // Copy over values
      ret.at(ctrl.pos).value = v;
      ret.at(ctrl.pos).ok.store(true);

      ctrl.run.store(false);
    } catch (std::runtime_error& e) {
      ctrl.error = e.what();
      ctrl.exit.store(true);
      error = true;
    }
  }

  return false;
}
}  // namespace PropmatClearskyAgendaGUI
#endif  // ARTS_GUI_ENABLED

void propmat_clearsky_agendaGUI(Workspace& ws [[maybe_unused]],
                                const Agenda& propmat_clearsky_agenda [[maybe_unused]],
                                const ArrayOfArrayOfSpeciesTag& abs_species [[maybe_unused]],
                                const Index& load [[maybe_unused]],
                                const Verbosity& verbosity [[maybe_unused]]) {
#ifdef ARTS_GUI_ENABLED
  ARTSGUI::PropmatClearsky::ResultsArray res;
  ARTSGUI::PropmatClearsky::Control ctrl;

  // Initialize values to something
  ArrayOfRetrievalQuantity jacobian_quantities{};
  ArrayOfSpeciesTag select_abs_species{};
  Vector f_grid(1e9, 1000, 1e9);
  Vector rtp_mag(3, 0);
  Vector rtp_los(2, 0);
  Numeric rtp_pressure{1000};
  Numeric rtp_temperature(300);
  EnergyLevelMap rtp_nlte{};
  Vector rtp_vmr(abs_species.nelem(), 1.0/Numeric(abs_species.nelem()));
  Numeric transmission_distance{1'000};

  // Set some defaults
  if (load) {
    if (auto* val = ws.get<Vector>("f_grid")) f_grid = *val;
    if (auto* val = ws.get<Vector>("rtp_mag")) rtp_mag = *val;
    if (auto* val = ws.get<Vector>("rtp_los")) rtp_los = *val;
    if (auto* val = ws.get<Vector>("rtp_vmr")) rtp_vmr = *val;
    if (auto* val = ws.get<Numeric>("rtp_pressure")) rtp_pressure = *val;
    if (auto* val = ws.get<Numeric>("rtp_temperature")) rtp_temperature = *val;
    if (auto* val = ws.get<EnergyLevelMap>("rtp_nlte")) rtp_nlte = *val;
  }

  auto success = std::async(std::launch::async,
                            &PropmatClearskyAgendaGUI::run,
                            std::ref(res),
                            std::ref(ctrl),
                            ws,
                            propmat_clearsky_agenda,
                            std::ref(jacobian_quantities),
                            std::ref(select_abs_species),
                            std::ref(f_grid),
                            std::ref(rtp_mag),
                            std::ref(rtp_los),
                            std::ref(rtp_pressure),
                            std::ref(rtp_temperature),
                            std::ref(rtp_nlte),
                            std::ref(rtp_vmr),
                            std::ref(transmission_distance));

  if (std::getenv("ARTS_HEADLESS")) {
    CREATE_OUT1;
    out1 << "Omitting GUI because ARTS_HEADLESS is set.\n";

    ctrl.run.store(true);
    while (not(res[0].ok.load() or ctrl.exit.load())) {
      std::this_thread::sleep_for(10ms);
    }
    ctrl.exit.store(true);
  } else {
    ARTSGUI::propmat(res,
                     ctrl,
                     jacobian_quantities,
                     select_abs_species,
                     f_grid,
                     rtp_mag,
                     rtp_los,
                     rtp_pressure,
                     rtp_temperature,
                     rtp_nlte,
                     rtp_vmr,
                     transmission_distance,
                     ArrayOfArrayOfSpeciesTag{abs_species});
  }

  bool invalid_state = not success.get();
  ARTS_USER_ERROR_IF(invalid_state, '\n', ctrl.error)

#else   // ARTS_GUI_ENABLED not defined
  ARTS_USER_ERROR("Did not compile with GUI")
#endif  // ARTS_GUI_ENABLED
}
