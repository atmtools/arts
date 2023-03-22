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
#include "atm.h"
#include "auto_md.h"
#include "debug.h"
#include "energylevelmap.h"
#include "jacobian.h"
#include "matpack_math.h"
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
                                 v.rtp_los,
                                 v.atm_point,
                                 propmat_clearsky_agenda);
}

bool run(ARTSGUI::PropmatClearsky::ResultsArray& ret,
         ARTSGUI::PropmatClearsky::Control& ctrl,
         Workspace& ws,
         const Agenda& propmat_clearsky_agenda,
         ArrayOfRetrievalQuantity& jacobian_quantities,
         ArrayOfSpeciesTag& select_abs_species,
         Vector& f_grid,
         Vector& rtp_los,
         AtmPoint& atm_point,
         Numeric& transmission_distance) {
  for (auto& v : ret) v.ok.store(false);
  ARTSGUI::PropmatClearsky::ComputeValues v;

  while (true) {
    std::this_thread::sleep_for(10ms);

    try {
      if (ctrl.exit.load()) return true;

      if (ctrl.error.load()) continue;

      if (ctrl.run.load()) {
        std::lock_guard allow_copy{ctrl.copy};

        v.jacobian_quantities = jacobian_quantities;
        v.select_abs_species = select_abs_species;
        v.f_grid = f_grid;
        v.rtp_los = rtp_los;
        v.atm_point = atm_point;
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
      ctrl.errmsg = e.what();
      ctrl.error.store(true);
      ctrl.run.store(false);
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
  Vector f_grid=uniform_grid(1e9, 1000, 1e9);
  Vector rtp_los(2, 0);
  Numeric transmission_distance{1'000};
  AtmPoint atm_point{
    Atm::Key::temperature, 300,
    Atm::Key::pressure, 1000
  };
  for (auto& spec: abs_species) atm_point[spec] = 1.0 / static_cast<Numeric>(abs_species.nelem());

  // Set some defaults
  if (load) {
    if (auto* val = ws.get<Vector>("f_grid")) f_grid = *val;
    if (auto* val = ws.get<Vector>("rtp_los")) rtp_los = *val;
    if (auto* val = ws.get<AtmPoint>("atm_point")) atm_point = *val;
  }

  auto success = std::async(std::launch::async,
                            &PropmatClearskyAgendaGUI::run,
                            std::ref(res),
                            std::ref(ctrl),
                            std::ref(ws),
                            std::cref(propmat_clearsky_agenda),
                            std::ref(jacobian_quantities),
                            std::ref(select_abs_species),
                            std::ref(f_grid),
                            std::ref(rtp_los),
                            std::ref(atm_point),
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
                     rtp_los,
                     atm_point,
                     transmission_distance,
                     ArrayOfArrayOfSpeciesTag{abs_species});
  }

  bool invalid_state = not success.get();
  ARTS_USER_ERROR_IF(invalid_state, '\n', ctrl.error)

#else   // ARTS_GUI_ENABLED not defined
  ARTS_USER_ERROR("Did not compile with GUI")
#endif  // ARTS_GUI_ENABLED
}
