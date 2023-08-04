#include <atomic>
#include <chrono>
#include <cstdlib>
#include <functional>
#include <future>
#include <iterator>
#include <mutex>
#include <stdexcept>
#include <thread>

#include <workspace.h>
#include "artstime.h"
#include "atm.h"
#include "debug.h"
#include "jacobian.h"
#include "matpack_data.h"
#include "matpack_math.h"
#include "rtepack.h"
#include "species_tags.h"

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

      v.tm = MuelmatVector(v.pm.nelem());
      v.aotm = MuelmatMatrix(v.jacobian_quantities.nelem(), v.pm.nelem());
      auto local_aotm = v.aotm;
      Vector dr(jacobian_quantities.nelem(), 0);
      for (Index i=0 ; i<v.tm.nelem(); i++) {
        auto& t = v.tm[i];
        auto&& dt1 = v.aotm[i];
        auto&& dt2 = local_aotm[i];
        auto& k1 = v.pm[i];
        auto& k2 = v.pm[i];
        auto&& dk1 = v.aopm[i];
        auto&& dk2 = v.aopm[i];
        auto& r = v.transmission_distance;
        rtepack::two_level_exp(t, dt1, dt2, k1, k2, dk1, dk2, r, dr, dr);
      }

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
                                const Index& load [[maybe_unused]]) {
#ifdef ARTS_GUI_ENABLED
  ARTSGUI::PropmatClearsky::ResultsArray res;
  ARTSGUI::PropmatClearsky::Control ctrl;

  // Initialize values to something
  ArrayOfRetrievalQuantity jacobian_quantities{};
  ArrayOfSpeciesTag select_abs_species{};
  Vector f_grid = uniform_grid(1e9, 1000, 1e9);
  Vector rtp_los(2, 0);
  Numeric transmission_distance{1'000};
  AtmPoint atm_point;
  atm_point[Atm::Key::t] = 300;
  atm_point[Atm::Key::p] = 1000;

  for (auto& spec: abs_species) atm_point[spec] = 1.0 / static_cast<Numeric>(abs_species.nelem());

  // Set some defaults
  if (load) {
    if (ws.contains("f_grid")) f_grid = ws.get<Vector>("f_grid");
    if (ws.contains("rtp_los")) rtp_los = ws.get<Vector>("rtp_los");
    if (ws.contains("atm_point")) atm_point = ws.get<AtmPoint>("atm_point");
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
