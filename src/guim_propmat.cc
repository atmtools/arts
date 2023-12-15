#include <workspace.h>

#include <atomic>
#include <cstdlib>
#include <functional>
#include <future>
#include <mutex>
#include <thread>

#include "atm.h"
#include "debug.h"
#include "matpack_data.h"
#include "matpack_math.h"
#include "rtepack.h"
#include "species_tags.h"

#ifdef ARTS_GUI_ENABLED
#include <gui_propmat.h>

using namespace std::chrono_literals;

namespace PropmatClearskyAgendaGUI {
void compute(const Workspace& ws,
             gui::PropmatClearsky::ComputeValues& v,
             const Agenda& propmat_clearsky_agenda) {
  PropmatMatrix empty_propmat{};
  StokvecMatrix empty_stokvec{};
  propmat_clearsky_agendaExecute(ws,
                                 v.pm,
                                 v.sv,
                                 empty_propmat,
                                 empty_stokvec,
                                 {},
                                 v.select_abs_species,
                                 v.f_grid,
                                 v.rtp_los,
                                 v.atm_point,
                                 propmat_clearsky_agenda);
}

bool run(gui::PropmatClearsky::ResultsArray& ret,
         gui::PropmatClearsky::Control& ctrl,
         const Workspace& ws,
         const Agenda& propmat_clearsky_agenda,
         ArrayOfSpeciesTag& select_abs_species,
         Vector& f_grid,
         Vector& rtp_los,
         AtmPoint& atm_point,
         Numeric& transmission_distance) {
  for (auto& v : ret) v.ok.store(false);
  gui::PropmatClearsky::ComputeValues v;

  while (true) {
    std::this_thread::sleep_for(10ms);

    try {
      if (ctrl.exit.load()) return true;

      if (ctrl.error.load()) continue;

      if (ctrl.run.load()) {
        std::lock_guard allow_copy{ctrl.copy};

        v.select_abs_species = select_abs_species;
        v.f_grid = f_grid;
        v.rtp_los = rtp_los;
        v.atm_point = atm_point;
        v.transmission_distance = transmission_distance;
      } else {
        continue;
      }
      compute(ws, v, propmat_clearsky_agenda);

      v.tm = MuelmatVector(v.pm.size());
      MuelmatMatrix empty_muelmat{};
      rtepack::two_level_exp(v.tm,
                             empty_muelmat,
                             empty_muelmat,
                             v.pm,
                             v.pm,
                             {},
                             {},
                             v.transmission_distance,
                             {},
                             {});

      // Lock after the compute to copy values
      std::lock_guard allow_copy{ctrl.copy};

      // Copy over values
      ret.at(ctrl.pos).value = v;
      ret.at(ctrl.pos).ok.store(true);

      ctrl.run.store(false);
    } catch (std::exception& e) {
      ctrl.errmsg = e.what();
      ctrl.error.store(true);
      ctrl.run.store(false);
    }
  }

  return false;
}
}  // namespace PropmatClearskyAgendaGUI
#endif  // ARTS_GUI_ENABLED

void propmat_clearsky_agendaGUI(const Workspace& ws [[maybe_unused]],
                                const Agenda& propmat_clearsky_agenda
                                [[maybe_unused]],
                                const ArrayOfArrayOfSpeciesTag& abs_species
                                [[maybe_unused]],
                                const Index& load [[maybe_unused]]) {
#ifdef ARTS_GUI_ENABLED
  gui::PropmatClearsky::ResultsArray res;
  gui::PropmatClearsky::Control ctrl;

  // Initialize values to something
  ArrayOfSpeciesTag select_abs_species{};
  Vector f_grid = uniform_grid(1e9, 1000, 1e9);
  Vector rtp_los(2, 0);
  Numeric transmission_distance{1'000};
  AtmPoint atm_point;
  atm_point.temperature = 300;
  atm_point.pressure = 1000;

  for (auto& spec : abs_species)
    atm_point[spec.Species()] = 1.0 / static_cast<Numeric>(abs_species.size());

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
                            std::cref(ws),
                            std::cref(propmat_clearsky_agenda),
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
    gui::propmat(res,
                 ctrl,
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
