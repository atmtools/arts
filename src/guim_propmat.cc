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
#include "path_point.h"
#include "rtepack.h"
#include "sorted_grid.h"
#include "species.h"
#include "species_tags.h"

#ifdef ARTS_GUI_ENABLED
#include <gui_propmat.h>

using namespace std::chrono_literals;

namespace PropmatClearskyAgendaGUI {
void compute(const Workspace& ws,
             gui::PropmatClearsky::ComputeValues& v,
             const Agenda& propagation_matrix_agenda) {
  PropmatMatrix empty_propmat{};
  StokvecMatrix empty_stokvec{};
  propagation_matrix_agendaExecute(ws,
                                   v.pm,
                                   v.sv,
                                   empty_propmat,
                                   empty_stokvec,
                                   {},
                                   v.select_species,
                                   AscendingGrid{v.f_grid},
                                   v.path_point,
                                   v.atm_point,
                                   propagation_matrix_agenda);
}

bool run(gui::PropmatClearsky::ResultsArray& ret,
         gui::PropmatClearsky::Control& ctrl,
         const Workspace& ws,
         const Agenda& propagation_matrix_agenda,
         SpeciesEnum& select_species,
         Vector& f_grid,
         PropagationPathPoint& path_point,
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

        v.select_species = select_species;
        v.f_grid = f_grid;
        v.path_point = path_point;
        v.atm_point = atm_point;
        v.transmission_distance = transmission_distance;
      } else {
        continue;
      }
      compute(ws, v, propagation_matrix_agenda);

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

void propagation_matrix_agendaGUI(const Workspace& ws [[maybe_unused]],
                                  const Agenda& propagation_matrix_agenda
                                  [[maybe_unused]],
                                  const ArrayOfArrayOfSpeciesTag& abs_species
                                  [[maybe_unused]],
                                  const Index& load [[maybe_unused]]) {
#ifdef ARTS_GUI_ENABLED
  gui::PropmatClearsky::ResultsArray res;
  gui::PropmatClearsky::Control ctrl;

  // Initialize values to something
  SpeciesEnum select_species{};
  Vector f_grid = uniform_grid(1e9, 1000, 1e9);
  PropagationPathPoint path_point{.pos = {0, 0, 0}, .los = {0, 0}};
  Numeric transmission_distance{1'000};
  AtmPoint atm_point;
  atm_point.temperature = 300;
  atm_point.pressure = 1000;

  for (auto& spec : abs_species)
    atm_point[spec.Species()] = 1.0 / static_cast<Numeric>(abs_species.size());

  // Set some defaults
  if (load) {
    if (ws.wsv_and_contains("frequency_grid")) f_grid = ws.get<AscendingGrid>("frequency_grid");
    if (ws.wsv_and_contains("atmospheric_point")) atm_point = ws.get<AtmPoint>("atmospheric_point");
    if (ws.wsv_and_contains("propagation_path_point"))
      path_point = ws.get<PropagationPathPoint>("propagation_path_point");
  }

  auto success = std::async(std::launch::async,
                            &PropmatClearskyAgendaGUI::run,
                            std::ref(res),
                            std::ref(ctrl),
                            std::cref(ws),
                            std::cref(propagation_matrix_agenda),
                            std::ref(select_species),
                            std::ref(f_grid),
                            std::ref(path_point),
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
                 select_species,
                 f_grid,
                 path_point,
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
