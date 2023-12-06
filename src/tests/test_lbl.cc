#include <lbl.h>

#include <cstdlib>
#include <exception>
#include <limits>
#include <stdexcept>

#include "atm.h"
#include "gui_plot.h"
#include "isotopologues.h"
#include "lbl_data.h"
#include "lbl_zeeman.h"
#include "math_funcs.h"
#include "matpack_math.h"
#include "new_jacobian.h"
#include "partfun.h"
#include "quantum_numbers.h"
#include "rtepack.h"
#include "species.h"
#include "wigner_functions.h"

AbsorptionBands bands(
    const bool one_by_one = false,
    const lbl::CutoffType cutoff_type = lbl::CutoffType::None,
    const Numeric cutoff = std::numeric_limits<Numeric>::infinity()) {
  lbl::line line;
  line.a = 9.123e-10;
  line.f0 = 61150568042.6513;
  line.e0 = 2.55055079232451e-21;
  line.gu = 19;
  line.gl = 21;
  line.z = lbl::zeeman::data{0.0223601, 0.2001645};
  line.ls.T0 = 296;
  line.ls.one_by_one = one_by_one;

  using enum lbl::temperature::model_type;
  line.ls.single_models.emplace_back(
      Species::Species::Oxygen,
      std::vector{
          std::pair{lbl::line_shape::variable::G0,
                    lbl::temperature::data{T1, {13314.2468393782, 0.72}}},
          std::pair{lbl::line_shape::variable::Y,
                    lbl::temperature::data{POLY,
                                           {-8.16811474067715e-05,
                                            8.50048986370598e-07,
                                            -2.94005747764421e-09,
                                            3.39745297443195e-12}}}});
  line.ls.single_models.emplace_back(
      Species::Species::Bath,
      std::vector{
          std::pair{lbl::line_shape::variable::G0,
                    lbl::temperature::data{T1, {13136.7235481865, 0.72}}},
          std::pair{lbl::line_shape::variable::Y,
                    lbl::temperature::data{POLY,
                                           {-7.04762100774315e-05,
                                            7.34353279586519e-07,
                                            -2.53727010350522e-09,
                                            2.92959627702933e-12}}}});
  line.qn = QuantumNumberLocalState{QuantumNumberValue{"J 10 9"},
                                    QuantumNumberValue{"N  9 9"}};

  lbl::band band;
  band.key = QuantumIdentifier{"O2-66 ElecStateLabel X X Lambda 0 0 S 1 1 v 0 0"};
  band.data.emplace_back(std::move(line));
  band.data.lineshape = lbl::Lineshape::VP_LTE;
  band.data.cutoff = cutoff_type;
  band.data.cutoff_value = cutoff;

  return {band};
}

AtmPoint atm() {
  AtmPoint atm;
  atm.pressure = 1e5;
  atm.temperature = 196;
  atm[Species::Species::Oxygen] = 0.21;
  atm[Species::Species::Nitrogen] = 0.79;
  atm.mag = {10e-6, 20e-6, 40e-6};
  atm[SpeciesIsotopeRecord("O2-66")] = .995262E+00;
  return atm;
}

void test_voigt_no_cutoff(const bool one_by_one = false) {
  const auto bands = ::bands(one_by_one,
                             lbl::CutoffType::None,
                             std::numeric_limits<Numeric>::infinity());
  const AtmPoint atm = ::atm();

  Vector f_grid;
  linspace(f_grid, 60e9, 63e9, 10e3);
  StokvecVector sv(f_grid.size());
  PropmatVector pm(f_grid.size());
  PropmatMatrix dpm(0, f_grid.size());
  StokvecMatrix dsv(0, f_grid.size());

  make_wigner_ready(100, 100, 6);
  lbl::calculate(pm, sv, dpm, dsv, f_grid, {}, bands, {}, atm);

  Vector vpm(f_grid.size());
  for (Index i = 0; i < f_grid.size(); ++i) {
    vpm[i] = pm[i][0];
  }

  gui::plot(f_grid, vpm);
}

int main() try {
  test_voigt_no_cutoff();
  return EXIT_SUCCESS;
} catch (std::exception& e) {
  std::cerr << e.what() << std::endl;
  return EXIT_FAILURE;
}
