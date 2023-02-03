#include "gui/plot.h"
#include "predefined_absorption_models.h"

int main(int argc, char** argv) try {
  constexpr Index minargs = 9;
  ARTS_USER_ERROR_IF(argc < minargs,
                     "Call as ",
                     argv[0],
                     " TAG F0 NF DF T P VMR_O2 VMR_H2O [GUI]");

  SpeciesTag tag(argv[1]);
  Numeric F0 = std::stod(argv[2]);
  Index NF = std::stoi(argv[3]);
  Numeric DF = std::stod(argv[4]);
  Numeric T = std::stod(argv[5]);
  Numeric P = std::stod(argv[6]);

  Absorption::PredefinedModel::VMRS vmr;
  vmr.O2 = std::stod(argv[7]);
  vmr.H2O = std::stod(argv[8]);

  bool gui = argc > minargs ? std::stoi(argv[minargs - 1]) : false;

  const Vector f_grid=uniform_grid(F0, NF, DF);
  PropagationMatrix propmat_clearsky(NF);
  ArrayOfPropagationMatrix x(0);

  Absorption::PredefinedModel::compute(
      propmat_clearsky, x, tag.Isotopologue(), f_grid, P, T, vmr, {}, {});

  if (not gui)
    std::cout << std::setprecision(15) << propmat_clearsky << '\n';
  else
    ARTSGUI::plot(f_grid, propmat_clearsky.Kjj());
  return EXIT_SUCCESS;
} catch (std::exception& e) {
  std::cerr << e.what() << '\n';
  return EXIT_FAILURE;
}
