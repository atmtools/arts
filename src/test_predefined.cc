#include <cstdlib>

#include "isotopologues.h"
#include "predefined_absorption_models.h"
#include "propagationmatrix.h"

int main() try {
  PropagationMatrix pm;

  for (auto spec : Species::Isotopologues) {
    if (Species::is_predefined_model(spec)) {
      if (not Absorption::PredefinedModel::compute_selection<true>(
              pm, spec, {}, {}, {}, {}, {}))
        throw spec;
    } else if (Absorption::PredefinedModel::compute_selection<true>(
              pm, spec, {}, {}, {}, {}, {}))
        throw spec.FullName();
  }

  return EXIT_SUCCESS;
} catch (const SpeciesIsotopeRecord& c) {
  std::cerr << "Missing implementation of computations of predefined model: "
            << c.FullName() << '\n';
  return EXIT_FAILURE;
} catch (const String& c) {
  std::cerr << "Extra implementation for computations of non-predefined model: "
            << c << '\n';
  return EXIT_FAILURE;
}