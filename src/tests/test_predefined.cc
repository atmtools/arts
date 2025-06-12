#include <cstdlib>
#include <iostream>

#include "isotopologues.h"
#include "predefined_absorption_models.h"

int main() try {
  for (auto spec : Species::Isotopologues) {
    if (spec.is_predefined()) {
      if (not Absorption::PredefinedModel::can_compute(spec)) throw spec;
      else std::cout << "Can compute: " << spec.FullName() << '\n';
    } else if (Absorption::PredefinedModel::can_compute(spec)) throw spec.FullName();
  }

  return EXIT_SUCCESS;
} catch (const SpeciesIsotope& c) {
  std::cerr << "Missing implementation of computations of predefined model: "
            << c.FullName() << '\n';
  return EXIT_FAILURE;
} catch (const String& c) {
  std::cerr << "Extra implementation for computations of non-predefined model: "
            << c << '\n';
  return EXIT_FAILURE;
}