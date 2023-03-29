#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include "species_tags.h"

int main(int argc, char ** argv) try {
  ARTS_USER_ERROR_IF(argc not_eq 2, "Call as ", argv[0], " SPECIES_STRING");
  
  const ArrayOfSpeciesTag tags(argv[1]);
  for (auto& x: tags) {
    std::cout << x << ' ' << x.Isotopologue() << ' ' << x.type << '\n';
  }

  return EXIT_SUCCESS;
} catch (std::runtime_error& e) {
  std::cerr << e.what() << '\n';
  return EXIT_FAILURE;
}
