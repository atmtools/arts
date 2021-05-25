#include <iostream>

#include "species_tags.h"

int main(int argc, char ** argv) {
  ARTS_USER_ERROR_IF(argc not_eq 2, "Call as ", argv[0], " SPECIES_STRING");
  
  const ArrayOfSpeciesTag tags(argv[1]);
  std::cout << tags << ' ' << tags.Species() << '\n';
  for (auto& x: tags) {
    std::cout << x << ' ' << x.Isotopologue() << ' ' << x.type << '\n';
  }
}
