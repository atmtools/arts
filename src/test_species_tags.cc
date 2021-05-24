#include <iostream>

#include "species_tags.h"

int main(int argc, char ** argv) {
  ARTS_USER_ERROR_IF(argc not_eq 2, "Call as ", argv[0], " SPECIES_STRING");
  std::cout << ArrayOfSpeciesTag(argv[1]) << '\n';
}
