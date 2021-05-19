#include <iostream>

#include "abs_species_tags.h"
#include "species_tags.h"

int main(int argc, char ** argv) {
  define_species_data();
  define_species_map();
  
  ARTS_USER_ERROR_IF(argc not_eq 2, "Call as ", argv[0], " SPECIES_STRING");
  String spec(argv[1]);
  ArrayOfSpeciesTag specs;
  array_species_tag_from_string(specs, spec);
  std::cout << specs << ' ' << ArrayOfSpeciesTag2(spec) << '\n';
}
