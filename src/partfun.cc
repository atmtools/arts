#include <iostream>
#include "partfun.h"

int main(int argc, char **argv) try {
  ARTS_USER_ERROR_IF(argc not_eq 4, "Call as ", argv[0], " SPECIES ISOTOPE TEMPERATURE");
  
  Species::Species SPECIES = Species::fromShortName(argv[1]);
  ARTS_USER_ERROR_IF(not good_enum(SPECIES), argv[1], " is not a recognized species")
  std::string_view ISOTOPE = argv[2];
  Numeric TEMPERATURE = std::stod(argv[3]);
  ARTS_USER_ERROR_IF(TEMPERATURE<=0, "Must have positive temperature")
  
  auto isotope = Species::IsotopeRecord(SPECIES, ISOTOPE);
  std::cout << std::setprecision(15) << isotope.FullName() << ' ' << "Q(" << TEMPERATURE << ") " << PartitionFunctions::Q(TEMPERATURE, isotope) << '\n';
  return EXIT_SUCCESS;
} catch(std::exception& e) {
  std::cerr << e.what() << '\n';
  return EXIT_FAILURE;
}
