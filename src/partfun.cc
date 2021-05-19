#include <iostream>
#include "partfun.h"

int main(int argc, char **argv) try {
  ARTS_USER_ERROR_IF(argc not_eq 4, "Call as ", argv[0], " SPECIES ISOTOPE TEMPERATURE");
  
  Species::Species spec = Species::fromShortName(argv[1]);
  ARTS_USER_ERROR_IF(not good_enum(spec), argv[1], " is not a recognized species")
  std::string_view isot = argv[2];
  Numeric T = std::stod(argv[3]);
  ARTS_USER_ERROR_IF(T<=0, "Must have positive temperature")
  
  auto isotope = Species::IsotopeRecord(spec, isot);
  std::cout << std::setprecision(15) << isotope.FullName() << ' ' << "Q(" << T << ") " << PartitionFunctions::Q(T, isotope) << '\n';
  return EXIT_SUCCESS;
} catch(std::exception& e) {
  std::cerr << e.what() << '\n';
  return EXIT_FAILURE;
}
