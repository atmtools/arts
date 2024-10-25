#include "debug.h"
#include "isotopologues.h"
#include "partfun.h"
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

int main(int argn, char** argv) try {
  ARTS_USER_ERROR_IF(argn < 3, "USAGE: PROG SPEC1... SPECN T")

  const std::vector<std::string> specs{argv+1, argv+argn-1};
  const Numeric T = std::strtod(argv[argn-1], nullptr);

  for (auto& spec: specs) {
    auto& v = Species::Isotopologues[Species::find_species_index(spec)];
    std::cout << v.FullName() << ' ';
    std::cout << PartitionFunctions::Q(T, v) << '\n';
  }

  return EXIT_SUCCESS;
} catch(std::exception& e) {
  std::cerr << "Failed with error:\n" << e.what() << '\n';
  return EXIT_FAILURE;
}