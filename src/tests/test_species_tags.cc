#include <debug.h>
#include <species_tags.h>

#include <cstdlib>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) try {
  ARTS_USER_ERROR_IF(argc not_eq 2, "Call as {} SPECIES_STRING", argv[0]);

  const SpeciesTag x(argv[1]);
  std::println("{} {} {}", x, x.Isotopologue(), x.type);

  return EXIT_SUCCESS;
} catch (std::runtime_error& e) {
  std::println(std::cerr, "{}", e.what());
  return EXIT_FAILURE;
}
