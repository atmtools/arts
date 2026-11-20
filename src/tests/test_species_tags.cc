#include <debug.h>
#include <species_tags.h>

#include <cstdlib>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) try {
  if (argc not_eq 2) arts::user_error("Call as {} SPECIES_STRING", argv[0]);

  const SpeciesTag x(argv[1]);
  arts::println("{} {} {}", x, x.Isotopologue(), x.type);

  return EXIT_SUCCESS;
} catch (std::runtime_error& e) {
  arts::println(std::cerr, "{}", e.what());
  return EXIT_FAILURE;
}
