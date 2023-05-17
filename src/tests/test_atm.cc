#include <exception>
#include <iostream>

#include "atm.h"
#include "gridded_fields.h"
#include "species.h"
#include "species_tags.h"

void point() {
  using namespace Atm;
  using enum Key;
  Point atm{pressure, 3e4, temperature, 250, ArrayOfSpeciesTag{"O2-66"}, 0.21};
  std::cout << atm << '\n' << '\n';
}

void field() {
  using namespace Atm;
  using enum Key;
  Field atm{temperature, 250.0, ArrayOfSpeciesTag{"O2-66"}, 0.21};
  atm.top_of_atmosphere = 100e3;
  
}

int main() try {
  std::cout << "\n"
               "##################################################"
               "\n"
               "START OF POINT"
               "\n"
               "\n";
  point();
  std::cout << "\n"
               "END   OF POINT"
               "\n"
               "##################################################"
               "\n"
               "\n";

  std::cout << "\n"
               "##################################################"
               "\n"
               "START OF FIELD"
               "\n"
               "\n";
  field();
  std::cout << "\n"
               "END   OF FIELD"
               "\n"
               "##################################################"
               "\n"
               "\n";
} catch (std::exception& e) {
  std::cerr << e.what() << '\n';
}