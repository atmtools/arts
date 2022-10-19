#include <iostream>

#include "atm.h"
#include "species.h"
#include "species_tags.h"

void point() {
  using namespace Atm;
  using enum Key;
  Point pnt {pressure, 3e4, temperature, 250, ArrayOfSpeciesTag{"O2-66"}, 0.21};
  std::cout << pnt << '\n'<< '\n';

  pnt.set(ArrayOfSpeciesTag{"H2O-161"},0.01);
  std::cout << pnt << '\n'<< '\n';

  std::cout << pnt[ArrayOfSpeciesTag{"O2-66"}] << '\n'<< '\n';
  std::cout << pnt[pressure] << '\n'<< '\n';
}

void field() {
}

int main() {
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
}