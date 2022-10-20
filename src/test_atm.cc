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

  atm.set(ArrayOfSpeciesTag{"H2O-161"}, 0.01);
  std::cout << atm << '\n' << '\n';

  std::cout << atm[ArrayOfSpeciesTag{"O2-66"}] << '\n' << '\n';
  std::cout << atm[pressure] << ' ' << atm.P() << '\n' << '\n';
}

void field() {
  using namespace Atm;
  using enum Key;
  Field atm{
      pressure, 3e4, temperature, 250.0, ArrayOfSpeciesTag{"O2-66"}, 0.21};
  std::cout << atm << '\n' << '\n';

  atm.set(ArrayOfSpeciesTag{"H2O-161"}, 0.01);
  std::cout << atm << '\n' << '\n';

  GriddedField4 data;
  data.set_grid(0, {Time{}});
  data.set_grid_name(0, "Time");
  data.set_grid(1, {0, 2, 4});
  data.set_grid_name(1, "Altitude");
  data.set_grid(2, {0, 2, 4});
  data.set_grid_name(2, "Latitude");
  data.set_grid(3, {0, 2, 4});
  data.set_grid_name(3, "Longitude");
  data.data.resize(1, 3, 3, 3);
  for (Index i1 = 0; i1 < data.get_grid_size(0); i1++) {
    for (Index i2 = 0; i2 < data.get_grid_size(1); i2++) {
      for (Index i3 = 0; i3 < data.get_grid_size(2); i3++) {
        for (Index i4 = 0; i4 < data.get_grid_size(3); i4++) {
          data.data(i1, i2, i3, i4) = static_cast<Numeric>(1 + i1 + 2 * i2 + 3 * i3 + 4 * i4);
        }
      }
    }
  }
  atm.set(mag_u, data);

  std::cout << atm.at(Time{}, 0.0, 3.0, 0.0) << '\n' << '\n';

  std::cout << atm.regularize({Time{}}, {0.0, 1.0}, {0.5}, {0.5}) << '\n'
            << '\n';
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