#include <array>
#include "hitran_species.h"

void test001() {
  define_species_data();
  define_species_map();
  
  const Index nmols = 100;
  std::array<char, 12> isos = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '0', 'A', 'B'};
  
  for (Index mol = 0; mol < nmols; mol++) {
    for (char iso : isos) {
      try {
        Hitran::from_lookup(mol, iso);
      } catch(std::runtime_error& e) {
        std::cout << e.what() << "\n\n";
      }
    }
  }
}

int main() {
  test001();
}
