#include <array>
#include "hitran_species.h"

void test001() {
  const Index nmols = 100;
  std::array<char, 12> isos = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '0', 'A', 'B'};
  
  for (Index mol = 0; mol < nmols; mol++) {
    for (char iso : isos) {
      try {
        Hitran::id_from_lookup(mol, iso, Hitran::Type::Newest);
      } catch(std::runtime_error& e) {
        std::cout << e.what() << "\n\n";
      }
    }
  }
}

int main() {
  test001();
}
