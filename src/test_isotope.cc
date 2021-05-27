#include "isotopologues.h"

int main() {
  using namespace Species;
  
  std::cout << "Test of order of Species fails for these pairs:\n";
  for (std::size_t i=0; i<Isotopologues.size()-1; i++) {
    auto& a = Isotopologues[i];
    auto& b = Isotopologues[i+1];
    if (std::size_t(a.spec) > std::size_t(b.spec)) {
      std::cout << a << ' ' << b << '\n';
    }
  }
  
  std::cout << "\n\nTest of order of Isotopes fails for these pairs:\n";
  for (std::size_t i=0; i<Isotopologues.size()-1; i++) {
    auto& a = Isotopologues[i];
    auto& b = Isotopologues[i+1];
    if (a.spec == b.spec) {
      if (a.isotname.compare(b.isotname) >= 0) {
        std::cout << a << ' ' << b << '\n';
      }
    }
  }
}
