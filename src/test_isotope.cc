#include "isotopologues.h"

int main() {
  using namespace Species;
  
  std::cout << "Test correctness of Species fails here:\n";
  for (std::size_t i=0; i<Isotopologues.size(); i++) {
    auto& x = Isotopologues[i];
    if (not good_enum(x.spec)) {
      std::cout << "Position: " << i << " has bad values\n";
    }
  }
  
  std::cout << "Test correctness of Species short-name conversion fails here:\n";
  for (Index i=0; i<Index(Species::Species::FINAL); i++) {
    auto a = Species::Species(i);
    auto b = Species::toShortName(a);
    auto c = Species::fromShortName(b);
    if (not good_enum(c) or c not_eq a) {
      std::cout << i << ' ' << a << ' ' << b << ' ' << c << '\n';
    }
  }
  
  std::cout << "\n\nTest of order of Species fails for these pairs:\n";
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
