#include "isotopologues.h"

int main() {
  std::cout << "Test correctness of Species fails here:\n";
  for (std::size_t i=0; i<Species::Isotopologues.size(); i++) {
    auto& x = Species::Isotopologues[i];
    if (not good_enum(x.spec)) {
      std::cout << "Position: " << i << " has bad values\n";
    }
  }
  
  std::cout << "\n\nTest correctness of Species short-name conversion fails here:\n";
  for (Index i=0; i<Index(Species::Species::FINAL); i++) {
    auto a = Species::Species(i);
    auto b = Species::toShortName(a);
    auto c = Species::fromShortName(b);
    if (not good_enum(c) or c not_eq a) {
      std::cout << i << ' ' << a << ' ' << b << ' ' << c << '\n';
    }
  }
  
  std::cout << "\n\nTest of order of Species fails for these pairs:\n";
  for (std::size_t i=0; i<Species::Isotopologues.size()-1; i++) {
    auto& a = Species::Isotopologues[i];
    auto& b = Species::Isotopologues[i+1];
    if (std::size_t(a.spec) > std::size_t(b.spec)) {
      std::cout << a.FullName() << ' ' << b.FullName() << '\n';
    }
  }
  
  std::cout << "\n\nTest of order of Isotopologues fails for these pairs:\n";
  for (std::size_t i=0; i<Species::Isotopologues.size()-1; i++) {
    auto& a = Species::Isotopologues[i];
    auto& b = Species::Isotopologues[i+1];
    if (a.spec == b.spec) {
      if (a.isotname.compare(b.isotname) >= 0) {
        std::cout << a.FullName() << ' ' << b.FullName() << '\n';
      }
    }
  }
  
  std::cout << "\n\nTest that all Isotopologues (that are explicit isotopes) has a ratio in the builtin IsotopologueRatios fails here:\n";
  auto iso_rat = Species::isotopologue_ratiosInitFromBuiltin();
  for (Index i=0; i<iso_rat.maxsize; i++) {
    if (not is_predefined_model(Species::Isotopologues[i]) and not Species::Isotopologues[i].joker() and nonstd::isnan(iso_rat.data[i])) {
      std::cout << Species::Isotopologues[i].FullName() << " has no isotopologue ratio by default: " << iso_rat.data[i] << '\n';
    }
  }
}
