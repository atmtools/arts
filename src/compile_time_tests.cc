#include "isotopologues.h"

#ifndef NDEBUG

namespace CompileTimeTests {
//! Don't call this manually, it only exists to catch a developer error
constexpr bool testIsotopologuesAllValid() noexcept {
  using namespace Species;
  for (auto& x: Isotopologues) if (not good_enum(x.spec)) return false;
  return true;
}

static_assert(testIsotopologuesAllValid(), "Error!\n\n"
              "Cannot have invalid species in the isotopologues list.\n"
              "One of your newly added or modified species cannot be understood\n");
}

//! Don't call this manually, it only exists to catch a developer error
constexpr bool testShortNames() noexcept {
  for (Index i=0; i<Index(Species::Species::FINAL); i++) {
    auto a = Species::Species(i);
    auto b = Species::toShortName(a);
    auto c = Species::fromShortName(b);
    if (not good_enum(c) or c not_eq a) {
      return false;
    }
  }
  return true;
}

static_assert(testShortNames(), "Error!\n\n"
              "Cannot convert some species from to its short name and back to its enum-value.\n"
              "All species in Species::Species must be translatable to and from short names in\n"
              "the two function signatures:\n\n"
              "\tstd::string_view toShortName(Species x)\n"
              "\tSpecies fromShortName(const std::string_view x)\n"
              "\n"
              "Please ensure that any new Species you've added is available in these functions\n"
              "and in all functions that makes use of a Species-switch.\n");

//! Don't call this manually, it only exists to catch a developer error
constexpr bool testSpeciesIncreasing() noexcept {
  using namespace Species;
  static_assert(Isotopologues.size() not_eq 0);
  for (std::size_t i=0; i<Isotopologues.size()-1; i++) {
    auto& a = Isotopologues[i];
    auto& b = Isotopologues[i+1];
    if (std::size_t(a.spec) > std::size_t(b.spec)) {
      //! First Species Must Be Lower Than Second
      return false;
    }
  }
  return true;
}

static_assert(testSpeciesIncreasing(), "Error!\n\n"
"Species in Isotopologues must be increasing.\n"
"One of your newly added isotopologues is not in increasing order\n");

//! Don't call this manually, it only exists to catch a developer error
constexpr bool testIsotopologuesIncreasing() noexcept {
  using namespace Species;
  static_assert(Isotopologues.size() not_eq 0);
  for (std::size_t i=0; i<Isotopologues.size()-1; i++) {
    auto& a = Isotopologues[i];
    auto& b = Isotopologues[i+1];
    if (a.spec == b.spec) {
      if (a.isotname.compare(b.isotname) >= 0) {
        //! First isotope Of A Species Must Be Different And Lower
        return false;
      }
    }
  }
  return true;
}

static_assert(testIsotopologuesIncreasing(), "Error!\n\n"
              "Isotopologues must be increasing.\n"
              "One of your newly added isotopologues is not in increasing order\n");

#endif
