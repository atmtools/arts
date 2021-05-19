#include "isotopologues.h"

#ifndef NDEBUG

namespace CompileTimeTests {
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
constexpr bool testIsotopologuesNotRepeating() noexcept {
  for (std::size_t i=0; i<Isotopologues.size(); i++) {
    for (std::size_t j=i+1; j<Isotopologues.size(); j++) {
      if (Isotopologues[i] == Isotopologues[j]) {
        return false;
      }
    }
  }
  return true;
}

static_assert(testIsotopologuesNotRepeating(), "Error!\n\n"
              "Cannot have repeating isotopologues in the isotopologues list.\n"
              "One of your newly added or modified species already exists.\n"
              "Ensure that there's no repetition!\n");

//! Don't call this manually, it only exists to catch a developer error
constexpr bool testIsotopologuesAllValid() noexcept {
  for (auto& x: Isotopologues) if (not good_enum(x.spec)) return false;
  return true;
}

static_assert(testIsotopologuesAllValid(), "Error!\n\n"
              "Cannot have invalid species in the isotopologues list.\n"
              "One of your newly added or modified species cannot be understood\n");
}

#endif
