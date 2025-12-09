#include <debug.h>
#include <isotopologues.h>
#include <partfun.h>

#include <iostream>

namespace {
//! Don't call this manually, it only exists to catch a developer error
constexpr std::size_t nonexistentPartfun() noexcept {
  for (std::size_t i = 0; i < Species::Isotopologues.size(); i++) {
    auto& ir = Species::Isotopologues[i];
    if (not ir.is_predefined() and not ir.is_joker()) {
      if (not PartitionFunctions::has_partfun(ir)) {
        return i;
      }
    }
  }
  return Species::Isotopologues.size();
}
}  // namespace

int main() try {
  const auto i = nonexistentPartfun();
  ARTS_USER_ERROR_IF(
      i not_eq Species::Isotopologues.size(),
      "A species without partition functions have been found.  It is the species: {}",
      Species::Isotopologues[i].FullName())
  return EXIT_SUCCESS;
} catch (std::runtime_error& e) {
  std::cerr << e.what() << '\n';
  return EXIT_FAILURE;
}
