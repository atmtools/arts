#include <isotopologues.h>

#include <algorithm>

int main() {
  const SpeciesIsotope* const first = Species::Isotopologues.begin();
  const SpeciesIsotope* const last  = Species::Isotopologues.end();

  const auto* const anybad = std::ranges::adjacent_find(
      first, last, [](SpeciesIsotope x, SpeciesIsotope y) { return x >= y; });

  if (anybad == last) return 0;

  std::vector<String> printing;
  if (anybad != first) printing.push_back((anybad - 1)->FullName());
  printing.push_back((anybad)->FullName());
  if ((anybad + 1) != last) printing.push_back((anybad + 1)->FullName());

  std::print(R"(The isotopologues are not sorted.

Several internal methods assume that they must be sorted.

The offending species is in the list:

{:B,qN}
)",
             printing);

  return 1;
}