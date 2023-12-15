#include "species.h"

namespace Species {
std::ostream& operator<<(std::ostream& os, const short_name& x) {
  return os << x.name << ' ' << x.species;
}

Species toSpeciesEnumOrThrow(const std::string_view x) {
  Species n = fromShortName(x);

  if (n == Species::FINAL) {
    n = toSpecies(x);
  }

  if (n == Species::FINAL) {
    std::ostringstream os;
    os << "Invalid species name: " << std::quoted(x)
       << "\nValid options are:\n[";
    for (const auto& s : short_names_name) {
      os << '[' << s << "], ";
    }
    os << ']';
    throw std::runtime_error(os.str());
  }

  return n;
}
}  // namespace Species