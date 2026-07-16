#include "jpl_species.h"

#include <isotopologues.h>

#include <algorithm>
#include <cassert>

//  Contributes the molparam_map
#include "auto_jpl_species_map.h"
#include "debug.h"
#include "quantum.h"

namespace Jpl {
namespace {
const JplSpeciesInfo& info_from_lookup(Index mol) {
  auto it = std::ranges::lower_bound(jpl_data, mol, {}, &JplSpeciesInfo::id);

  if (it == jpl_data.end() or it->id != mol) {
    throw std::runtime_error(std::format("No species found in JPL data with ID {}", mol));
  }

  return *it;
}

QuantumIdentifier qid_from_data(const JplSpeciesInfo& info) { return QuantumIdentifier(info.spec); }
}  // namespace

LineDataMod data_lookup(Index mol) try {
  const auto& info = info_from_lookup(mol);
  return LineDataMod{.qid = qid_from_data(info), .T0 = info.T0, .QT0 = info.QT0};
} catch (const std::exception& e) {
  ARTS_USER_ERROR("Error occurred while looking up species index {}:\n{}", mol, e.what());
}
}  // namespace Jpl
