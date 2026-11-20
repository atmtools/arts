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
const JplSpeciesInfo& info_from_lookup(Index mol, std::span<const JplSpeciesInfo> data) {
  auto it = stdr::lower_bound(data, mol, {}, &JplSpeciesInfo::id);

  if (it == data.end() or it->id != mol) {
    throw std::runtime_error(std::format("No species found in JPL data with ID {}", mol));
  }

  return *it;
}

QuantumIdentifier qid_from_data(const JplSpeciesInfo& info) { return QuantumIdentifier(info.spec); }
}  // namespace

std::span<const JplSpeciesInfo> builtin_data() { return jpl_data; }

ArrayOfJplSpeciesInfo sorted_data(ArrayOfJplSpeciesInfo data) {
  stdr::sort(data, {}, &JplSpeciesInfo::id);

  if (auto it = stdr::adjacent_find(data, {}, &JplSpeciesInfo::id); it != data.end())
    arts::user_error("The JPL species tag {} is given more than once", it->id);

  for (const auto& info : data) {
    if (not(info.QT0 > 0 and info.T0 > 0))
      arts::user_error(
          "The JPL species tag {} must have a positive partition function and reference temperature, got QT0 = {} and T0 = {}",
          info.id,
          info.QT0,
          info.T0);
  }

  return data;
}

LineDataMod data_lookup(Index mol, std::span<const JplSpeciesInfo> data) try {
  const auto& info = info_from_lookup(mol, data);
  return LineDataMod{.qid = qid_from_data(info), .T0 = info.T0, .QT0 = info.QT0};
} catch (const std::exception& e) {
  arts::user_error("Error occurred while looking up species index {}:\n{}", mol, e.what());
}
}  // namespace Jpl
