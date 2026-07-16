#pragma once

#include <quantum.h>

namespace Jpl {
struct LineDataMod {
  QuantumIdentifier qid;
  Numeric           T0;
  Numeric           QT0;
};

/** Finds the ID of the ARTS species from JPL */
LineDataMod data_lookup(Index mol);
}  // namespace Jpl
