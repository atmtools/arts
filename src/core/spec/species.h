#pragma once

#include <array.h>

#include "enums.h"

using ArrayOfSpeciesEnum = Array<SpeciesEnum>;
using ArrayOfArrayOfSpeciesEnum = Array<ArrayOfSpeciesEnum>;

std::ostream& operator<<(std::ostream& os, const ArrayOfSpeciesEnum& a);
std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfSpeciesEnum& a);
