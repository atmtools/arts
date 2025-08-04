#pragma once

#include <vector>

#include "configtypes.h"

/** An array. */
template <typename base>
using Array = std::vector<base>;

/** An array of Index. */
using ArrayOfIndex = Array<Index>;

using ArrayOfArrayOfIndex = Array<ArrayOfIndex>;

/** An array of Numeric. */
using ArrayOfNumeric = Array<Numeric>;