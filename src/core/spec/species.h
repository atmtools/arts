#pragma once

#include <array.h>
#include <format_tags.h>
#include <enumsSpeciesEnum.h>

using ArrayOfSpeciesEnum        = Array<SpeciesEnum>;
using ArrayOfArrayOfSpeciesEnum = Array<ArrayOfSpeciesEnum>;

std::ostream& operator<<(std::ostream& os, const ArrayOfSpeciesEnum& a);
std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfSpeciesEnum& a);
