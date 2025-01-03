#pragma once

#include <array.h>
#include <enumsSpeciesEnum.h>
#include <format_tags.h>

using ArrayOfSpeciesEnum        = Array<SpeciesEnum>;
using ArrayOfArrayOfSpeciesEnum = Array<ArrayOfSpeciesEnum>;

consteval SpeciesEnum operator""_spec(const char* x, std::size_t) {
  return to<SpeciesEnum>(x);
}
