#pragma once

#include <enumsPartitionFunctionsType.h>
#include <matpack.h>

namespace PartitionFunctions {
struct Data {
  PartitionFunctionsType type;
  Matrix data;
};
}  // namespace PartitionFunctions

using PartitionFunctionsData = PartitionFunctions::Data;
