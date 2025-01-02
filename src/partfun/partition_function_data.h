#pragma once

#include <matpack.h>
#include <enumsPartitionFunctionsType.h>

namespace PartitionFunctions {
struct Data {
  PartitionFunctionsType type;
  Matrix data;
};
}  // namespace PartitionFunctions

using PartitionFunctionsData = PartitionFunctions::Data;
