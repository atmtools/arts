#include "test_source_location.h"

#include <cstdio>
#include <exception>
#include <stdexcept>
#include <string_view>

src_location first_source_location() {
  source_location_record_use(1);
  return {};
}

int main() try {
  auto first = first_source_location();
  auto other = other_source_location();
  if (not first.get().contains("test_source_location.cc") or not first.getfunc().contains("first_source_location"))
    throw std::runtime_error("First translation unit lost its source location");
  if (not other.get().contains("test_source_location_other.cc") or
      not other.getfunc().contains("other_source_location"))
    throw std::runtime_error("Second translation unit reused the first source location");
  return 0;
} catch (const std::exception& error) {
  std::fprintf(stderr, "%s\n", error.what());
  return 1;
}
