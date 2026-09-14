#pragma once

#include <debug.h>

#include <mutex>
#include <vector>

// The same weak template-static objects must be emitted by both translation
// units. GCC/Darwin previously coalesced their adjacent source-location records,
// making error reporting use another function's metadata or invalid pointers.
template <typename T> void source_location_record_use(T value) {
  static std::mutex mutex;
  static std::vector<T> values;
  std::lock_guard lock(mutex);
  values.push_back(value);
}

src_location other_source_location();
