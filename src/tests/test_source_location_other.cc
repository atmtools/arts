#include "test_source_location.h"

src_location other_source_location() {
  source_location_record_use(2);
  return {};
}
