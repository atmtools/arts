
#include <string>

#include "config.h"

std::string_view arts_get_version_string() {
  static const std::string arts_version{"arts-" ARTS_VERSION};
  return arts_version;
}
