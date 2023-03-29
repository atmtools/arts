list (APPEND CMAKE_MODULE_PATH "${ARTS_SOURCE_DIR}/cmake/modules")

include (ArtsVersion)

ARTS_WRITE_VERSION_FILE ("auto_version.cc"
"
#include <string>

std::string_view arts_get_version_string() {
  static const std::string arts_version{\"arts-\${ARTS_VERSION}\"};
  return arts_version;
}
"
 )

