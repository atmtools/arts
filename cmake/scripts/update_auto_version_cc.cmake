list (APPEND CMAKE_MODULE_PATH "${ARTS_SOURCE_DIR}/cmake/modules")

include (ArtsVersion)

if (EXISTS "auto_version_cc.txt")
  file (READ "auto_version_cc.txt" ARTS_PREVIOUS_VERSION)
endif (EXISTS "auto_version_cc.txt")

ARTS_GET_VERSION(ARTS_VERSION)

if (NOT "${ARTS_VERSION}" STREQUAL "${ARTS_PREVIOUS_VERSION}")
  file(WRITE "auto_version_cc.txt" "${ARTS_VERSION}")
  ARTS_WRITE_VERSION_FILE ("auto_version.cc"
  "
  #include <string>

  std::string_view arts_get_version_string() {
    static const std::string arts_version{\"arts-\${ARTS_VERSION}\"};
    return arts_version;
  }
  "
   )
endif()

