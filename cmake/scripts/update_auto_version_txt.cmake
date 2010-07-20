set (CMAKE_MODULE_PATH "${ARTS_SOURCE_DIR}/cmake/modules")

include (ArtsVersion)

if (EXISTS "auto_version.txt")
  ARTS_GET_VERSION_FROM_AUTOVERSION (ARTS_PREVIOUS_VERSION "auto_version.txt")
endif (EXISTS "auto_version.txt")

ARTS_GET_VERSION_FROM_CHANGELOG (ARTS_VERSION "${ARTS_SOURCE_DIR}/ChangeLog")

if (NOT "${ARTS_VERSION}" STREQUAL "${ARTS_PREVIOUS_VERSION}")
  file(WRITE "auto_version.txt" "${ARTS_VERSION}")
endif (NOT "${ARTS_VERSION}" STREQUAL "${ARTS_PREVIOUS_VERSION}")

message (STATUS "ARTS version: ${ARTS_VERSION}")

