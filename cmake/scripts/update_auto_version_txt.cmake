set (CMAKE_MODULE_PATH "${ARTS_SOURCE_DIR}/cmake/modules")

include (ArtsVersion)

if (EXISTS "auto_version.txt")
  file (READ "auto_version.txt" ARTS_PREVIOUS_VERSION)
endif (EXISTS "auto_version.txt")

ARTS_GET_VERSION (ARTS_VERSION)

if (NOT "${ARTS_VERSION}" STREQUAL "${ARTS_PREVIOUS_VERSION}")
  file(WRITE "auto_version.txt" "${ARTS_VERSION}")
endif (NOT "${ARTS_VERSION}" STREQUAL "${ARTS_PREVIOUS_VERSION}")

message (STATUS "ARTS version: ${ARTS_VERSION}")

