list (APPEND CMAKE_MODULE_PATH "${ARTS_SOURCE_DIR}/cmake/modules")

include (ArtsVersion)

ARTS_WRITE_VERSION_FILE ("auto_version.h" "#define ARTS_FULL_VERSION \"arts-\${ARTS_VERSION}\"\n")

