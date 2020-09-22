list (APPEND CMAKE_MODULE_PATH "${ARTS_SOURCE_DIR}/cmake/modules")

include (ArtsVersion)

file (READ "${ARTS_BINARY_DIR}/auto_version.txt" VERSION)
configure_file (${ARTS_SOURCE_DIR}/doc/index.html.in
                ${ARTS_BINARY_DIR}/doc/index.html @ONLY)

