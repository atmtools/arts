set (CMAKE_MODULE_PATH "${ARTS_SOURCE_DIR}/cmake/modules")

include (ArtsVersion)

set (top_srcdir ${ARTS_SOURCE_DIR})

file (READ "${ARTS_BINARY_DIR}/auto_version.txt" VERSION)
configure_file (${ARTS_SOURCE_DIR}/doc/doxygen/Doxyfile.in
                ${ARTS_BINARY_DIR}/doc/doxygen/Doxyfile @ONLY)

