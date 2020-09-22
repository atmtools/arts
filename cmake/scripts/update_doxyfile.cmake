list (APPEND CMAKE_MODULE_PATH "${ARTS_SOURCE_DIR}/cmake/modules")

include (ArtsVersion)

set (top_srcdir ${ARTS_SOURCE_DIR})
set (top_builddir ${ARTS_BINARY_DIR})

set (HAVE_DOT NO)
if (ENABLE_DOT)
  find_program(FOUND_DOT dot)

  if (FOUND_DOT MATCHES "/dot$")
    set (HAVE_DOT YES)
  endif()
endif()

file (READ "${ARTS_BINARY_DIR}/auto_version.txt" VERSION)
configure_file (${ARTS_SOURCE_DIR}/doc/doxygen/Doxyfile.in
                ${ARTS_BINARY_DIR}/doc/doxygen/Doxyfile @ONLY)

