find_package (Doxygen)

if (DOXYGEN_FOUND)

set (DOXYGEN_INPUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
set (DOXYGEN_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/doxified)

set (top_srcdir ${CMAKE_SOURCE_DIR})
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
                ${DOXYGEN_INPUT}
                @ONLY IMMEDIATE)

add_custom_command (
  OUTPUT ${DOXYGEN_OUTPUT}
  COMMAND ${CMAKE_COMMAND} -E echo_append "Building API Documentation..."
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT} > doxygen.log
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
                           ${CMAKE_CURRENT_SOURCE_DIR}/../arts-splash.gif html/
  COMMAND ${CMAKE_COMMAND} -E touch doxified
  COMMAND ${CMAKE_COMMAND} -E echo "Done."
  DEPENDS ${DOXYGEN_INPUT}
)

add_custom_target (apidoc ALL DEPENDS ${DOXYGEN_OUTPUT})
add_dependencies (apidoc BuildAutoSources)

add_custom_target (apidoc_forced
  COMMAND ${CMAKE_COMMAND} -E echo_append "Building API Documentation..."
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
                           ${CMAKE_CURRENT_SOURCE_DIR}/../arts-splash.gif html/
  COMMAND ${CMAKE_COMMAND} -E touch doxified
  COMMAND ${CMAKE_COMMAND} -E echo "Done."
  DEPENDS ${DOXYGEN_INPUT}
)

set_directory_properties (PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "html")

endif (DOXYGEN_FOUND)

