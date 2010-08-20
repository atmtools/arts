find_package (Doxygen)
find_package (LATEX)

if (DOXYGEN_FOUND AND LATEX_COMPILER AND DVIPS_CONVERTER)

set (DOXYGEN_INPUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
set (DOXYGEN_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/doxified)

add_custom_target (doxyfile_in ALL DEPENDS Doxyfile)
add_dependencies (doxyfile_in UpdateAutoVersion)

add_custom_command (
  OUTPUT Doxyfile
  COMMAND ${CMAKE_COMMAND} -D "ARTS_BINARY_DIR:STRING=\"${CMAKE_BINARY_DIR}\"" -D "ARTS_SOURCE_DIR:STRING=\"${CMAKE_SOURCE_DIR}\"" -P ${CMAKE_SOURCE_DIR}/cmake/scripts/update_doxyfile.cmake
  DEPENDS ${CMAKE_BINARY_DIR}/auto_version.txt
)

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

install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html
         DESTINATION share/doc/arts/doxygen)

else (DOXYGEN_FOUND AND LATEX_COMPILER AND DVIPS_CONVERTER)

  if (NOT DOXYGEN_FOUND)
    message(STATUS "Doxygen documentation disabled due to missing doxygen")
  endif (NOT DOXYGEN_FOUND)

  if (NOT LATEX_COMPILER)
    message(STATUS "Doxygen documentation disabled due to missing latex")
  endif (NOT LATEX_COMPILER)

  if (NOT DVIPS_CONVERTER)
    message(STATUS "Doxygen documentation disabled due to missing dvips")
  endif (NOT DVIPS_CONVERTER)

endif (DOXYGEN_FOUND AND LATEX_COMPILER AND DVIPS_CONVERTER)

