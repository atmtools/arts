find_package (Doxygen)

if (DOXYGEN_FOUND)

add_custom_command (
  OUTPUT Doxyfile
  COMMAND ${CMAKE_COMMAND} -D "ARTS_BINARY_DIR:STRING=\"${ARTS_BINARY_DIR}\"" -D "ARTS_SOURCE_DIR:STRING=\"${ARTS_SOURCE_DIR}\"" -P ${ARTS_SOURCE_DIR}/cmake/scripts/update_doxyfile.cmake
  COMMAND ${CMAKE_COMMAND} -E touch Doxyfile
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
  )

file(GLOB ARTS_SOURCE_FILES ${ARTS_SOURCE_DIR}/src/*.h ${ARTS_SOURCE_DIR}/src/*.cc)
add_custom_command (
  OUTPUT doxified
  COMMAND ${CMAKE_COMMAND} -E echo_append "Building API Documentation..."
  COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile > doxygen.log
  #COMMAND ${CMAKE_COMMAND} -E make_directory html
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
                           ${CMAKE_CURRENT_SOURCE_DIR}/../arts-splash.png html/
  COMMAND ${CMAKE_COMMAND} -E touch doxified
  COMMAND ${CMAKE_COMMAND} -E echo "Done."
  DEPENDS Doxyfile ${ARTS_BINARY_DIR}/auto_version.txt ${ARTS_SOURCE_FILES}
)

add_custom_target (apidoc ALL DEPENDS doxified)
add_dependencies (apidoc UpdateAutoVersion BuildAutoSources)

set_directory_properties (PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "html")

install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html
         DESTINATION share/doc/arts/doxygen)

else (DOXYGEN_FOUND)

  if (NOT DOXYGEN_FOUND)
    message(STATUS "Doxygen documentation disabled due to missing doxygen")
  endif (NOT DOXYGEN_FOUND)

endif (DOXYGEN_FOUND)

