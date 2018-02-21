########### subdirectories ###############

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/arts.css
                ${CMAKE_CURRENT_BINARY_DIR}/arts.css COPYONLY)

add_custom_target (index_html ALL DEPENDS index.html)
add_dependencies (index_html UpdateAutoVersion)

add_custom_command (
  OUTPUT index.html
  COMMAND ${CMAKE_COMMAND} -D "ARTS_BINARY_DIR:STRING=\"${CMAKE_BINARY_DIR}\"" -D "ARTS_SOURCE_DIR:STRING=\"${CMAKE_SOURCE_DIR}\"" -P ${CMAKE_SOURCE_DIR}/cmake/scripts/update_index_html.cmake
  DEPENDS ${CMAKE_BINARY_DIR}/auto_version.txt
)

install (FILES ${CMAKE_CURRENT_BINARY_DIR}/index.html
               ${CMAKE_CURRENT_BINARY_DIR}/arts.css
         DESTINATION share/doc/arts)

add_subdirectory (doxygen)
add_subdirectory (uguide)

