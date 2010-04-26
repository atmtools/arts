########### subdirectories ###############

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/index.html.in
                ${CMAKE_CURRENT_BINARY_DIR}/index.html @ONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/arts.css
                ${CMAKE_CURRENT_BINARY_DIR}/arts.css COPYONLY)

install (FILES ${CMAKE_CURRENT_BINARY_DIR}/index.html
               ${CMAKE_CURRENT_BINARY_DIR}/arts.css
         DESTINATION share/doc/arts)

add_subdirectory (doxygen)
add_subdirectory (uguide)


