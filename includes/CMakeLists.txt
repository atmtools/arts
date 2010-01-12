
set (ARTS_INCLUDE_DIR share/arts/includes)

foreach (EXTENSION arts xml xml.gz m)
  file (GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.${EXTENSION}")
  install (FILES ${files} DESTINATION ${ARTS_INCLUDE_DIR})
endforeach (EXTENSION)

########### subdirectories ###############

add_subdirectory (amsu)
add_subdirectory (hirs)
add_subdirectory (mhs)
add_subdirectory (odin)

