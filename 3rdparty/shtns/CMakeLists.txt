file(GLOB_RECURSE shtns_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *)
foreach(f IN LISTS shtns_files)
  configure_file(${f} ${f} COPYONLY)
endforeach()

set(SHTNS_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/libshtns.a CACHE FILEPATH INTERNAL FORCE)
add_custom_command(
  OUTPUT libshtns.a
  COMMAND ./configure CFLAGS=-fpic ENABLE_MEM=no
  COMMAND make
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
add_custom_target(libshtns DEPENDS ${SHTNS_LIBRARY})
