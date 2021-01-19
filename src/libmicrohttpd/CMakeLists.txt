
########### next target ###############

if (ENABLE_DOCSERVER)
  if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    file (GLOB MICROHTTPDSOURCES "*.c")
    set_source_files_properties (${MICROHTTPDSOURCES} PROPERTIES
      COMPILE_FLAGS "-Wno-shadow -Wno-conversion")
  endif (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

  add_library (microhttpd STATIC
    connection.c
    connection.h
    daemon.c
    internal.c
    internal.h
    memorypool.c
    memorypool.h
    microhttpd.h
    platform.h
    reason_phrase.c
    reason_phrase.h
    response.c
    response.h)
    
  target_compile_options(microhttpd PRIVATE -Wno-unused-result)
endif (ENABLE_DOCSERVER)

