
########### next target ###############

if (ENABLE_DOCSERVER)
  if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
    set_source_files_properties (daemon.c PROPERTIES
      COMPILE_FLAGS "-Wno-shadow")
  endif (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)

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
endif (ENABLE_DOCSERVER)

