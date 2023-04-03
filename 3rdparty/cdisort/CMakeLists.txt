if(NOT ENABLE_ARTS_LGPL)
add_library (cdisort STATIC
  cdisort.c
  locate.c
)

add_executable(disotest
  cdisort.c
  locate.c
  disotest.c)

  target_link_libraries(disotest m)
  
  set(CDISORT_COMPILE_FLAGS "-Wno-all -Wno-unused-parameter -Wno-conversion -Wno-shadow -Wno-float-conversion")
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set(CDISORT_COMPILE_FLAGS "${CDISORT_COMPILE_FLAGS} -Wno-format-overflow")
  endif()
  
  set_target_properties (disotest PROPERTIES
      COMPILE_FLAGS "${CDISORT_COMPILE_FLAGS} -DENABLE_ORIGINAL_OUTPUT_HANDLING=1")
  
  set_target_properties (cdisort PROPERTIES
      COMPILE_FLAGS ${CDISORT_COMPILE_FLAGS})
else()
  add_library(cdisort INTERFACE)
endif()

