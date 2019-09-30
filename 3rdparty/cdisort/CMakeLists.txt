add_library (cdisort STATIC
  cdisort.c
  locate.c
)

add_executable(disotest
  cdisort.c
  locate.c
  disotest.c)

set_target_properties (disotest PROPERTIES
  COMPILE_FLAGS "-Wno-all -Wno-unused-parameter -Wno-conversion -Wno-shadow -Wno-float-conversion -DENABLE_ORIGINAL_OUTPUT_HANDLING=1")

set_target_properties (cdisort PROPERTIES
  COMPILE_FLAGS "-Wno-all -Wno-unused-parameter -Wno-conversion -Wno-shadow -Wno-float-conversion")
