set (ENABLE_CDISORT true)

add_library (cdisort STATIC
  cdisort.c
  locate.c
)

add_executable(disotest
  disotest.c)
target_link_libraries(disotest cdisort)

set_target_properties (disotest PROPERTIES
  COMPILE_FLAGS "-Wno-all -Wno-unused-parameter -Wno-conversion -Wno-shadow -Wno-float-conversion")

set_target_properties (cdisort PROPERTIES
  COMPILE_FLAGS "-Wno-all -Wno-unused-parameter -Wno-conversion -Wno-shadow -Wno-float-conversion")
