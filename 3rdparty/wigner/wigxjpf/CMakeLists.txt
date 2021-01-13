add_custom_target(automatic_wigxjpf
    COMMAND $(MAKE) gen/wigxjpf_auto_config.h
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

if(NOT FASTWIGNER)
  add_library(wigxjpf STATIC
              src/prime_factor.c
              src/calc.c
              src/trivial_zero.c
              src/c_wrap.c
              src/fortran_wrap.c)

  add_dependencies(wigxjpf automatic_wigxjpf)

  include_directories(inc/ cfg/ gen/ src/)

  set_target_properties(wigxjpf PROPERTIES COMPILE_FLAGS "-fPIC ${CMAKE_C_FLAGS}")
  
  target_compile_options(wigxjpf PRIVATE -Wno-shadow)
endif()
