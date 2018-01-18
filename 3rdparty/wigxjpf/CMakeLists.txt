add_custom_target(automatic_wigxjpf
    COMMAND make gen/wigxjpf_auto_config.h
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    
add_library(wigner STATIC 
            src/prime_factor.c
            src/calc.c
            src/trivial_zero.c
            src/c_wrap.c
            src/fortran_wrap.c)
            
add_dependencies(wigner automatic_wigxjpf)

include_directories(inc/ cfg/ gen/ src/)
    
set_target_properties(wigner PROPERTIES COMPILE_FLAGS "${CMAKE_C_FLAGS}")
