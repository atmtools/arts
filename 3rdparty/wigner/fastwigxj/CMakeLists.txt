if(FASTWIGNER LESS 0)
  message (SEND_ERROR "ERROR:  FASTWIGNER must be positive if active")
  return()
endif()

add_custom_target(automatic_fastwigxj
    COMMAND make gen/fastwigxj_auto_config.h
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    
include_directories(inc/ cfg/ gen/ src/ ../wigxjpf/inc/ ../wigxjpf/cfg/ ../wigxjpf/gen/ ../wigxjpf/src/)

add_definitions("-DFASTWIGXJ_USE_FLOAT128=0")

add_library(fastwigxj STATIC
            ../wigxjpf/src/prime_factor.c
            ../wigxjpf/src/calc.c
            ../wigxjpf/src/trivial_zero.c
            ../wigxjpf/src/c_wrap.c
            ../wigxjpf/src/fortran_wrap.c
            ../wigxjpf/src/calc_float128.c
            ../wigxjpf/src/c_wrap_float128.c
            src/fastwigxj.c
            src/canonicalise.c
            src/fortran_wrap.c)
            
add_dependencies(fastwigxj automatic_fastwigxj automatic_wigxjpf)
    
set_target_properties(fastwigxj PROPERTIES COMPILE_FLAGS "-fPIC ${CMAKE_C_FLAGS} -march=native")

add_library(hashing_fastwigxj STATIC
            src/unique_js.cc
            src/combine_js.cc
            src/gen_header.cc)

target_link_libraries(hashing_fastwigxj fastwigxj)

set_target_properties(hashing_fastwigxj PROPERTIES COMPILE_FLAGS "-fPIC ${CMAKE_CXX_FLAGS} -march=native")

add_executable(hash_fastwigxj src/hash_js.cc)

target_link_libraries(hash_fastwigxj hashing_fastwigxj fastwigxj)

set_target_properties(hash_fastwigxj PROPERTIES COMPILE_FLAGS "-fPIC ${CMAKE_CXX_FLAGS} -march=native")

add_custom_command(TARGET hash_fastwigxj POST_BUILD
  COMMAND "hash_fastwigxj" 
  ARGS "--max-E-3j=${FASTWIGNER}" "/dev/null" "./fasttable.3j"
  COMMAND "hash_fastwigxj" 
  ARGS "--max-E-6j=${FASTWIGNER}" "/dev/null" "./fasttable.6j"
  WORKING_DIRECTORY "${ARTS_BINARY_DIR}/3rdparty/wigner/fastwigxj/")
