find_program(LEGACY_MAKE make REQUIRED)

set(WIGXJPF_CPP_ERROR 1)

add_subdirectory(wigxjpf)

if(FASTWIGNER)
  add_subdirectory(fastwigxj)
endif()

if(FASTWIGNER)
  add_library(wigner ALIAS fastwigxj)
else()
  add_library(wigner ALIAS wigxjpf)
endif()

if (FASTWIGNER)
  message (STATUS "Fast Wigner calculations enabled")
else()
  message (STATUS "Fast Wigner calculations disabled (use -DFASTWIGNER=1 to enable if your compiler supports it)")
endif()
