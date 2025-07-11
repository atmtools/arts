add_library(util STATIC
  arts_omp.cc
  debug.cpp
  format_tags.cpp
  planet_data.cpp
)
target_include_directories(util PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if (OpenMP_CXX_FOUND)
  target_link_libraries (util PUBLIC OpenMP::OpenMP_CXX)
endif ()
