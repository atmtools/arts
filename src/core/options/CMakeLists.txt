add_library(arts_options STATIC arts_options.cc)
target_include_directories(arts_options PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(arts_options PUBLIC strings species_enum_info_gen)

add_executable(make_enums make_enums.cc)
target_link_libraries(make_enums PUBLIC arts_options util strings)

arts_generate(make_enums
  OUTPUT enums.h enums.cpp
  COMMENT "Generating automatic enumeration options"
)

add_library(arts_enum_options STATIC enums.h enums.cpp)
target_include_directories(arts_enum_options PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(arts_enum_options PUBLIC util strings)
