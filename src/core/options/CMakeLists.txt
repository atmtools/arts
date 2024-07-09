add_library(arts_options STATIC arts_options.cc)
target_include_directories(arts_options PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(make_enums make_enums.cc)
target_link_libraries(make_enums PUBLIC arts_options util)

add_custom_command(OUTPUT enums.h enums.cpp
  COMMAND make_enums
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS make_enums
  COMMENT "Generating automatic enumeration options"
)

add_library(arts_enum_options STATIC enums.h enums.cpp)
target_include_directories(arts_enum_options PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(arts_enum_options PUBLIC util)
