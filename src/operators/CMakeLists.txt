add_library(operators STATIC operators.cc)
target_link_libraries(operators PUBLIC matpack)
target_include_directories(operators PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(make_array_ostream make_array_ostream.cc)

add_custom_command(
  OUTPUT array_ostream.cc array_ostream.h
  COMMAND make_array_ostream > array_ostream.h 2> array_ostream.cc
  DEPENDS make_array_ostream 
  COMMENT "Generating array ostream"
  VERBATIM
)

add_library(array_ostream STATIC array_ostream.cc)
target_include_directories(array_ostream PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(array_ostream PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/../workspace)
target_link_libraries(array_ostream PUBLIC workspace_internals)

