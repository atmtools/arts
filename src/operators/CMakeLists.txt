add_library(operators STATIC operators.cc)
target_link_libraries(operators PUBLIC matpack)
target_include_directories(operators PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
