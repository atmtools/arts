add_library(operators INTERFACE)
target_link_libraries(operators INTERFACE matpack)
target_include_directories(operators INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
