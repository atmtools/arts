add_library(rng INTERFACE)

target_include_directories(rng INTERFACE  "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(rng INTERFACE matpack artstime)
