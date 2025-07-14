
add_library(strings STATIC mystring.cpp)
target_include_directories(strings PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(strings PUBLIC "${ARTS_SOURCE_DIR}/3rdparty/")
target_link_libraries(strings PUBLIC util)
