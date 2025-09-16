add_library(operators STATIC spectral_radiance_transform_operator.cc)
target_link_libraries(operators PUBLIC rtepack path)
target_include_directories(operators PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
