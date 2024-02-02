add_library(physics STATIC
physics_funcs.cc
wigner_functions.cc
)
target_link_libraries(physics PUBLIC matpack wigner)
target_include_directories(physics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
