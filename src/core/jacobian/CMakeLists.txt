add_library(jacobian STATIC
jacobian.cc
)

target_link_libraries(jacobian PUBLIC atm sensor scattering lbl)
target_include_directories(jacobian PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
