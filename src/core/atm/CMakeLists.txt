add_library(atm STATIC
  atm_field.cpp
  xml_atm.cpp
)

target_link_libraries(atm PUBLIC matpack operators functional_numeric_ternary quantum species scattering physics)
target_include_directories(atm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
