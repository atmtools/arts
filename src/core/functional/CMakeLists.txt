add_library(functional_numeric_ternary STATIC
  functional_gravity.cpp
  functional_atm_field_interp.cpp
  functional_atm_mag_field.cpp
  functional_atm_pressure_field.cpp
  functional_numeric_ternary.cpp
  functional_atm_field.cpp
)

target_link_libraries(functional_numeric_ternary PUBLIC matpack geodesy arts_enum_options artstime igrf)
target_include_directories(functional_numeric_ternary PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(functional_numeric_unary STATIC
  functional_numeric_unary.cpp
)

target_link_libraries(functional_numeric_unary PUBLIC matpack)
target_include_directories(functional_numeric_unary PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

### META-package

add_library(functional INTERFACE)
target_link_libraries(functional INTERFACE functional_numeric_ternary functional_numeric_unary)
target_include_directories(functional_numeric_ternary PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
