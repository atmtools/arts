add_library(atm STATIC
  atm_field.cpp
  atm_mag_field.cpp
  atm_field_interp.cpp
  xml_atm.cpp
  atm_pressure_field.cpp
)

target_link_libraries(atm PUBLIC matpack operators quantum species scattering physics)
target_include_directories(atm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})