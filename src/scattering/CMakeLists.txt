include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${ARTS_SOURCE_DIR}/3rdparty/shtns
  ${ARTS_SOURCE_DIR}/3rdparty/Eigen/unsupported
)
add_library(
  scattering
  STATIC
  sht.cpp
  single_scattering_data.cpp
  arts_ssdb.cpp
  particle.cpp
  particle_habit.cpp
)
add_dependencies(scattering shtns)
target_link_libraries(scattering shtns)

add_executable (test_scattering_math test_math.cpp)

add_executable (test_sht test_sht.cpp)
target_link_libraries(test_sht scattering shtns fftw3)
