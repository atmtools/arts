include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${FFTW_INCLUDE_DIR}
  ${ARTS_SOURCE_DIR}/3rdparty/shtns
  ${ARTS_SOURCE_DIR}/3rdparty/Eigen
  ${ARTS_SOURCE_DIR}/3rdparty/Eigen/unsupported
  )

# scattering library
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
target_link_libraries(scattering PUBLIC ${SHTNS_LIBRARY}  ${FFTW_LIBRARIES})

# test executables
add_executable (test_scattering_math test_math.cpp)
add_executable (test_sht test_sht.cpp)
target_link_libraries(test_sht scattering shtns fftw3)
