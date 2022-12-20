include_directories (
  ${FFTW_INCLUDE_DIR}
  ${ARTS_BINARY_DIR}/src
  ${ARTS_SOURCE_DIR}/src
  ${ARTS_SOURCE_DIR}/3rdparty/shtns
  ${ARTS_SOURCE_DIR}/3rdparty/Eigen
  ${ARTS_SOURCE_DIR}/3rdparty/Eigen/unsupported
  )

# scattering library
add_library(
  scattering
  STATIC
  sht.cpp
  scattering_data_field.cpp
  single_scattering_data.cpp
  arts_ssdb.cpp
  particle.cpp
  particle_habit.cpp
  bulk_particle_habit.cpp
)
add_dependencies(scattering shtns BuildAutoSources)
target_link_libraries(scattering PUBLIC shtns  ${FFTW_LIBRARIES})

# test executables
add_executable (test_scattering_math test_math.cpp)

add_executable (test_scattering_interpolation test_interpolation.cpp)

add_executable (test_sht test_sht.cpp)

add_executable (test_mie test_mie.cpp)

add_executable (test_scattering_data_field test_scattering_data_field.cpp)
target_link_libraries(test_scattering_data_field PUBLIC scattering shtns  ${FFTW_LIBRARIES})

add_executable(test_single_scattering_data test_single_scattering_data.cpp)
target_link_libraries(test_single_scattering_data PUBLIC scattering shtns  ${FFTW_LIBRARIES})

add_executable(test_scattering_particle test_particle.cpp)
target_link_libraries(test_scattering_particle PUBLIC scattering shtns  ${FFTW_LIBRARIES} ${NETCDF_LIBRARIES})

add_executable(test_scattering_particle_habit test_particle_habit.cpp)
target_link_libraries(test_scattering_particle_habit PUBLIC scattering shtns  ${FFTW_LIBRARIES} ${NETCDF_LIBRARIES})

target_link_libraries(test_sht scattering shtns ${FFTW_LIBRARIES})
