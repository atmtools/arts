include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${ARTS_SOURCE_DIR}/3rdparty/shtns
  ${ARTS_SOURCE_DIR}/3rdparty/Eigen/unsupported
)
add_library(scattering STATIC sht.cxx)
add_dependencies(scattering libshtns)
target_link_libraries(scattering ${SHTNS_LIBRARY})

add_executable (test_scattering_math test_math.cc)

add_executable (test_sht test_sht.cc)
target_link_libraries(test_sht scattering ${SHTNS_LIBRARY} fftw3)
