add_library(scattering STATIC
   scattering_species.cc
   properties.cc
   psd.cc
   sht.cc
   phase_matrix.cc
)

target_include_directories(scattering PUBLIC  "${CMAKE_CURRENT_SOURCE_DIR}" "${FFTW_INCLUDE_DIR}" "${CMAKE_SOURCE_DIR}/3rdparty/shtns")
target_link_libraries(scattering PUBLIC auto_partfun matpack artscore shtns)
