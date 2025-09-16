add_library(scattering STATIC
   integration.cc
   scattering_species.cc
   general_tro_spectral.cc
   properties.cc
   psd.cc
   sht.cc
   utils.cc
   phase_matrix.cc
   henyey_greenstein.cc
   particle_habit.cc
   scattering_habit.cc
   bulk_scattering_properties.cc
   optproperties.cc
   gas_scattering.cc
)

target_include_directories(scattering PUBLIC  "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(scattering PUBLIC auto_partfun matpack atm rtepack)

if (DISABLE_SHTNS)
target_compile_definitions(scattering PUBLIC ARTS_NO_SHTNS)
else()
target_include_directories(scattering PRIVATE "${FFTW_INCLUDE_DIR}" "${CMAKE_SOURCE_DIR}/3rdparty/shtns")
target_link_libraries(scattering PRIVATE fftw3 shtns)
endif()
