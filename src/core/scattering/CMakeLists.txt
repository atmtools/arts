add_library(scattering STATIC
   scattering_species.cc
   properties.cc
   psd.cc
)

target_link_libraries(scattering PUBLIC matpack artscore)
target_include_directories(scattering PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
