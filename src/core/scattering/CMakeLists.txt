add_library(scattering STATIC
   species.cc
   properties.cc
)

target_link_libraries(scattering PUBLIC matpack artscore)
target_include_directories(scattering PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/..")
