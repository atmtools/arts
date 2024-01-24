add_library(scattering STATIC
   species.cc
)

target_link_libraries(scattering PUBLIC matpack)
target_include_directories(scattering PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/..")
