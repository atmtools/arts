add_library(surface_scattering STATIC
    lambertian.cc
    surface_scattering_model.cc
)

target_include_directories(surface_scattering PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(surface_scattering PUBLIC matpack surface)

