# The licence texts of the bundled code of this build, compiled in.  This
# directory is processed after 3rdparty/, where the components register.
arts_write_bundled_license_source(${CMAKE_CURRENT_BINARY_DIR}/bundled_licenses.cpp)
add_library(arts_licenses STATIC ${CMAKE_CURRENT_BINARY_DIR}/bundled_licenses.cpp)
target_include_directories(arts_licenses PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
