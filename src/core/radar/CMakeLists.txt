add_library(radar STATIC radar.cc)

target_link_libraries(radar PUBLIC jacobian path rtepack scattering)
target_include_directories(radar PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
