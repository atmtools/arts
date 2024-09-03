# ########## testcases ###############
add_executable(test_phase_matrix test_phase_matrix.cc)

target_link_libraries(test_phase_matrix scattering shtns fftw3)
target_include_directories(test_phase_matrix PRIVATE ${ARTS_SOURCE_DIR}/src/core)
