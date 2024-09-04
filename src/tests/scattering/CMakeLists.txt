# ########## testcases ###############
add_executable(test_phase_matrix test_phase_matrix.cc)
target_link_libraries(test_phase_matrix scattering shtns fftw3)
target_include_directories(test_phase_matrix PRIVATE ${ARTS_SOURCE_DIR}/src/core)
add_test(NAME "cpp.fast.scattering.test_phase_matrix" COMMAND test_phase_matrix)

add_executable(test_extinction_matrix test_extinction_matrix.cc)
target_link_libraries(test_extinction_matrix scattering shtns fftw3)
target_include_directories(test_extinction_matrix PRIVATE ${ARTS_SOURCE_DIR}/src/core)
add_test(NAME "cpp.fast.scattering.test_extinction_matrix" COMMAND test_extinction_matrix)

add_executable(test_absorption_vector test_absorption_vector.cc)
target_link_libraries(test_absorption_vector scattering shtns fftw3)
target_include_directories(test_absorption_vector PRIVATE ${ARTS_SOURCE_DIR}/src/core)
add_test(NAME "cpp.fast.scattering.test_absorption_vector" COMMAND test_absorption_vector)

add_executable(test_sht test_sht.cc)
target_link_libraries(test_sht scattering shtns fftw3)
target_include_directories(test_sht PRIVATE ${ARTS_SOURCE_DIR}/src/core)
add_test(NAME "cpp.fast.scattering.test_sht" COMMAND test_sht)

add_executable(test_angle_integration test_integration.cc)
target_link_libraries(test_angle_integration scattering shtns fftw3)
target_include_directories(test_angle_integration PRIVATE ${ARTS_SOURCE_DIR}/src/core)
add_test(NAME "cpp.fast.scattering.angle_integration" COMMAND test_angle_integration)

add_executable(test_single_scattering_data test_single_scattering_data.cc)
target_link_libraries(test_single_scattering_data scattering shtns fftw3 artsworkspace)
target_include_directories(test_single_scattering_data PRIVATE ${ARTS_SOURCE_DIR}/src/core)
add_dependencies(test_single_scattering_data arts_testdata_folder_download_xml)
add_test(NAME "cpp.fast.scattering.test_single_scattering_data" COMMAND test_single_scattering_data)
