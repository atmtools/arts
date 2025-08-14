
add_executable(test_core_xml test_core_xml.cpp)
target_link_libraries(test_core_xml PUBLIC xml_io matpack)
add_test(NAME "cpp.fast.core.test_core_xml" COMMAND test_core_xml)
add_dependencies(check-deps test_core_xml)


add_executable(test_isosorting test_isosorting.cpp)
target_link_libraries(test_isosorting PUBLIC species)
add_test(NAME "cpp.fast.core.test_isosorting" COMMAND test_isosorting)
add_dependencies(check-deps test_isosorting)


add_executable(test_laginterp test_laginterp.cpp)
target_link_libraries(test_laginterp PUBLIC matpack rng)
add_test(NAME "cpp.fast.core.test_laginterp" COMMAND test_laginterp)
add_dependencies(check-deps test_laginterp)
