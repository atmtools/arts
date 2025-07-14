
add_executable(test_core_xml test_core_xml.cpp)
target_link_libraries(test_core_xml PUBLIC xml_io matpack)
add_test(NAME "cpp.fast.core.test_core_xml" COMMAND test_core_xml)
add_dependencies(check-deps test_core_xml)
