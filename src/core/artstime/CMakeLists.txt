add_library(artstime STATIC artstime.cc)

target_include_directories(artstime PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(artstime PUBLIC time_report xml_io FastFloat::fast_float)

####################################################################################

add_library(time_report STATIC time_report.cc)

target_include_directories(time_report PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(time_report PUBLIC util strings)
