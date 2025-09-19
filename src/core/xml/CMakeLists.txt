add_library(xml_io STATIC
 xml_io_stream_core.cpp
 xml_io_base.cpp
 xml_io.cpp
)

target_link_libraries(xml_io PUBLIC binio util strings coretypes arts_enum_options time_report)
target_include_directories(xml_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
