add_library(artstime STATIC artstime.cc time_report.cc)

target_include_directories(artstime PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(artstime PUBLIC matpack)