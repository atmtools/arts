add_library(sensor STATIC 
  obsel.cpp
  sensor_meta_info.cpp
  frequency_bandpass_filters.cpp
  frequency_channel_selection.cpp
  frequency_range_selection.cpp
)

target_link_libraries(sensor PUBLIC matpack rtepack arts_enum_options)

target_include_directories(sensor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
