add_library(xml_partfun STATIC xml_io_partfun.cc)
target_link_libraries (xml_partfun PUBLIC matpack xml_io coretypes)
target_include_directories(xml_partfun PUBLIC ${ARTS_SOURCE_DIR}/src)

add_executable(partfun_generator make_auto_partfuns.cc)
target_link_libraries (partfun_generator PRIVATE xml_partfun species)

file(GLOB_RECURSE partfun_xml_files
         CONFIGURE_DEPENDS
         ${ARTS_CAT_DATA_DIR}/partition-functions/*.xml)

set(partfun_cpp_files "")
foreach(xml_file IN LISTS partfun_xml_files)
    get_filename_component(name_we ${xml_file} NAME_WE)
    list(APPEND partfun_cpp_files ${CMAKE_CURRENT_BINARY_DIR}/${name_we}.cpp)
endforeach()

list(LENGTH partfun_xml_files partfun_xml_files_count)
message(STATUS "Found ${partfun_xml_files_count} partition function XML files")

if (partfun_xml_files_count EQUAL 0)
  message(FATAL_ERROR "No partition function XML files found in ${ARTS_CAT_DATA_DIR}/partition-functions/
You need to pass -DARTS_CAT_DATA_PATH=<path> to cmake to point to the directory containing your copy of arts-cat-data.")
endif()

# Generate auto_partfun.h for the ARTS interface
add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/auto_partfun.h ${partfun_cpp_files}
        COMMAND partfun_generator ${partfun_xml_files}
        DEPENDS partfun_generator ${partfun_xml_files}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating automatic partition function interface"
)

add_library(auto_partfun STATIC ${partfun_cpp_files} ${CMAKE_CURRENT_BINARY_DIR}/auto_partfun.h)
target_link_libraries(auto_partfun PUBLIC matpack binio xml_partfun)
target_include_directories(auto_partfun PUBLIC ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(auto_partfun PUBLIC ${ARTS_SOURCE_DIR}/src)

if (IPO_SUPPORTED)
    set_property(TARGET auto_partfun PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

add_library(partfun STATIC partfun.cc ${CMAKE_BINARY_DIR}/src/core/spec/auto_partition_functions.h)
target_link_libraries(partfun PUBLIC auto_partfun species)
