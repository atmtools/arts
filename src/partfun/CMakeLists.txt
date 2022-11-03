add_executable(partfun_generator
    make_auto_partfuns.cc
    xml_io_partfun.cc
)

target_link_libraries (partfun_generator PRIVATE matpack xmliobase)
target_include_directories(partfun_generator PRIVATE ..)

file(GLOB_RECURSE partfun_xml_files LIST_DIRECTORIES FALSE RELATIVE ${ARTS_PARTFUN_DATA_PATH} "*.xml")
set(partfun_cpp_files ${partfun_xml_files})
list(TRANSFORM partfun_cpp_files REPLACE ".xml" ".cc")
list(TRANSFORM partfun_xml_files PREPEND "${ARTS_PARTFUN_DATA_PATH}/")

# Generate auto_partfun.h for the ARTS interface
add_custom_command(OUTPUT auto_partfun.h ${partfun_cpp_files}
                   COMMAND partfun_generator ${partfun_xml_files}
                   DEPENDS partfun_generator)

add_library(auto_partfun STATIC ${partfun_cpp_files} xml_io_partfun.cc)
target_include_directories(auto_partfun PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

if (IPO_SUPPORTED)
    set_property(TARGET auto_partfun PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()