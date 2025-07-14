add_library (binio STATIC
        bifstream.cc
        bofstream.cc
        binio.cc
        gzstream.cc 
)

set_source_files_properties (binio.cc PROPERTIES
                             COMPILE_FLAGS "-fno-strict-aliasing")
target_link_libraries (binio PUBLIC util ZLIB::ZLIB)
target_include_directories(binio PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
