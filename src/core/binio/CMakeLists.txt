add_library (binio STATIC
        bifstream.cc
        bofstream.cc
        binio.cc
        gzstream.cc 
)

target_link_libraries (binio PUBLIC util ZLIB::ZLIB)
target_include_directories(binio PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
