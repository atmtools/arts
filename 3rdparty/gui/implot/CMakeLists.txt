add_library(implot STATIC
            implot.cpp
            implot_items.cpp
            )
target_include_directories(implot PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(implot PUBLIC ImDrawIdx=unsigned\ int)
target_compile_options(implot PRIVATE -Wno-conversion)
target_link_libraries(implot imgui)
