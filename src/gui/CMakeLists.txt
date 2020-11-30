# GUI Library: For basic GUI creation
add_library (artsgui STATIC
            gui.cc
            plot.cc
            )
target_link_libraries(artsgui PUBLIC imgui implot matpack)
target_include_directories(artsgui PUBLIC ${ARTS_SOURCE_DIR}/src
                                          ${ARTS_SOURCE_DIR}/3rdparty/gui/imgui/
                                          ${ARTS_SOURCE_DIR}/3rdparty/gui/implot/
                                          ${ARTS_SOURCE_DIR}/3rdparty/gui/imgui-filebrowser/
                                          ${ARTS_SOURCE_DIR}/3rdparty/)
