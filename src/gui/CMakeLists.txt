# GUI Library: For basic GUI creation
if (ENABLE_GUI)
  add_library (artsgui STATIC
              error.cc
              gui.cc
              menu.cc
              plot.cc
              propmat.cc
              )
  target_link_libraries(artsgui PUBLIC imgui implot artscore)
  target_include_directories(artsgui PUBLIC ${ARTS_SOURCE_DIR}/src
                                            ${ARTS_SOURCE_DIR}/3rdparty/gui/imgui/
                                            ${ARTS_SOURCE_DIR}/3rdparty/gui/implot/
                                            ${ARTS_SOURCE_DIR}/3rdparty/gui/imgui-filebrowser/
                                            ${ARTS_SOURCE_DIR}/3rdparty/)
  add_dependencies(artsgui auto_generated_workspace_headers)

  if (IPO_SUPPORTED)
    set_property(TARGET artsgui PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()
endif()
