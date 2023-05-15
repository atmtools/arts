# GUI Library: For basic GUI creation
if (ENABLE_GUI)
  add_library (artsgui STATIC
              error.cc
              gui.cc
              menu.cc
              plot.cc
              propmat.cc
              )
  add_dependencies(artsgui BuildAutoSources)
  target_link_libraries(artsgui PUBLIC imgui implot matpack artscore species)
  target_include_directories(artsgui PUBLIC ${ARTS_SOURCE_DIR}/src
                                            ${ARTS_SOURCE_DIR}/3rdparty/gui/imgui/
                                            ${ARTS_SOURCE_DIR}/3rdparty/gui/implot/
                                            ${ARTS_SOURCE_DIR}/3rdparty/gui/imgui-filebrowser/
                                            ${ARTS_SOURCE_DIR}/3rdparty/)

  if (IPO_SUPPORTED)
    set_property(TARGET artsgui PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()
endif()