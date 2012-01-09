if (QT4_FOUND)
  set (artsgui_SOURCES ag_main.cc mainwindow_qt.cc)
  set (artsgui_HEADERS mainwindow_qt.h)

  file (GLOB HEADERFILES "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
  add_custom_target(UtilityHeadersGui SOURCES ${HEADERFILES})

  qt4_wrap_cpp (artsgui_HEADERS_MOC ${artsgui_HEADERS})

  include (${QT_USE_FILE})
  add_definitions (${QT_DEFINITIONS})

  add_library (artsgui ${artsgui_SOURCES} ${artsgui_HEADERS_MOC})
endif (QT4_FOUND)

