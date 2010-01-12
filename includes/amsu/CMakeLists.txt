
set (ARTS_INCLUDE_DIR share/arts/includes/amsu)

foreach (EXTENSION arts xml xml.gz m)
  file (GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.${EXTENSION}")
  install (FILES ${files} DESTINATION ${ARTS_INCLUDE_DIR})
endforeach (EXTENSION)

