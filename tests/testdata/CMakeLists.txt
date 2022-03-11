list(APPEND DOWNLOAD_XMLDATA
    planets/Earth/Fascod/tropical/tropical.z.xml
    planets/Earth/Fascod/tropical/tropical.t.xml
)
list(APPEND DOWNLOAD_CATDATA
    lines/O2-66.xml
)

add_custom_command (
  OUTPUT arts-xml-data arts-cat-data
  COMMAND python ${ARTS_SOURCE_DIR}/tests/testdata/get_testdata.py ${CMAKE_CURRENT_BINARY_DIR} cat ${DOWNLOAD_CATDATA} &&
          python ${ARTS_SOURCE_DIR}/tests/testdata/get_testdata.py ${CMAKE_CURRENT_BINARY_DIR} xml ${DOWNLOAD_XMLDATA}
  DEPENDS pyarts
)

add_custom_target(arts_testdata_folder_download DEPENDS arts-xml-data arts-cat-data)
