if (ARTS_CAT_DATA_DIR)
    set (DOWNLOAD_CATDATA "")  # empty
else()
    list(APPEND DOWNLOAD_CATDATA
        lines/O2-66.xml
    )
endif()

set (DOWNLOAD_CATDATA_DEPS "")
list(LENGTH DOWNLOAD_CATDATA NCAT)
list(LENGTH DOWNLOAD_CATDATA_DEPS ICAT)
while (ICAT LESS NCAT)
    list(GET DOWNLOAD_CATDATA ${ICAT} DOWNLOAD_CATDATA_ITEM)
    list(APPEND DOWNLOAD_CATDATA_DEPS "${CMAKE_CURRENT_BINARY_DIR}/arts-cat-data/${DOWNLOAD_CATDATA_ITEM}")
    math(EXPR ICAT "${ICAT}+1")
endwhile()

add_custom_command (
  OUTPUT DOWNLOAD_CATDATA_DEPS
  COMMAND python ${ARTS_SOURCE_DIR}/tests/testdata/get_testdata.py ${CMAKE_CURRENT_BINARY_DIR} cat ${DOWNLOAD_CATDATA}
  DEPENDS get_testdata.py
)

if (ARTS_XML_DATA_DIR)
    set(DOWNLOAD_XMLDATA "")  # empty
else()
    list(APPEND DOWNLOAD_XMLDATA
        planets/Earth/Fascod/tropical/tropical.O2.xml
    )
endif()

set (DOWNLOAD_XMLDATA_DEPS "")
list(LENGTH DOWNLOAD_XMLDATA NXML)
list(LENGTH DOWNLOAD_XMLDATA_DEPS IXML)
while (IXML LESS NXML)
    list(GET DOWNLOAD_XMLDATA ${IXML} DOWNLOAD_XMLDATA_ITEM)
    list(APPEND DOWNLOAD_XMLDATA_DEPS "${CMAKE_CURRENT_BINARY_DIR}/arts-xml-data/${DOWNLOAD_XMLDATA_ITEM}")
    math(EXPR IXML "${IXML}+1")
endwhile()

add_custom_command (
  OUTPUT DOWNLOAD_XMLDATA_DEPS
  COMMAND python ${ARTS_SOURCE_DIR}/tests/testdata/get_testdata.py ${CMAKE_CURRENT_BINARY_DIR} xml ${DOWNLOAD_XMLDATA}
  DEPENDS get_testdata.py
)

add_custom_target(arts_testdata_folder_download DEPENDS DOWNLOAD_CATDATA_DEPS DOWNLOAD_XMLDATA_DEPS)

add_dependencies(check-deps arts_testdata_folder_download)
