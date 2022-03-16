if (ARTS_CAT_DATA_DIR)
    set(DOWNLOAD_XMLDATA)  # empty
else()
    list(APPEND DOWNLOAD_XMLDATA
        planets/Earth/Fascod/tropical/tropical.CO2.xml
        planets/Earth/Fascod/tropical/tropical.H2O.xml
        planets/Earth/Fascod/tropical/tropical.N2.xml
        planets/Earth/Fascod/tropical/tropical.O2.xml
        planets/Earth/Fascod/tropical/tropical.t.xml
        planets/Earth/Fascod/tropical/tropical.z.xml
        spectroscopy/Perrin/C2H2.xml.gz
    )
endif()

if (ARTS_CAT_DATA_DIR)
    set (DOWNLOAD_CATDATA)  # empty
else()
    list(APPEND DOWNLOAD_CATDATA
        lines/N2-44.xml
        lines/O2-66.xml
        lines/H2O-161.xml
        lines/CO2-626.xml
    )
endif()

add_custom_command (
  OUTPUT arts-xml-data arts-cat-data
  COMMAND python ${ARTS_SOURCE_DIR}/tests/testdata/get_testdata.py ${CMAKE_CURRENT_BINARY_DIR} cat ${DOWNLOAD_CATDATA} &&
          python ${ARTS_SOURCE_DIR}/tests/testdata/get_testdata.py ${CMAKE_CURRENT_BINARY_DIR} xml ${DOWNLOAD_XMLDATA}
)

add_custom_target(arts_testdata_folder_download DEPENDS arts-xml-data arts-cat-data)
