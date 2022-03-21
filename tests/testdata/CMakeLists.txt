if (ARTS_XML_DATA_DIR)
    set(DOWNLOAD_XMLDATA "")  # empty
else()
    list(APPEND DOWNLOAD_XMLDATA
        planets/Earth/Fascod/tropical/tropical.CO2.xml
        planets/Earth/Fascod/tropical/tropical.H2O.xml
        planets/Earth/Fascod/tropical/tropical.N2.xml
        planets/Earth/Fascod/tropical/tropical.O2.xml
        planets/Earth/Fascod/tropical/tropical.O3.xml
        planets/Earth/Fascod/tropical/tropical.t.xml
        planets/Earth/Fascod/tropical/tropical.z.xml
        spectroscopy/Perrin/C2H2.xml.gz
        spectroscopy/Perrin/He.xml.gz
        spectroscopy/cia/hitran2011/hitran_cia2012_adapted.xml.gz
    )
endif()

if (ARTS_CAT_DATA_DIR)
    set (DOWNLOAD_CATDATA "")  # empty
else()
    list(APPEND DOWNLOAD_CATDATA
        lines/N2-44.xml
        lines/O2-66.xml
        lines/O3-666.xml
        lines/H2O-161.xml
        lines/CO2-626.xml
    )
endif()

add_custom_command (
  OUTPUT arts_testdata_folders
  COMMAND python ${ARTS_SOURCE_DIR}/tests/testdata/get_testdata.py ${CMAKE_CURRENT_BINARY_DIR} cat ${DOWNLOAD_CATDATA} &&
          python ${ARTS_SOURCE_DIR}/tests/testdata/get_testdata.py ${CMAKE_CURRENT_BINARY_DIR} xml ${DOWNLOAD_XMLDATA}
  DEPENDS get_testdata.py
)


add_custom_command (
  OUTPUT abs_lines.xml
  COMMAND python ${CMAKE_CURRENT_BINARY_DIR} full ${ARTS_CAT_DATA_DIR}/lines/O2-66.xml ${CMAKE_CURRENT_BINARY_DIR}/arts-cat-data/lines/O2-66.xml 60e9 61e9 0 none
  DEPENDS prep_testdata.py pyarts
  )

add_custom_target(arts_testdata_folder_download DEPENDS arts_testdata_folders)

add_custom_target(prepare_abs_lines DEPENDS arts_testdata_folders abs_lines.xml)

add_dependencies(check-deps arts_testdata_folder_download)
