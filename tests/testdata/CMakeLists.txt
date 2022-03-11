list(APPEND DOWNLOAD_XMLDATA
    planets/Earth/Fascod/tropical/tropical.C2H2.xml
    planets/Earth/Fascod/tropical/tropical.C2H6.xml
    planets/Earth/Fascod/tropical/tropical.CH3Cl.xml
    planets/Earth/Fascod/tropical/tropical.CH4.xml
    planets/Earth/Fascod/tropical/tropical.ClO.xml
    planets/Earth/Fascod/tropical/tropical.CO.xml
    planets/Earth/Fascod/tropical/tropical.CO2.xml
    planets/Earth/Fascod/tropical/tropical.COF2.xml
    planets/Earth/Fascod/tropical/tropical.H2CO.xml
    planets/Earth/Fascod/tropical/tropical.H2O.xml
    planets/Earth/Fascod/tropical/tropical.H2O2.xml
    planets/Earth/Fascod/tropical/tropical.HBr.xml
    planets/Earth/Fascod/tropical/tropical.HCl.xml
    planets/Earth/Fascod/tropical/tropical.HCN.xml
    planets/Earth/Fascod/tropical/tropical.HF.xml
    planets/Earth/Fascod/tropical/tropical.HI.xml
    planets/Earth/Fascod/tropical/tropical.HNO3.xml
    planets/Earth/Fascod/tropical/tropical.HOCl.xml
    planets/Earth/Fascod/tropical/tropical.N2.xml
    planets/Earth/Fascod/tropical/tropical.N2O.xml
    planets/Earth/Fascod/tropical/tropical.NH3.xml
    planets/Earth/Fascod/tropical/tropical.NO.xml
    planets/Earth/Fascod/tropical/tropical.NO2.xml
    planets/Earth/Fascod/tropical/tropical.O2.xml
    planets/Earth/Fascod/tropical/tropical.O3.xml
    planets/Earth/Fascod/tropical/tropical.OH.xml
    planets/Earth/Fascod/tropical/tropical.PH3.xml
    planets/Earth/Fascod/tropical/tropical.SF6.xml
    planets/Earth/Fascod/tropical/tropical.SO2.xml
    planets/Earth/Fascod/tropical/tropical.t.xml
    planets/Earth/Fascod/tropical/tropical.z.xml
    spectroscopy/Perrin/C2H2.xml.gz
)
list(APPEND DOWNLOAD_CATDATA
    lines/O2-66.xml
    lines/O3-666.xml
)

add_custom_command (
  OUTPUT arts-xml-data arts-cat-data
  COMMAND python ${ARTS_SOURCE_DIR}/tests/testdata/get_testdata.py ${CMAKE_CURRENT_BINARY_DIR} cat ${DOWNLOAD_CATDATA} &&
          python ${ARTS_SOURCE_DIR}/tests/testdata/get_testdata.py ${CMAKE_CURRENT_BINARY_DIR} xml ${DOWNLOAD_XMLDATA}
  DEPENDS pyarts
)

add_custom_target(arts_testdata_folder_download DEPENDS arts-xml-data arts-cat-data)
