if (ARTS_CAT_DATA_DIR)
    set (DOWNLOAD_CATDATA "")  # empty
else()
    list(APPEND DOWNLOAD_CATDATA
        lines/O2-66.xml
        model/mt_ckd_4.0/H2O.xml
    )
endif()

string(REGEX REPLACE ";" " arts-cat-data/" DOWNLOAD_CATDATA_DEPS "${DOWNLOAD_CATDATA}")

if(DOWNLOAD_CATDATA_DEPS)
    string(PREPEND DOWNLOAD_CATDATA_DEPS "arts-cat-data/")
    add_custom_command(
        OUTPUT ${DOWNLOAD_CATDATA_DEPS}
        COMMAND ${Python3_EXECUTABLE} ${ARTS_SOURCE_DIR}/tests/testdata/get_testdata.py ${CMAKE_CURRENT_BINARY_DIR} cat ${DOWNLOAD_CATDATA}
        DEPENDS get_testdata.py
    )

    add_custom_target(arts_testdata_folder_download_cat DEPENDS ${DOWNLOAD_CATDATA_DEPS})

    add_dependencies(check-deps arts_testdata_folder_download_cat)

    #
    # Clean
    #
    set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES arts-cat-data)
endif()


if (ARTS_XML_DATA_DIR)
    set(DOWNLOAD_XMLDATA "")  # empty
else()
    list(APPEND DOWNLOAD_XMLDATA
        planets/Earth/Fascod/tropical/tropical.t.xml
        planets/Earth/Fascod/tropical/tropical.z.xml
        planets/Earth/Fascod/tropical/tropical.O2.xml
        planets/Earth/Fascod/tropical/tropical.N2.xml
        star/Sun/solar_spectrum_May_2004.xml
        scattering/H2O_ice
        scattering/H2O_frozen
    )
endif()

string(REGEX REPLACE ";" " arts-xml-data/" DOWNLOAD_XMLDATA_DEPS "${DOWNLOAD_XMLDATA}")

if(DOWNLOAD_XMLDATA_DEPS)
    string(PREPEND DOWNLOAD_XMLDATA_DEPS "arts-xml-data/")
    add_custom_command(
        OUTPUT ${DOWNLOAD_XMLDATA_DEPS}
        COMMAND ${Python3_EXECUTABLE} ${ARTS_SOURCE_DIR}/tests/testdata/get_testdata.py ${CMAKE_CURRENT_BINARY_DIR} xml ${DOWNLOAD_XMLDATA}
        DEPENDS get_testdata.py
    )

    add_custom_target(arts_testdata_folder_download_xml DEPENDS ${DOWNLOAD_XMLDATA_DEPS})

    add_dependencies(check-deps arts_testdata_folder_download_xml)

    #
    # Clean
    #
    set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES arts-xml-data)



endif()



list(APPEND DOWNLOAD_SSDBDATA
  Ice/SingleCrystals/Pristine/ThinPlate_Id16/TotallyRandom/Dmax00025um_Dveq00010um_Mass4.79983e-13kg.nc
  Ice/SingleCrystals/Pristine/ThinPlate_Id16/TotallyRandom/Dmax00367um_Dveq00145um_Mass1.46402e-09kg.nc
  Liquid/SingleCrystals/Pristine/LiquidSphere_Id25/TotallyRandom/Dmax00010um_Dveq00010um_Mass4.94568e-13kg.nc Liquid/SingleCrystals/Pristine/LiquidSphere_Id25/TotallyRandom/Dmax00023um_Dveq00023um_Mass6.72428e-12kg.nc
  Liquid/SingleCrystals/Pristine/LiquidSphere_Id25/TotallyRandom/Dmax01459um_Dveq01459um_Mass1.62701e-06kg.nc
)
string(REGEX REPLACE ";" " " DOWNLOAD_SSDBDATA_DEPS "${DOWNLOAD_SSDBDATA}")

if(DOWNLOAD_SSDBDATA_DEPS)
  add_custom_command(
    OUTPUT ${DOWNLOAD_SSDBDATA_DEPS}
    COMMAND ${Python3_EXECUTABLE} ${ARTS_SOURCE_DIR}/tests/testdata/get_testdata.py ${CMAKE_CURRENT_BINARY_DIR} ssdb ${DOWNLOAD_SSDBDATA}
    DEPENDS get_testdata.py
    )
  add_custom_target(arts_testdata_folder_download_ssdb DEPENDS ${DOWNLOAD_SSDBDATA_DEPS})
endif()
