if(NOT ARTS_CAT_DATA_DIR MATCHES "^${ARTS_BINARY_DIR}")
    set(DOWNLOAD_CATDATA "") # empty
else()
    list(APPEND DOWNLOAD_CATDATA
        cia/O2-CIA-O2.xml
        cia/O2-CIA-N2.xml
        lines/O2-66.xml
        lines/N2-44.xml
        lines/CO2-626.xml
        lines/H2O-161.xml
        predef/H2O-SelfContCKDMT400.xml
        predef/H2O-ForeignContCKDMT400.xml
        xsec/O3-XFIT.xml
        xsec/O3-XFIT.xml.bin
    )
endif()

string(REGEX REPLACE ";" " arts-cat-data/" DOWNLOAD_CATDATA_DEPS "${DOWNLOAD_CATDATA}")

if(DOWNLOAD_CATDATA_DEPS)
    string(PREPEND DOWNLOAD_CATDATA_DEPS "arts-cat-data/")
    add_custom_target(
        arts_testdata_folder_download_cat
        COMMAND ${Python_EXECUTABLE} ${ARTS_SOURCE_DIR}/testdata/get_testdata.py cat ${DOWNLOAD_CATDATA}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS get_testdata.py
        COMMENT "Downloading test data"
    )

    add_dependencies(check-deps arts_testdata_folder_download_cat)

    #
    # Clean
    #
    set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_CLEAN_FILES arts-cat-data)
endif()

if(NOT ARTS_XML_DATA_DIR MATCHES "^${ARTS_BINARY_DIR}")
    set(DOWNLOAD_XMLDATA "") # empty
else()
    list(APPEND DOWNLOAD_XMLDATA
        planets/Earth/afgl/tropical/t.xml
        planets/Earth/afgl/tropical/p.xml
        planets/Earth/afgl/tropical/O2.xml
        planets/Earth/afgl/tropical/CO2.xml
        planets/Earth/afgl/tropical/H2O.xml
        planets/Earth/afgl/tropical/O3.xml
        planets/Earth/afgl/tropical/N2.xml
        star/Sun/solar_spectrum_QUIET.xml
        planets/Earth/Optimized-Flux-Frequencies/SW-flux-optimized-f_grid.xml
        planets/Earth/Optimized-Flux-Frequencies/LW-flux-optimized-f_grid.xml
        planets/Earth/Optimized-Flux-Frequencies/SW-flux-optimized-quadrature_weights.xml
        planets/Earth/Optimized-Flux-Frequencies/LW-flux-optimized-quadrature_weights.xml
        scattering/H2O_ice/ScatteringMetaFile_allH2Oice.xml
        scattering/H2O_ice/SingleScatteringFile_allH2Oice.xml
        scattering/H2O_ice/MieSphere_R5.00000e-01um.meta.xml
        scattering/H2O_ice/MieSphere_R5.00000e-01um.xml
    )
endif()

string(REGEX REPLACE ";" " arts-xml-data/" DOWNLOAD_XMLDATA_DEPS "${DOWNLOAD_XMLDATA}")

if(DOWNLOAD_XMLDATA_DEPS)
    string(PREPEND DOWNLOAD_XMLDATA_DEPS "arts-xml-data/")
    add_custom_target(
        arts_testdata_folder_download_xml
        COMMAND ${Python_EXECUTABLE} ${ARTS_SOURCE_DIR}/testdata/get_testdata.py xml ${DOWNLOAD_XMLDATA}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS get_testdata.py
    )

    add_dependencies(check-deps arts_testdata_folder_download_xml)

    #
    # Clean
    #
    set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_CLEAN_FILES arts-xml-data)
endif()
