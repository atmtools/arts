if(NOT ARTS_CAT_DATA_DIR MATCHES "^${ARTS_BINARY_DIR}")
    set(DOWNLOAD_CATDATA "") # empty
else()
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/cat_data_files.txt" DOWNLOAD_CATDATA_RAW)
    foreach(line IN LISTS DOWNLOAD_CATDATA_RAW)
        string(STRIP "${line}" line)
        if(NOT line STREQUAL "" AND NOT line MATCHES "^#")
            list(APPEND DOWNLOAD_CATDATA "${line}")
        endif()
    endforeach()
endif()

string(REGEX REPLACE ";" " arts-cat-data/" DOWNLOAD_CATDATA_DEPS "${DOWNLOAD_CATDATA}")

if(DOWNLOAD_CATDATA_DEPS)
    string(PREPEND DOWNLOAD_CATDATA_DEPS "arts-cat-data/")
    add_custom_target(
        arts_testdata_folder_download_cat
        COMMAND ${Python_EXECUTABLE} ${ARTS_SOURCE_DIR}/testdata/get_testdata.py ${CMAKE_CURRENT_SOURCE_DIR}/cat_data_files.txt
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS get_testdata.py
        COMMENT "Downloading cat testdata"
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
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/xml_data_files.txt" DOWNLOAD_XMLDATA_RAW)
    foreach(line IN LISTS DOWNLOAD_XMLDATA_RAW)
        string(STRIP "${line}" line)
        if(NOT line STREQUAL "" AND NOT line MATCHES "^#")
            list(APPEND DOWNLOAD_XMLDATA "${line}")
        endif()
    endforeach()
endif()

string(REGEX REPLACE ";" " arts-xml-data/" DOWNLOAD_XMLDATA_DEPS "${DOWNLOAD_XMLDATA}")

if(DOWNLOAD_XMLDATA_DEPS)
    string(PREPEND DOWNLOAD_XMLDATA_DEPS "arts-xml-data/")
    add_custom_target(
        arts_testdata_folder_download_xml
        COMMAND ${Python_EXECUTABLE} ${ARTS_SOURCE_DIR}/testdata/get_testdata.py ${CMAKE_CURRENT_SOURCE_DIR}/xml_data_files.txt
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS get_testdata.py
        COMMENT "Downloading xml testdata"
    )

    add_dependencies(check-deps arts_testdata_folder_download_xml)

    #
    # Clean
    #
    set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_CLEAN_FILES arts-xml-data)
endif()
