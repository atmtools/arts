if (ARTS_XML_DATA_DIR MATCHES "^${ARTS_BINARY_DIR}")
  add_custom_target(
    get-testdata-xml
    COMMAND ${Python_EXECUTABLE}
            ${ARTS_SOURCE_DIR}/testdata/get_testdata.py
            ${ARTS_SOURCE_DIR}/testdata/xml_data_files.txt
    WORKING_DIRECTORY ${ARTS_BINARY_DIR}/testdata
    COMMENT "Downloading arts-xml-data"
  )
else()
  add_custom_target(get-testdata-xml
    COMMENT "arts-xml-data is user-supplied"
  )
endif()

if (ARTS_CAT_DATA_DIR MATCHES "^${ARTS_BINARY_DIR}")
  add_custom_target(
    get-testdata-cat
    COMMAND ${Python_EXECUTABLE}
            ${ARTS_SOURCE_DIR}/testdata/get_testdata.py
            ${ARTS_SOURCE_DIR}/testdata/cat_data_files.txt
    WORKING_DIRECTORY ${ARTS_BINARY_DIR}/testdata
    COMMENT "Downloading arts-cat-data"
  )
else()
  add_custom_target(get-testdata-cat
    COMMENT "arts-cat-data is user-supplied"
  )
endif()

if (ARTS_XML_DATA_DIR MATCHES "^${ARTS_BINARY_DIR}")
  add_custom_target(
    refresh-testdata-xml
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${ARTS_BINARY_DIR}/testdata/arts-xml-data
    COMMAND ${Python_EXECUTABLE}
            ${ARTS_SOURCE_DIR}/testdata/get_testdata.py
            ${ARTS_SOURCE_DIR}/testdata/xml_data_files.txt
    WORKING_DIRECTORY ${ARTS_BINARY_DIR}/testdata
    COMMENT "Refreshing arts-xml-data"
  )
else()
  add_custom_target(refresh-testdata-xml
    COMMENT "arts-xml-data is user-supplied, skipping refresh"
  )
endif()

if (ARTS_CAT_DATA_DIR MATCHES "^${ARTS_BINARY_DIR}")
  add_custom_target(
    refresh-testdata-cat
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${ARTS_BINARY_DIR}/testdata/arts-cat-data
    COMMAND ${Python_EXECUTABLE}
            ${ARTS_SOURCE_DIR}/testdata/get_testdata.py
            ${ARTS_SOURCE_DIR}/testdata/cat_data_files.txt
    WORKING_DIRECTORY ${ARTS_BINARY_DIR}/testdata
    COMMENT "Refreshing arts-cat-data"
  )
else()
  add_custom_target(refresh-testdata-cat
    COMMENT "arts-cat-data is user-supplied, skipping refresh"
  )
endif()

add_custom_target(
  get-testdata
  DEPENDS get-testdata-xml get-testdata-cat
  COMMENT "Testdata ready"
)

add_custom_target(
  refresh-testdata
  DEPENDS refresh-testdata-xml refresh-testdata-cat
  COMMENT "Testdata refreshed"
)

add_dependencies(check-deps get-testdata)
