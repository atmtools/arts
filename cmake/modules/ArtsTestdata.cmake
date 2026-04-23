
macro (SETUP_ARTS_XML_DATA)
  if (ARTS_XML_DATA_PATH)
    get_filename_component (ARTS_XML_DATA_DIR "${ARTS_XML_DATA_PATH}" ABSOLUTE)
    if (NOT IS_DIRECTORY "${ARTS_XML_DATA_DIR}")
      message(FATAL_ERROR "ARTS_XML_DATA_PATH is not a valid directory: ${ARTS_XML_DATA_PATH}")
    endif()
    message(STATUS "ARTS_XML_DATA_PATH is set to ${ARTS_XML_DATA_DIR}")
  else ()
    get_filename_component (ARTS_XML_DATA_DIR "${CMAKE_SOURCE_DIR}/../arts-xml-data" ABSOLUTE)
    if (NOT IS_DIRECTORY "${ARTS_XML_DATA_DIR}")
      set (ARTS_XML_DATA_DIR "${ARTS_BINARY_DIR}/testdata/arts-xml-data")
      message (STATUS "arts-xml-data not found, will be partially downloaded to ${ARTS_XML_DATA_DIR}")

      file(MAKE_DIRECTORY ${ARTS_BINARY_DIR}/testdata)

      execute_process(
          COMMAND ${Python_EXECUTABLE}
                  ${ARTS_SOURCE_DIR}/testdata/get_testdata.py
                  ${ARTS_SOURCE_DIR}/testdata/xml_data_files.txt
          WORKING_DIRECTORY ${ARTS_BINARY_DIR}/testdata
          RESULT_VARIABLE XMLDATA_DOWNLOAD_RESULT
      )
      if (XMLDATA_DOWNLOAD_RESULT)
        message(FATAL_ERROR "Failed to download arts-xml-data.")
      endif()
    else()
      message(STATUS "Found arts-xml-data in ${ARTS_XML_DATA_DIR}")
    endif()
  endif ()
endmacro ()

macro (SETUP_ARTS_CAT_DATA)
  if (ARTS_CAT_DATA_PATH)
    get_filename_component (ARTS_CAT_DATA_DIR "${ARTS_CAT_DATA_PATH}" ABSOLUTE)
    if (NOT IS_DIRECTORY "${ARTS_CAT_DATA_DIR}")
      message(FATAL_ERROR "ARTS_CAT_DATA_PATH is not a valid directory: ${ARTS_CAT_DATA_PATH}")
    endif()
    message(STATUS "ARTS_CAT_DATA_PATH is set to ${ARTS_CAT_DATA_DIR}")
  else ()
    get_filename_component (ARTS_CAT_DATA_DIR "${CMAKE_SOURCE_DIR}/../arts-cat-data" ABSOLUTE)
    if (NOT IS_DIRECTORY "${ARTS_CAT_DATA_DIR}")
      set (ARTS_CAT_DATA_DIR "${ARTS_BINARY_DIR}/testdata/arts-cat-data")
      message (STATUS "arts-cat-data not found, will be partially downloaded to ${ARTS_CAT_DATA_DIR}")

      file(MAKE_DIRECTORY ${ARTS_BINARY_DIR}/testdata)

      execute_process(
          COMMAND ${Python_EXECUTABLE}
                  ${ARTS_SOURCE_DIR}/testdata/get_testdata.py
                  ${ARTS_SOURCE_DIR}/testdata/cat_data_files.txt
          WORKING_DIRECTORY ${ARTS_BINARY_DIR}/testdata
          RESULT_VARIABLE CATDATA_DOWNLOAD_RESULT
      )
      if (CATDATA_DOWNLOAD_RESULT)
        message(FATAL_ERROR "Failed to download arts-cat-data.")
      endif()

    else()
      message(STATUS "Found arts-cat-data in ${ARTS_CAT_DATA_DIR}")
    endif()
  endif ()
endmacro ()


