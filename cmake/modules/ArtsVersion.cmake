macro (ARTS_GET_VERSION VERSION_NUMBER)

set (VERSION_FILE "${ARTS_SOURCE_DIR}/VERSION")
if (EXISTS ${VERSION_FILE})
  file (STRINGS "${VERSION_FILE}" ARTS_FULL_VERSION)
  string (REGEX REPLACE "([^;]+).*" "\\1"
    ARTS_FULL_VERSION "${ARTS_FULL_VERSION}")
else()
  file (READ "${ARTS_SOURCE_DIR}/ChangeLog" ARTS_CHANGELOG)
  string (REGEX MATCH "\\* arts-([0-9]+-[0-9]+-[0-9a-z]+)"
    ARTS_FULL_VERSION "${ARTS_CHANGELOG}")
  string (REGEX REPLACE ".*\\* arts-([0-9]+)-([0-9]+)-([0-9]+.*)" "\\1.\\2.\\3"
    ARTS_FULL_VERSION "${ARTS_FULL_VERSION}")
endif()

set (GIT_DIR "${ARTS_SOURCE_DIR}/.git")

if (IS_DIRECTORY "${GIT_DIR}")
  execute_process (COMMAND
                   git describe --always --abbrev=8 --dirty
                   WORKING_DIRECTORY "${ARTS_SOURCE_DIR}"
                   OUTPUT_VARIABLE GIT_COMMIT_HASH
                   ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

if (GIT_COMMIT_HASH)
  set (${VERSION_NUMBER} "${ARTS_FULL_VERSION} (git: ${GIT_COMMIT_HASH})")
else()
  set (${VERSION_NUMBER} "${ARTS_FULL_VERSION}")
endif()

endmacro ()


macro (ARTS_WRITE_VERSION_FILE VERSION_FILE VERSION_STRING)
  file (STRINGS "${ARTS_BINARY_DIR}/auto_version.txt" ARTS_VERSION)
  file (WRITE "${VERSION_FILE}" "${VERSION_STRING}")
endmacro ()


