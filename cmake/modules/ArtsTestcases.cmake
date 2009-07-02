macro (ARTS_ADD_PREPARE_TESTCASES TARGETS)
  set (ADD_TARGETS "")
  foreach(TESTCASE ${ARGN})
    set (TARGET PrepareTest${TESTCASE})
    file (GLOB COPYFILES
      "${CMAKE_CURRENT_SOURCE_DIR}/${TESTCASE}/*.arts"
      "${CMAKE_CURRENT_SOURCE_DIR}/${TESTCASE}/*.xml"
      "${CMAKE_CURRENT_SOURCE_DIR}/${TESTCASE}/*.xml.bin"
      "${CMAKE_CURRENT_SOURCE_DIR}/${TESTCASE}/*.xml.gz"
    )
    add_custom_command(
      OUTPUT ${TARGET}
      COMMAND cmake -E
        ARGS make_directory ${CMAKE_CURRENT_BINARY_DIR}/${TESTCASE}
      COMMAND cp
        ARGS -u ${COPYFILES}
             "${CMAKE_CURRENT_BINARY_DIR}/${TESTCASE}/"
    )
    set (ADD_TARGETS ${ADD_TARGETS} ${TARGET})
  endforeach (TESTCASE)
  set(${TARGETS} ${ADD_TARGETS})
endmacro (ARTS_ADD_PREPARE_TESTCASES)

macro (ARTS_ADD_TESTCASES TARGETS DEPENDENCIES)
  set (ADD_TARGETS "")
  foreach(TESTCASE ${ARGN})
    set (TARGET Test${TESTCASE})
    add_custom_command(
      OUTPUT ${TARGET}
      COMMAND ${CMAKE_TEST_COMMAND} if TOPSRCDIR=${CMAKE_SOURCE_DIR} python
              testall.py ${TARGET} > ${TARGET}.log 2>&1\; then
              echo "${TARGET} ok."\; else
              echo "${TARGET} failed. Find details in ${TARGET}.log"\; exit 1\;
              fi
      DEPENDS ../src/arts ${DEPENDENCIES}
    )
    set(ADD_TARGETS ${ADD_TARGETS} ${TARGET})
  endforeach(TESTCASE)
  set(${TARGETS} ${ADD_TARGETS})
endmacro (ARTS_ADD_TESTCASES)

