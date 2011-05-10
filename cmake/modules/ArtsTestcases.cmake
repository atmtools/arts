macro (ARTS_ADD_PREPARE_TESTCASES TARGETS)
  set (ADD_TARGETS "")
  foreach(TESTCASE ${ARGN})
    set (TARGET PrepareTest${TESTCASE}Directory)
    add_custom_command(
      OUTPUT ${TARGET}
      COMMAND ${CMAKE_COMMAND} -E
        ARGS make_directory ${CMAKE_CURRENT_BINARY_DIR}/${TESTCASE}
    )
    set (ADD_TARGETS ${ADD_TARGETS} ${TARGET})

    file (GLOB COPYFILES
      "${CMAKE_CURRENT_SOURCE_DIR}/${TESTCASE}/*.arts"
      "${CMAKE_CURRENT_SOURCE_DIR}/${TESTCASE}/*.xml"
      "${CMAKE_CURRENT_SOURCE_DIR}/${TESTCASE}/*.xml.bin"
      "${CMAKE_CURRENT_SOURCE_DIR}/${TESTCASE}/*.xml.gz"
    )
    foreach (COPYFILE ${COPYFILES})
      get_filename_component (BASECOPYFILE ${COPYFILE} NAME)
      set (TARGET PrepareTest${TESTCASE}File${BASECOPYFILE})
      add_custom_command(
        OUTPUT ${TARGET}
        COMMAND ${CMAKE_COMMAND} -E
          ARGS copy_if_different ${COPYFILE}
               "${CMAKE_CURRENT_BINARY_DIR}/${TESTCASE}/"
        DEPENDS PrepareTest${TESTCASE}Directory
      )
    set (ADD_TARGETS ${ADD_TARGETS} ${TARGET})
    endforeach (COPYFILE)
  endforeach (TESTCASE)
  set(${TARGETS} ${ADD_TARGETS})
endmacro (ARTS_ADD_PREPARE_TESTCASES)

macro (ARTS_ADD_TESTCASES TARGETS DEPENDENCIES)
  set (ADD_TARGETS "")
  foreach(TESTCASE ${ARGN})
    set (TARGET Test${TESTCASE})
    add_custom_command(
      OUTPUT ${TARGET}
      COMMAND if CMAKE_GENERATOR=${CMAKE_GENERATOR} 
              python testall.py ${TARGET} > ${TARGET}.log 2>&1\; then
                  echo "${TARGET} ok."\;
	      else
                  echo "${TARGET} failed. Find details in ${TARGET}.log"\;
	          exit 1\;
              fi
      DEPENDS ../src/arts PrepareTestScripts ${DEPENDENCIES}
    )
    set(ADD_TARGETS ${ADD_TARGETS} ${TARGET})

    get_directory_property (TESTCLEANFILES ADDITIONAL_MAKE_CLEAN_FILES)
    set (TESTCLEANFILES ${TESTCLEANFILES} ${TARGET}.log)
    set_directory_properties (PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                              "${TESTCLEANFILES}")
  endforeach(TESTCASE)
  set(${TARGETS} ${ADD_TARGETS})
endmacro (ARTS_ADD_TESTCASES)

