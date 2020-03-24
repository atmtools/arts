# - Check for Python modules
#
# Copyright (c) 2020, Oliver Lemke, <oliver.lemke@uni-hamburg.de>

macro (CHECK_PYTHON_MODULES MODULES)
  foreach(MODULENAME ${MODULES})
    execute_process(
      COMMAND "${PYTHON_EXECUTABLE}" "-c" "import ${MODULENAME}"
      RESULT_VARIABLE MODULE_FOUND
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    if(NOT MODULE_FOUND EQUAL 0)
      message(FATAL_ERROR
        "Required Python module ${MODULENAME} not found. Please install it.\n"
        "ARTS requires: ${MODULES}\n"
        "Note that the lark module is provided by the lark-parser package.")
    endif()
  endforeach()

endmacro()
