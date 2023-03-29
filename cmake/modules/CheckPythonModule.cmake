# - Check for Python modules
#
# Copyright (c) 2020, Oliver Lemke, <oliver.lemke@uni-hamburg.de>

macro (CHECK_PYTHON_MODULES)
  set(REQUIRED_MODULES
    docutils
    lark.parse_tree_builder
    matplotlib
    netCDF4
    numpy
    pytest
    scipy
    setuptools
    xarray)

  set(PYPI_NAMES
    docutils
    lark-parser
    matplotlib
    netCDF4
    numpy
    pytest
    scipy
    setuptools
    xarray)

  list(LENGTH REQUIRED_MODULES len1)
  math(EXPR len2 "${len1} - 1")

  set(PYPIERROR 0)
  foreach(i RANGE ${len2})
    list(GET REQUIRED_MODULES ${i} MODULENAME)
    list(GET PYPI_NAMES ${i} PYPINAME)
    execute_process(
      COMMAND "${Python3_EXECUTABLE}" "-c" "import ${MODULENAME}"
      RESULT_VARIABLE MODULE_FOUND
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    if(NOT MODULE_FOUND EQUAL 0)
      string(REPLACE ";" " " PYPILIST "${PYPI_NAMES}")
      message(STATUS "ERROR: Required Python package not found: ${PYPINAME}\n")
      set(PYPIERROR 1)
    endif()
  endforeach()

  if(PYPIERROR)
    message(FATAL_ERROR "Please install the missing Python package(s)")
  endif()
endmacro()
