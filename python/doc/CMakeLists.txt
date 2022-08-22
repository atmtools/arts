find_package(Sphinx)

if (SPHINX_FOUND)
  set(SPHINX_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/source)
  set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/build)
  if (NOT SPHINX_JOBS)
    set(SPHINX_JOBS 1)
  endif()

  add_custom_target(pyarts_docs
    COMMAND
    ${SPHINX_EXECUTABLE} -j${SPHINX_JOBS} -b html
    ${SPHINX_SOURCE} ${SPHINX_BUILD}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating documentation with Sphinx")
  add_dependencies (pyarts_docs pyarts gen_pyarts_arts_rst)

  file(GLOB_RECURSE PYTHON_SOURCE_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.py")
  foreach(file ${PYTHON_SOURCE_FILES})
    configure_file(${file} ${file} @ONLY)
  endforeach()

  file(GLOB_RECURSE RST_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.rst")
  foreach(file ${RST_FILES})
    configure_file(${file} ${file} @ONLY)
  endforeach()

  file(GLOB_RECURSE RST_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.css")
  foreach(file ${RST_FILES})
    configure_file(${file} ${file} @ONLY)
  endforeach()

  configure_file(build/.nojekyll build/.nojekyll @ONLY)
  
  add_custom_target(gen_pyarts_arts_rst
  COMMAND PYTHONPATH=../${CMAKE_CURRENT_BINARY_DIR} ${Python3_EXECUTABLE} gen_pyarts_arts_rst.py > ${CMAKE_CURRENT_BINARY_DIR}/source/pyarts_arts.rst
  )
  add_dependencies (gen_pyarts_arts_rst pyarts)
endif()
