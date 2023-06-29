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
  add_dependencies (pyarts_docs pyarts gen_pyarts_arts_rst gen_arts_overview_rst)

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
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_pyarts_arts_rst.py > ${CMAKE_CURRENT_BINARY_DIR}/source/pyarts_arts.rst
  )
  add_dependencies (gen_pyarts_arts_rst pyarts)
  
  add_custom_target(gen_arts_overview_rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_overview_list.py Agendas > ${CMAKE_CURRENT_BINARY_DIR}/source/auto_wsa_overview.rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_overview_list.py Variables > ${CMAKE_CURRENT_BINARY_DIR}/source/auto_wsv_overview.rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_overview_list.py Groups > ${CMAKE_CURRENT_BINARY_DIR}/source/auto_wsg_overview.rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_overview_list.py Methods > ${CMAKE_CURRENT_BINARY_DIR}/source/auto_wsm_overview.rst
  )
  add_dependencies (gen_arts_overview_rst pyarts)

  #
  # Clean
  #
  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES
    build source/pyarts_arts.rst)
endif()
