find_package(Sphinx)

if(SPHINX_FOUND)
  set(SPHINX_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/source)
  set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/build)

  if(NOT SPHINX_JOBS)
    set(SPHINX_JOBS 1)
  endif()

  add_custom_target(pyarts_docs
    ${CMAKE_COMMAND} -E env PYTHONPATH="${ARTS_BINARY_DIR}/python/src" ARTS_HEADLESS=1 ARTS_CAT_DATA_DIR=${ARTS_CAT_DATA_DIR} ARTS_XML_DATA_DIR=${ARTS_XML_DATA_DIR} -- ${SPHINX_EXECUTABLE} -j${SPHINX_JOBS} -b html ${SPHINX_SOURCE} ${SPHINX_BUILD}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating documentation with Sphinx")

  add_dependencies(pyarts_docs
    pyarts 
    gen_pyarts_arts_rst 
    gen_arts_overview_rst
    gen_pyarts3_plots
    pyarts_docs_copy_guide
    pyarts_docs_copy_python
    gen_examples_rst
    check-deps
  )

  add_custom_target(
    pyarts_docs_copy_examples
    COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different ${ARTS_SOURCE_DIR}/examples/ ${ARTS_BINARY_DIR}/examples/
    COMMENT "Copying examples directory"
  )

  file(GLOB_RECURSE PYTHON_SOURCE_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/source/*")
  add_custom_target(
    pyarts_docs_copy_python
    COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different source/ ${CMAKE_CURRENT_BINARY_DIR}/source/
    DEPENDS ${PYTHON_SOURCE_FILES}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Copying python documentation structure"
  )

  file(GLOB PYTHON_DOC_GENERATORS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
  add_custom_target(
    pyarts_docs_gen_copy_python
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PYTHON_DOC_GENERATORS} ${CMAKE_CURRENT_BINARY_DIR}/
    DEPENDS ${PYTHON_DOC_GENERATORS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Copying python docmentation generation structure"
  )

  configure_file(build/.nojekyll build/.nojekyll @ONLY)

  add_custom_target(gen_pyarts_arts_rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_pyarts_arts_autoclass_rst.py source
    COMMENT "Generating ARTS interface RST files"
  )
  add_dependencies(gen_pyarts_arts_rst pyarts pyarts_docs_gen_copy_python pyarts_docs_copy_python)

  add_custom_target(gen_arts_overview_rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_overview_list.py Agendas > ${CMAKE_CURRENT_BINARY_DIR}/source/workspace.agendas.auto.rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_overview_list.py Variables > ${CMAKE_CURRENT_BINARY_DIR}/source/workspace.variables.auto.rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_overview_list.py Groups > ${CMAKE_CURRENT_BINARY_DIR}/source/workspace.groups.auto.rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_overview_list.py Methods > ${CMAKE_CURRENT_BINARY_DIR}/source/workspace.methods.auto.rst
    COMMENT "Generating overview RST files"
  )
  add_dependencies(gen_arts_overview_rst pyarts pyarts_docs_gen_copy_python pyarts_docs_copy_python)

  add_custom_target(gen_pyarts3_plots
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_pyarts3_plots.py ${CMAKE_CURRENT_BINARY_DIR}/source/
    COMMENT "Generating pyarts3 plots RST files"
  )
  add_dependencies(gen_pyarts3_plots pyarts pyarts_docs_gen_copy_python pyarts_docs_copy_python)

  if (ENABLE_EXAMPLE_PLOTS)
    set(WITH_PLOTS "plots")
  else()
    set(WITH_PLOTS "")
  endif()

  add_custom_target(gen_examples_rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_examples.py ${ARTS_BINARY_DIR} ${CMAKE_CURRENT_BINARY_DIR}/source/ ${WITH_PLOTS}
    COMMENT "Generating examples RST files"
  )
  add_dependencies(gen_examples_rst
    pyarts_docs_gen_copy_python
    pyarts_docs_copy_python
    pyarts_docs_copy_examples
  )

  #
  # Clean
  #
  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_CLEAN_FILES
    build source/pyarts_arts.rst)
endif()
