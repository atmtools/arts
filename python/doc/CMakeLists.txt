find_package(Sphinx)

if(SPHINX_FOUND)
  set(SPHINX_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/source)
  set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/build)

  if(NOT SPHINX_JOBS)
    set(SPHINX_JOBS 1)
  endif()

  add_custom_target(pyarts_docs
    COMMAND
    ${SPHINX_EXECUTABLE} -j${SPHINX_JOBS} -b html
    ${SPHINX_SOURCE} ${SPHINX_BUILD}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating documentation with Sphinx")
  add_dependencies(pyarts_docs
    pyarts 
    gen_pyarts_arts_rst 
    gen_arts_overview_rst
    guide_develop_doc
    copy_python_source_files_pydocs
    copy_rst_files_pydocs
    copy_misc_files_pydocs
    gen_examples_rst
  )

  file(GLOB_RECURSE PYTHON_SOURCE_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.py")
  add_custom_target(
    copy_python_source_files_pydocs
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PYTHON_SOURCE_FILES} ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${PYTHON_SOURCE_FILES}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

  file(GLOB_RECURSE RST_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.rst")
  add_custom_target(
    copy_rst_files_pydocs
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${RST_FILES} ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${RST_FILES}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

  file(GLOB_RECURSE MISC_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.css" "*.svg" "*.png" "*.ico")
  add_custom_target(
    copy_misc_files_pydocs
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${MISC_FILES} ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${MISC_FILES}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )

  configure_file(build/.nojekyll build/.nojekyll @ONLY)

  add_custom_target(gen_pyarts_arts_rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_pyarts_arts_rst.py > ${CMAKE_CURRENT_BINARY_DIR}/source/pyarts.arts.rst
  )
  add_dependencies(gen_pyarts_arts_rst pyarts)

  add_custom_target(gen_arts_overview_rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_overview_list.py Agendas > ${CMAKE_CURRENT_BINARY_DIR}/source/workspace.agendas.auto.rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_overview_list.py Variables > ${CMAKE_CURRENT_BINARY_DIR}/source/workspace.variables.auto.rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_overview_list.py Groups > ${CMAKE_CURRENT_BINARY_DIR}/source/workspace.groups.auto.rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_overview_list.py Methods > ${CMAKE_CURRENT_BINARY_DIR}/source/workspace.methods.auto.rst
  )
  add_dependencies(gen_arts_overview_rst pyarts)

  add_custom_target(gen_examples_rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_examples.py ${ARTS_SOURCE_DIR}/examples examples source
  )
  add_dependencies(gen_examples_rst pyarts)

  #
  # Clean
  #
  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES
    build source/pyarts_arts.rst)
endif()
