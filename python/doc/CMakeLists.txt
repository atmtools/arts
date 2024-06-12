find_package(Sphinx)

if (SPHINX_FOUND)
  set(SPHINX_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/source)
  set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/build)
  if (NOT SPHINX_JOBS)
    set(SPHINX_JOBS 1)
  endif()

  add_custom_target(pyarts_build_docs
    COMMAND
    ${SPHINX_EXECUTABLE} -j${SPHINX_JOBS} -b html
    ${SPHINX_SOURCE} ${SPHINX_BUILD}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating documentation with Sphinx")
  add_dependencies (pyarts_build_docs pyarts gen_pyarts_arts_rst gen_arts_overview_rst gen_examples_rst)

  if (PDFLATEX_COMPILER AND BIBTEX_COMPILER AND MAKEINDEX_COMPILER)
    add_custom_target(pyarts_copy_docs
      COMMAND ${CMAKE_COMMAND} -E
                      copy_if_different
                      "${ARTS_BINARY_DIR}/doc/uguide/arts_user.pdf"
                      "${ARTS_BINARY_DIR}/doc/uguide/arts_theory.pdf"
                      "${ARTS_BINARY_DIR}/doc/uguide/arts_developer.pdf"
                      "${ARTS_BINARY_DIR}/python/doc/build/_downloads/"
      COMMENT "Copying guides to documentation build directory")
    add_dependencies (pyarts_copy_docs pyarts_build_docs arts_docs)

    add_custom_target(pyarts_docs DEPENDS pyarts_copy_docs)
  else()
    add_custom_target(pyarts_docs DEPENDS pyarts_build_docs)
  endif()


  file(GLOB_RECURSE PYTHON_SOURCE_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.py")
  foreach(file ${PYTHON_SOURCE_FILES})
    configure_file(${file} ${file} @ONLY)
  endforeach()

  file(GLOB_RECURSE RST_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.rst")
  foreach(file ${RST_FILES})
    configure_file(${file} ${file} @ONLY)
  endforeach()

  file(GLOB_RECURSE MISC_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.css" "*.svg" "*.png" "*.ico")
  foreach(file ${MISC_FILES})
    configure_file(${file} ${file} COPYONLY)
  endforeach()

  configure_file(build/.nojekyll build/.nojekyll @ONLY)
  
  add_custom_target(gen_pyarts_arts_rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_pyarts_arts_rst.py > ${CMAKE_CURRENT_BINARY_DIR}/source/pyarts.arts.rst
  )
  add_dependencies (gen_pyarts_arts_rst pyarts)
  
  add_custom_target(gen_arts_overview_rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_overview_list.py Agendas > ${CMAKE_CURRENT_BINARY_DIR}/source/workspace.agendas.auto.rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_overview_list.py Variables > ${CMAKE_CURRENT_BINARY_DIR}/source/workspace.variables.auto.rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_overview_list.py Groups > ${CMAKE_CURRENT_BINARY_DIR}/source/workspace.groups.auto.rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_overview_list.py Methods > ${CMAKE_CURRENT_BINARY_DIR}/source/workspace.methods.auto.rst
  )
  add_dependencies (gen_arts_overview_rst pyarts)
  
  add_custom_target(gen_examples_rst
    COMMAND ${ARTS_PYTHON_INTERPRETER} gen_examples.py ${ARTS_SOURCE_DIR}/examples examples source
  )
  add_dependencies (gen_examples_rst pyarts)

  #
  # Clean
  #
  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES
    build source/pyarts_arts.rst)
endif()
