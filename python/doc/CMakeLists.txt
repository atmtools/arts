find_package(Sphinx)
set(SPHINX_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/source)
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/build)

add_custom_target(Sphinx
  COMMAND
  ${SPHINX_EXECUTABLE} -b html
  ${SPHINX_SOURCE} ${SPHINX_BUILD}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating documentation with Sphinx")



file(GLOB_RECURSE PYTHON_SOURCE_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.py")
foreach(file ${PYTHON_SOURCE_FILES})
  configure_file(${file} ${file} @ONLY)
endforeach()

file(GLOB_RECURSE RST_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.rst")
foreach(file ${RST_FILES})
  configure_file(${file} ${file} @ONLY)
endforeach()
