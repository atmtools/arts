set(files
  guide.rst
  guide.theory.rst
  guide.theory.radiative_transfer.rst
  guide.user.rst
  develop.rst
  develop.classes.rst
  develop.classes.formatter.rst
)

add_custom_target(
  guide_develop_doc 
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/python/doc/source/
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${files} ${CMAKE_BINARY_DIR}/python/doc/source/
  DEPENDS ${files}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
