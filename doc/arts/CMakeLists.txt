file(GLOB_RECURSE files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.rst")
file(GLOB_RECURSE bibfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.bib")

add_custom_target(
  pyarts_docs_copy_figs
  COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different Figs/ ${CMAKE_BINARY_DIR}/python/doc/source/Figs/
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Copying documentation figures"
)

add_custom_target(
  pyarts_docs_copy_bib
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/python/doc/source/
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${bibfiles} ${CMAKE_BINARY_DIR}/python/doc/source/
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Copying BIB files"
)

add_custom_target(
  pyarts_docs_copy_guide
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/python/doc/source/
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${files} ${CMAKE_BINARY_DIR}/python/doc/source/
  DEPENDS ${files} pyarts_docs_copy_figs pyarts_docs_copy_bib
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Copying RST files"
)
