file(GLOB_RECURSE files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.rst")

add_custom_target(
  pyarts_docs_copy_guide
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/python/doc/source/
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${files} ${CMAKE_BINARY_DIR}/python/doc/source/
  DEPENDS ${files}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Copying RST files"
)
