set(files
  guide.rst
  guide.theory.rst
  guide.user.rst
  develop.rst
)

foreach(fil ${files})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${fil} ${CMAKE_BINARY_DIR}/python/doc/source/${fil} COPYONLY)
endforeach(fil ${files})
